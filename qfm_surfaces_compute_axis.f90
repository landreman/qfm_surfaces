subroutine qfm_surfaces_compute_axis

  use qfm_surfaces_variables

  implicit none

  integer :: iteration, j_line_search, vector_size, n, j
  real(dp) :: residual_norm, last_residual_norm, initial_residual_norm, step_scale
  real(dp), dimension(:), allocatable :: state, state0, residual, residual0, step_direction, R_residual, Z_residual

  real(dp), dimension(:), allocatable :: phi_axis, phi_weights, R_axis, Z_axis
  real(dp), dimension(:), allocatable :: B_R, B_phi, B_Z
  real(dp), dimension(:), allocatable :: R_axis_reconstruction, Z_axis_reconstruction
  real(dp), dimension(:,:), allocatable :: differentiation_matrix, DD, Jacobian
  real(dp) :: factor, sinangle, cosangle

  ! Variables needed by LAPACK:                                                                                            
  integer :: INFO
  integer, dimension(:), allocatable :: IPIV

  if (proc0) print "(a)"," Beginning solve for magnetic axis."

  ! For N_phi_axis to be odd.
  if (mod(N_phi_axis,2)==0) N_phi_axis = N_phi_axis + 1

  allocate(phi_axis(N_phi_axis))
  allocate(R_axis(N_phi_axis))
  allocate(Z_axis(N_phi_axis))
  allocate(R_residual(N_phi_axis))
  allocate(Z_residual(N_phi_axis))
  allocate(phi_weights(N_phi_axis))
  allocate(differentiation_matrix(N_phi_axis,N_phi_axis))
  allocate(DD(N_phi_axis,N_phi_axis))
  allocate(B_R(N_phi_axis))
  allocate(B_phi(N_phi_axis))
  allocate(B_Z(N_phi_axis))
  call qfm_surfaces_differentiation_matrix(N_phi_axis, 0.0d+0, 2*pi/nfp, 20, 1, phi_axis, phi_weights, differentiation_matrix, DD)
  deallocate(phi_weights,DD)

  if (stellarator_symmetry) then
     vector_size = N_phi_axis
  else
     vector_size = N_phi_axis * 2
  end if

  allocate(state(vector_size))
  allocate(state0(vector_size))
  allocate(IPIV(vector_size))
  allocate(residual(vector_size))
  allocate(residual0(vector_size))
  allocate(step_direction(vector_size))
  allocate(Jacobian(vector_size,vector_size))

  ! Initialize state
  state = 0
  if (stellarator_symmetry) then
     state(1:((N_phi_axis+1)/2)) = 1 ! Set R to a constant of 1
  else
     state(1:N_phi_axis) = 1
  end if

  call qfm_surfaces_axis_residual()
  initial_residual_norm = sqrt(sum(residual * residual))
  residual_norm = initial_residual_norm
  if (proc0) print "(a,es10.3)","                 Initial residual L2 norm:",residual_norm

  ! Here is the main Newton iteration:
  Newton: do iteration = 1, N_iterations_axis
     last_residual_norm = residual_norm
     !if (residual_norm / initial_residual_norm < Newton_tolerance) then
     if (residual_norm < Newton_tolerance_axis) then
        exit Newton
     end if

     call qfm_surfaces_axis_Jacobian()

     state0 = state
     if (proc0) print "(a,i3)","  Newton iteration ",iteration
     ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
     ! step_direction = - matrix \ residual
     step_direction = -residual ! Note that LAPACK will over-write step_direction and with the solution, and over-write Jacobian with the LU factorization.
     call DGESV(vector_size, 1, Jacobian, vector_size, IPIV, step_direction, vector_size, INFO)
     if (INFO /= 0) then
        print *, "Error in LAPACK call DGESV: info = ", INFO
        stop
     end if

     step_scale = 1
     line_search: do j_line_search = 1, N_line_search_axis
        state = state0 + step_scale * step_direction

        call qfm_surfaces_axis_residual()
        residual_norm = sqrt(sum(residual * residual))
        !if (verbose) print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Relative residual L2 norm:",residual_norm / initial_residual_norm,"  iota:",iota
        if (proc0) print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Residual L2 norm:",residual_norm
        if (residual_norm < last_residual_norm) exit line_search

        step_scale = step_scale / 2
     end do line_search

     if (residual_norm > last_residual_norm) then
        if (proc0) print *,"Line search failed to reduce residual."
        exit Newton
     end if
  end do Newton
  ! End of Newton solve.

  !print *,"Final state:",state

  call qfm_surfaces_unpack_axis_state()
  !print *,"Final R_axis:",R_axis
  !print *,"Final Z_axis:",Z_axis

  ! ------------------------------------------------
  ! Now that we have the axis shape on a grid in phi,
  ! Fourier transform it.

  nmax_axis = (N_phi_axis-1)/2
  allocate(R0c(nmax_axis+1))
  allocate(R0s(nmax_axis+1))
  allocate(Z0c(nmax_axis+1))
  allocate(Z0s(nmax_axis+1))
  R0c = 0
  R0s = 0
  Z0c = 0
  Z0s = 0
  R0c(1) = sum(R_axis) / N_phi_axis
  Z0c(1) = sum(Z_axis) / N_phi_axis
  factor = 2.0d+0 / N_phi_axis
  do n = 1, nmax_axis
     do j = 1, N_phi_axis
        sinangle = sin(n * nfp * phi_axis(j))
        cosangle = cos(n * nfp * phi_axis(j))
        R0c(n+1) = R0c(n+1) + R_axis(j) * cosangle * factor
        R0s(n+1) = R0s(n+1) + R_axis(j) * sinangle * factor
        Z0c(n+1) = Z0c(n+1) + Z_axis(j) * cosangle * factor
        Z0s(n+1) = Z0s(n+1) + Z_axis(j) * sinangle * factor
     end do
  end do
  if (.not. stellarator_symmetry) then
     R0s = 0
     Z0c = 0
  end if

  allocate(R_axis_reconstruction(N_phi_axis))
  allocate(Z_axis_reconstruction(N_phi_axis))
  R_axis_reconstruction = 0
  Z_axis_reconstruction = 0
  do n = 0, nmax_axis
     do j = 1, N_phi_axis
        sinangle = sin(n * nfp * phi_axis(j))
        cosangle = cos(n * nfp * phi_axis(j))
        R_axis_reconstruction(j) = R_axis_reconstruction(j) + R0c(n+1) * cosangle + R0s(n+1) * sinangle
        Z_axis_reconstruction(j) = Z_axis_reconstruction(j) + Z0c(n+1) * cosangle + Z0s(n+1) * sinangle
     end do
  end do

  factor = maxval(abs(R_axis - R_axis_reconstruction))
  if (proc0) print *,"Residual from in R transform:",factor
  if (factor > 1.0d-12) then
     print *,"R_axis before transform:",R_axis
     print *,"R_axis after transform: ",R_axis_reconstruction
     stop "Residual is too large in R transform"
  end if

  factor = maxval(abs(Z_axis - Z_axis_reconstruction))
  if (proc0) print *,"Residual from in Z transform:",factor
  if (factor > 1.0d-12) then
     print *,"Z_axis before transform:",Z_axis
     print *,"Z_axis after transform: ",Z_axis_reconstruction
     stop "Zesidual is too large in Z transform"
  end if

  deallocate(state, state0, IPIV, step_direction)
  deallocate(residual, residual0, Jacobian, R_axis, Z_axis, R_residual, Z_residual)
  deallocate(B_R, B_phi, B_Z)
  deallocate(phi_axis, differentiation_matrix)

  if (proc0) print "(a)"," Done solving for magnetic axis."
contains

! -----------------------------------------------------------------------------------

  subroutine qfm_surfaces_unpack_axis_state()
    ! Unpack state vector into R_axis and Z_axis

    implicit none

    if (stellarator_symmetry) then
       R_axis(1:((N_phi_axis+1)/2)) = state(1:((N_phi_axis+1)/2))
       !print *,"before:",R_axis
       R_axis(N_phi_axis:((N_phi_axis+3)/2):-1) = state(2:((N_phi_axis+1)/2))
       !print *,"after: ",R_axis
       Z_axis(1) = 0
       Z_axis(2:((N_phi_axis+1)/2)) = state(((N_phi_axis+3)/2):N_phi_axis)
       Z_axis(N_phi_axis:((N_phi_axis+3)/2):-1) = -state(((N_phi_axis+3)/2):N_phi_axis)
    else
       R_axis = state(1:N_phi_axis)
       Z_axis = state(N_phi_axis+1:vector_size)
    end if

  end subroutine qfm_surfaces_unpack_axis_state

  ! -----------------------------------------------------------------------------------

  subroutine qfm_surfaces_axis_residual()

    implicit none

    call qfm_surfaces_unpack_axis_state()

!!$    print *,"state:",state
!!$    print *,"R_axis:",R_axis
!!$    print *,"Z_axis:",Z_axis

    call qfm_surfaces_compute_B(N_phi_axis, R_axis, phi_axis, Z_axis, B_R, B_phi, B_Z)

    R_residual = matmul(differentiation_matrix, R_axis) / R_axis - B_R / B_phi
    Z_residual = matmul(differentiation_matrix, Z_axis) / R_axis - B_Z / B_phi

!!$    print *,"matmul(differentiation_matrix, R_axis):",matmul(differentiation_matrix, R_axis)
!!$    print *,"B_R:",B_R
!!$    print *,"B_phi:",B_phi
!!$    print *,"B_Z:",B_Z
!!$    print *,"R_residual:",R_residual
!!$    print *,"Z_residual:",Z_residual

    if (stellarator_symmetry) then
       residual(1:((N_phi_axis+1)/2)) = Z_residual(1:((N_phi_axis+1)/2))
       residual(((N_phi_axis+3)/2):N_phi_axis) = R_residual(((N_phi_axis+3)/2):N_phi_axis)
    else
       residual(1:N_phi_axis) = R_residual
       residual(N_phi_axis+1:vector_size) = Z_residual
    end if

!    print *,"Residual:",residual

  end subroutine qfm_surfaces_axis_residual

  ! -----------------------------------------------------------------------------------

  subroutine qfm_surfaces_axis_Jacobian()

    implicit none

    real(dp) :: epsilon = 1.0d-7
    integer :: j

    residual0 = residual
    state0 = state

    do j = 1, vector_size
       state = state0
       state(j) = state(j) + epsilon
       call qfm_surfaces_axis_residual
       Jacobian(:,j) = (residual - residual0) / epsilon
    end do

    state = state0
    residual = residual0

  end subroutine qfm_surfaces_axis_Jacobian

end subroutine qfm_surfaces_compute_axis
