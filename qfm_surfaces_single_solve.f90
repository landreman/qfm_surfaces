subroutine qfm_surfaces_single_solve(j_volume)

  use qfm_surfaces_variables

  implicit none

  integer, intent(in) :: j_volume
  integer :: j_resolution, m, n, j, k
  integer :: N_theta, N_phi, mpol, ntor
  real(dp), dimension(:), allocatable :: theta, phi, temp_1D
  real(dp), dimension(:,:), allocatable :: ddtheta, ddphi, temp_2D, Jacobian
  real(dp), dimension(:), allocatable :: R0, d_R0_d_phi, d2_R0_d_phi2
  real(dp), dimension(:), allocatable :: Z0, d_Z0_d_phi, d2_Z0_d_phi2
  real(dp), dimension(:), allocatable :: state_vector, state_vector0, last_state_vector, residual, residual0, step_direction
  real(dp) :: angle, sinangle, cosangle, volume_target, dtheta, dphi
  integer :: mnmax, vector_size, last_vector_size=0, index, jm, jn
  integer, dimension(:), allocatable :: xm, xn, last_xm, last_xn
  logical :: found_match
  real(dp), dimension(:,:), allocatable :: sin_m_theta, cos_m_theta, sin_n_phi, cos_n_phi
  real(dp), dimension(:), allocatable :: sin_phi, cos_phi, sin_theta, cos_theta
  real(dp), dimension(:,:), allocatable :: Z, R, dXdtheta, dYdtheta, dZdtheta, dXdphi, dYdphi, dZdphi
  real(dp), dimension(:,:), allocatable :: d2Xdtheta2, d2Ydtheta2, d2Zdtheta2, d2Xdphi2, d2Ydphi2, d2Zdphi2
  real(dp), dimension(:,:), allocatable :: d2Xdthetadphi, d2Ydthetadphi, d2Zdthetadphi
  real(dp), dimension(:,:), allocatable :: X_base, Y_base, Z_base, R_base, dXdtheta_base, dYdtheta_base, dZdtheta_base, dXdphi_base, dYdphi_base, dZdphi_base
  real(dp), dimension(:,:), allocatable :: d2Xdtheta2_base, d2Ydtheta2_base, d2Zdtheta2_base, d2Xdphi2_base, d2Ydphi2_base, d2Zdphi2_base
  real(dp), dimension(:,:), allocatable :: d2Xdthetadphi_base, d2Ydthetadphi_base, d2Zdthetadphi_base
  real(dp), parameter :: epsilon = 1.0d-7
  integer, parameter :: normal_vector_sign = -1
  real(dp), dimension(:,:), allocatable :: B_R, B_phi, B_Z, B_X, B_Y, NX, NY, NZ, NR, normal_X, normal_Y, normal_Z, Bnormal, norm_normal
  real(dp), dimension(:,:), allocatable :: fundamental_form_E, fundamental_form_F, fundamental_form_G, fundamental_form_L, fundamental_form_M, fundamental_form_N
  real(dp), dimension(:,:), allocatable :: mean_curvature, shape_gradient, B_dot_e_theta, B_dot_e_phi, N_dot_B_cross_e_phi, N_dot_e_theta_cross_B
  real(dp), dimension(:,:), allocatable :: d_Bnormal_d_theta, d_Bnormal_d_phi, integrand
  real(dp), dimension(:), allocatable :: phi_copied
  integer :: iteration, j_line_search
  real(dp) :: residual_norm, initial_residual_norm, last_residual_norm, step_scale
  logical :: verbose = .true.

  ! Variables needed by LAPACK:                                                                                            
  integer :: INFO
  integer, dimension(:), allocatable :: IPIV

  volume_target = volumes(j_volume)
  print "(a,i5,a,i5,a)"," Solving for volume",j_volume," of",N_volumes,"."
  print *,"volume_target=",volume_target

  do j_resolution = 1, N_resolutions
     mpol = mpols(j_resolution)
     ntor = ntors(j_resolution)
     N_theta = max(min_N_theta,mpol*4)
     N_phi = max(min_N_phi,ntor*4)
     print "(a,i5,a,i5,a,i6,a,i6,a)","  Using resolution mpol=",mpol,", ntor=",ntor,", N_theta=",N_theta,", N_phi=",N_phi,"."

     allocate(theta(N_theta))
     allocate(temp_1D(N_theta))
     allocate(ddtheta(N_theta,N_theta))
     allocate(temp_2D(N_theta,N_theta))
     call qfm_surfaces_differentiation_matrix(N_theta, 0.0d+0, 2*pi, 20, 1, theta, temp_1D, ddtheta, temp_2D)
     deallocate(temp_1D)
     deallocate(temp_2D)
     dtheta = theta(2) - theta(1)

     allocate(phi(N_phi))
     allocate(temp_1D(N_phi))
     allocate(ddphi(N_phi,N_phi))
     allocate(temp_2D(N_phi,N_phi))
     call qfm_surfaces_differentiation_matrix(N_phi, 0.0d+0, 2*pi, 20, 1, phi, temp_1D, ddphi, temp_2D)
     deallocate(temp_1D)
     deallocate(temp_2D)
     dphi = phi(2) - phi(1)

     allocate(sin_m_theta(N_theta,0:mpol))
     allocate(cos_m_theta(N_theta,0:mpol))
     allocate(sin_n_phi(N_phi,-ntor:ntor))
     allocate(cos_n_phi(N_phi,-ntor:ntor))
     do m = 0, mpol
        do j = 1, N_theta
           sin_m_theta(j,m) = sin(m*theta(j))
           cos_m_theta(j,m) = cos(m*theta(j))
        end do
     end do
     do n = -ntor, ntor
        do j = 1, N_phi
           sin_n_phi(j,n) = sin(n*nfp*phi(j))
           cos_n_phi(j,n) = cos(n*nfp*phi(j))
        end do
     end do

     allocate(sin_phi(N_phi))
     allocate(cos_phi(N_phi))
     do j = 1, N_phi
        sin_phi(j) = sin(phi(j))
        cos_phi(j) = cos(phi(j))
     end do

     allocate(sin_theta(N_theta))
     allocate(cos_theta(N_theta))
     do j = 1, N_theta
        sin_theta(j) = sin(theta(j))
        cos_theta(j) = cos(theta(j))
     end do

     !allocate(X(N_theta,N_phi))
     !allocate(Y(N_theta,N_phi))
     allocate(Z(N_theta,N_phi))
     allocate(R(N_theta,N_phi))
     allocate(dXdtheta(N_theta,N_phi))
     allocate(dYdtheta(N_theta,N_phi))
     allocate(dZdtheta(N_theta,N_phi))
     allocate(dXdphi(N_theta,N_phi))
     allocate(dYdphi(N_theta,N_phi))
     allocate(dZdphi(N_theta,N_phi))
     allocate(d2Xdtheta2(N_theta,N_phi))
     allocate(d2Ydtheta2(N_theta,N_phi))
     allocate(d2Zdtheta2(N_theta,N_phi))
     allocate(d2Xdthetadphi(N_theta,N_phi))
     allocate(d2Ydthetadphi(N_theta,N_phi))
     allocate(d2Zdthetadphi(N_theta,N_phi))
     allocate(d2Xdphi2(N_theta,N_phi))
     allocate(d2Ydphi2(N_theta,N_phi))
     allocate(d2Zdphi2(N_theta,N_phi))

     !allocate(X_base(N_theta,N_phi))
     !allocate(Y_base(N_theta,N_phi))
     allocate(Z_base(N_theta,N_phi))
     allocate(R_base(N_theta,N_phi))
     allocate(dXdtheta_base(N_theta,N_phi))
     allocate(dYdtheta_base(N_theta,N_phi))
     allocate(dZdtheta_base(N_theta,N_phi))
     allocate(dXdphi_base(N_theta,N_phi))
     allocate(dYdphi_base(N_theta,N_phi))
     allocate(dZdphi_base(N_theta,N_phi))
     allocate(d2Xdtheta2_base(N_theta,N_phi))
     allocate(d2Ydtheta2_base(N_theta,N_phi))
     allocate(d2Zdtheta2_base(N_theta,N_phi))
     allocate(d2Xdthetadphi_base(N_theta,N_phi))
     allocate(d2Ydthetadphi_base(N_theta,N_phi))
     allocate(d2Zdthetadphi_base(N_theta,N_phi))
     allocate(d2Xdphi2_base(N_theta,N_phi))
     allocate(d2Ydphi2_base(N_theta,N_phi))
     allocate(d2Zdphi2_base(N_theta,N_phi))

     allocate(B_R(N_theta,N_phi))
     allocate(B_phi(N_theta,N_phi))
     allocate(B_Z(N_theta,N_phi))
     allocate(B_X(N_theta,N_phi))
     allocate(B_Y(N_theta,N_phi))
     allocate(NX(N_theta,N_phi))
     allocate(NY(N_theta,N_phi))
     allocate(NZ(N_theta,N_phi))
     allocate(NR(N_theta,N_phi))
     allocate(normal_X(N_theta,N_phi))
     allocate(normal_Y(N_theta,N_phi))
     allocate(normal_Z(N_theta,N_phi))
     allocate(Bnormal(N_theta,N_phi))
     allocate(norm_normal(N_theta,N_phi))
     allocate(fundamental_form_E(N_theta,N_phi))
     allocate(fundamental_form_F(N_theta,N_phi))
     allocate(fundamental_form_G(N_theta,N_phi))
     allocate(fundamental_form_L(N_theta,N_phi))
     allocate(fundamental_form_M(N_theta,N_phi))
     allocate(fundamental_form_N(N_theta,N_phi))
     allocate(mean_curvature(N_theta,N_phi))
     allocate(shape_gradient(N_theta,N_phi))
     allocate(B_dot_e_theta(N_theta,N_phi))
     allocate(B_dot_e_phi(N_theta,N_phi))
     allocate(N_dot_B_cross_e_phi(N_theta,N_phi))
     allocate(N_dot_e_theta_cross_B(N_theta,N_phi))
     allocate(d_Bnormal_d_theta(N_theta,N_phi))
     allocate(d_Bnormal_d_phi(N_theta,N_phi))
     allocate(integrand(N_theta,N_phi))
     allocate(phi_copied(N_theta))


     ! Evaluate the contribution to the surface shape from the axis, and its derivatives:
     allocate(R0(N_phi))
     allocate(d_R0_d_phi(N_phi))
     allocate(d2_R0_d_phi2(N_phi))
     allocate(Z0(N_phi))
     allocate(d_Z0_d_phi(N_phi))
     allocate(d2_Z0_d_phi2(N_phi))
     R0 = 0
     d_R0_d_phi = 0
     d2_R0_d_phi2 = 0
     Z0 = 0
     d_Z0_d_phi = 0
     d2_Z0_d_phi2 = 0
     do n = 0, nmax_axis
        do j = 1, N_phi
           angle = n * nfp * phi(j)
           sinangle = sin(angle)
           cosangle = cos(angle)
           R0(j) = R0(j) + R0c(n+1) * cosangle + R0s(n+1) * sinangle
           d_R0_d_phi(j) = d_R0_d_phi(j) + n * nfp * (R0c(n+1) * (-sinangle) + R0s(n+1) * cosangle)
           d2_R0_d_phi2(j) = d2_R0_d_phi2(j) + n * n * nfp * nfp * (R0c(n+1) * (-cosangle) + R0s(n+1) * (-sinangle))
           Z0(j) = Z0(j) + Z0c(n+1) * cosangle + Z0s(n+1) * sinangle
           d_Z0_d_phi(j) = d_Z0_d_phi(j) + n * nfp * (Z0c(n+1) * (-sinangle) + Z0s(n+1) * cosangle)
           d2_Z0_d_phi2(j) = d2_Z0_d_phi2(j) + n * n * nfp * nfp * (Z0c(n+1) * (-cosangle) + Z0s(n+1) * (-sinangle))
        end do
     end do

     ! Initialize mnmax, xm, and xn:
     mnmax = (ntor*2 + 1) * mpol + ntor + 1
     vector_size = mnmax + 1 ! +1 for volume constraint.
     allocate(xm(mnmax))
     allocate(xn(mnmax))
     xm = 0
     xn = 0
     xn(2:(ntor+1)) = [( j, j=1, ntor )]
     index = ntor+1
     do jm = 1,mpol
        do jn = -ntor, ntor
           index = index + 1
           xn(index) = jn
           xm(index) = jm
        end do
     end do
     !xn = xn * nfp
     !print "(a,*(i3))","xm:",xm
     !print "(a,*(i3))","xn:",xn

     ! Initialize state vector
     allocate(state_vector(vector_size))
     allocate(state_vector0(vector_size))
     allocate(residual(vector_size))
     allocate(residual0(vector_size))
     allocate(Jacobian(vector_size,vector_size))
     allocate(step_direction(vector_size))
     state_vector = 0
     if (j_resolution==1) then
        state_vector(1) = sqrt(volume_target/(2*pi*(sum(R0)/N_phi)*pi))
     else
        ! Copy previous state vector
!!$        print *,"Copying previous state vector."
!!$        print *,"size(last_xm):",size(last_xm)
!!$        print *,"size(last_xn):",size(last_xn)
!!$        print *,"size(xm):",size(xm)
!!$        print *,"size(xn):",size(xn)
!!$        print *,"last_xm:",last_xm
!!$        print *,"last_xn:",last_xn
!!$        print *,"xm:",xm
!!$        print *,"xn:",xn
        do j = 1, (last_vector_size-1)
           !print *,"j=",j
           found_match = .false.
           do k = 1, mnmax
              !print *,"   k=",k
              if (last_xm(j)==xm(k) .and. last_xn(j)==xn(k)) then
                 found_match = .true.
                 state_vector(k) = last_state_vector(j)
                 cycle
              end if
           end do
           if (.not. found_match) then
              print *,"Error copying state vector!"
              stop
           end if
        end do
        state_vector(vector_size) = last_state_vector(last_vector_size)
        print *,"last_xm:",last_xm
        print *,"     xm:",xm
        print *,"last_xn:",last_xn
        print *,"     xn:",xn
        print *,"last_state_vector:",last_state_vector
        print *,"     state_vector:",state_vector
        
        deallocate(last_xm, last_xn, last_state_vector)
     end if

     ! As a prelude to the Newton iteration, get the initial residual:
     call qfm_surfaces_residual(0)
     initial_residual_norm = sqrt(sum(residual * residual))
     residual_norm = initial_residual_norm
     print "(a,es10.3)","                 Initial residual L2 norm:",residual_norm
     print *,"residual:",residual
     return

     ! Here is the main Newton iteration:
     Newton: do iteration = 1, N_iterations
        last_residual_norm = residual_norm
        !if (residual_norm / initial_residual_norm < Newton_tolerance) then
        if (residual_norm < Newton_tolerance) then
           exit Newton
        end if
        
        call qfm_surfaces_Jacobian()
        
        state_vector0 = state_vector
        if (verbose) print "(a,i3)","  Newton iteration ",iteration
        ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
        ! step_direction = - matrix \ residual
        step_direction = -residual ! Note that LAPACK will over-write step_direction and with the solution, and over-write Jacobian with the LU factorization.
        call DGESV(vector_size, 1, Jacobian, vector_size, IPIV, step_direction, vector_size, INFO)
        if (INFO /= 0) then
           print *, "Error in LAPACK call DGESV: info = ", INFO
           stop
        end if
        
        step_scale = 1
        line_search: do j_line_search = 1, N_line_search
           state_vector = state_vector0 + step_scale * step_direction
           
           call qfm_surfaces_residual(0)
           residual_norm = sqrt(sum(residual * residual))
           !if (verbose) print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Relative residual L2 norm:",residual_norm / initial_residual_norm,"  iota:",iota
           if (verbose) print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Residual L2 norm:",residual_norm
           if (residual_norm < last_residual_norm) exit line_search
           
           step_scale = step_scale / 2
        end do line_search
        
        if (residual_norm > last_residual_norm) then
           print *,"Line search failed to reduce residual."
           exit Newton
        end if
     end do Newton
     ! End of Newton solve.
     
     !print *,"Final state_vector:",state_vector


     allocate(last_state_vector(vector_size))
     allocate(last_xm(mnmax))
     allocate(last_xn(mnmax))
     last_vector_size = vector_size
     last_state_vector = state_vector
     last_xm = xm
     last_xn = xn

     ! Done with this resolution. Deallocate everything.
     deallocate(phi, ddphi, theta, ddtheta)
     deallocate(R0, d_R0_d_phi, d2_R0_d_phi2, Z0, d_Z0_d_phi, d2_Z0_d_phi2)
     deallocate(xm, xn, state_vector, state_vector0, residual, residual0, Jacobian, step_direction)
     deallocate(sin_m_theta, cos_m_theta, sin_n_phi, cos_n_phi, sin_phi, cos_phi, sin_theta, cos_theta)
     deallocate(Z,R)
     deallocate(dXdtheta,dYdtheta,dZdtheta,dXdphi,dYdphi,dZdphi)
     deallocate(d2Xdtheta2,d2Ydtheta2,d2Zdtheta2,d2Xdphi2,d2Ydphi2,d2Zdphi2)
     deallocate(d2Xdthetadphi,d2Ydthetadphi,d2Zdthetadphi)
     deallocate(X_base,Y_base,Z_base,R_base)
     deallocate(dXdtheta_base,dYdtheta_base,dZdtheta_base,dXdphi_base,dYdphi_base,dZdphi_base)
     deallocate(d2Xdtheta2_base,d2Ydtheta2_base,d2Zdtheta2_base,d2Xdphi2_base,d2Ydphi2_base,d2Zdphi2_base)
     deallocate(d2Xdthetadphi_base,d2Ydthetadphi_base,d2Zdthetadphi_base)
     deallocate(B_R, B_phi, B_Z, B_X, B_Y, NX, NY, NZ, NR, normal_X, normal_Y, normal_Z, Bnormal, norm_normal)
     deallocate(fundamental_form_E, fundamental_form_F, fundamental_form_G, fundamental_form_L, fundamental_form_M, fundamental_form_N)
     deallocate(mean_curvature, shape_gradient, B_dot_e_theta, B_dot_e_phi, N_dot_B_cross_e_phi, N_dot_e_theta_cross_B)
     deallocate(d_Bnormal_d_theta, d_Bnormal_d_phi, integrand, phi_copied)

  end do ! Loop over resolution
  deallocate(last_xm, last_xn, last_state_vector)


contains

  ! -----------------------------------------------------------------------------------

  subroutine qfm_surfaces_Jacobian()

    implicit none

    integer :: j

    residual0 = residual
    state_vector0 = state_vector

    do j = 1, vector_size
       if (j==vector_size) state_vector(vector_size) = state_vector(vector_size) + epsilon
       ! Aside from the last column, the perturbations to the state vector are handled inside qfm_surfaces_residual for the sake of speed.
       call qfm_surfaces_residual(j)
       Jacobian(:,j) = (residual - residual0) / epsilon
    end do

    residual = residual0
    state_vector = state_vector0

  end subroutine qfm_surfaces_Jacobian

  ! -----------------------------------------------------------------------------------

  subroutine qfm_surfaces_residual(Jacobian_column)
    
    implicit none

    integer, intent(in) :: Jacobian_column
    integer :: mnmax_to_use, imn, m, n
    logical :: computing_Jacobian
    integer :: itheta, iphi
    real(dp) :: amnc, d_cosangle_dtheta, d_cosangle_dphi, d2_cosangle_dtheta2, d2_cosangle_dphi2, d2_cosangle_dtheta_dphi

    print "(a,i3)","Entering qfm_surfaces_residual. Jacobian_column=",Jacobian_column
    print *,"state_vector:",state_vector
    computing_Jacobian = (Jacobian_column > 0)

    if (computing_Jacobian) then
       ! Initialize arrays to the base case
       !X = X_base
       !Y = Y_base
       Z = Z_base
       R = R_base
       dXdtheta = dXdtheta_base
       dYdtheta = dYdtheta_base
       dZdtheta = dZdtheta_base
       dXdphi = dXdphi_base
       dYdphi = dYdphi_base
       dZdphi = dZdphi_base
       d2Xdtheta2 = d2Xdtheta2_base
       d2Ydtheta2 = d2Ydtheta2_base
       d2Zdtheta2 = d2Zdtheta2_base
       d2Xdthetadphi = d2Xdthetadphi_base
       d2Ydthetadphi = d2Ydthetadphi_base
       d2Zdthetadphi = d2Zdthetadphi_base
       d2Xdphi2 = d2Xdphi2_base
       d2Ydphi2 = d2Ydphi2_base
       d2Zdphi2 = d2Zdphi2_base
       mnmax_to_use = 1
       if (Jacobian_column == vector_size) mnmax_to_use = 0 ! If perturbing the Lagrange multiplier lambda
    else
       ! We are not computing the finite-difference Jacobian
       !X = 0
       !Y = 0
       Z = 0
       R = 0
       dXdtheta = 0
       dYdtheta = 0
       dZdtheta = 0
       dXdphi = 0
       dYdphi = 0
       dZdphi = 0
       d2Xdtheta2 = 0
       d2Ydtheta2 = 0
       d2Zdtheta2 = 0
       d2Xdthetadphi = 0
       d2Ydthetadphi = 0
       d2Zdthetadphi = 0
       d2Xdphi2 = 0
       d2Ydphi2 = 0
       d2Zdphi2 = 0
       ! Add contribution from the axis
       do itheta = 1, N_theta
          !X(itheta,:) = R0 * cos_phi
          !Y(itheta,:) = R0 * sin_phi
          Z(itheta,:) = Z0
          R(itheta,:) = R0
          dXdphi(itheta,:) = d_R0_d_phi * cos_phi + R0 * (-sin_phi)
          dYdphi(itheta,:) = d_R0_d_phi * sin_phi + R0 * cos_phi
          dXdphi(itheta,:) = d_Z0_d_phi
          d2Xdphi2(itheta,:) = d2_R0_d_phi2 * cos_phi + 2 * d_R0_d_phi * (-sin_phi) + R0 * (-cos_phi)
          d2Ydphi2(itheta,:) = d2_R0_d_phi2 * sin_phi + 2 * d_R0_d_phi * cos_phi + R0 * (-sin_phi)
       end do
       mnmax_to_use = mnmax
    end if
    
    print *,"R0:",R0
    print *,"Z0:",Z0
    print *,"cos_theta:",cos_theta

    do imn = 1, mnmax_to_use
       if (computing_Jacobian) then
          m = xm(Jacobian_column)
          n = xn(Jacobian_column)
          amnc = epsilon
       else
          m = xm(imn)
          n = xn(imn)
          amnc = state_vector(imn)
       end if
       print *,"imn=",imn,"amnc=",amnc
       do iphi = 1, N_phi
          do itheta = 1, N_theta
             !sinangle = sin(m*theta-n*phi) = sin(m*theta) * cos(n*phi) - cos(m*theta) * sin(n*phi)
             sinangle = sin_m_theta(itheta,m) * cos_n_phi(iphi,n) - cos_m_theta(itheta,m) * sin_n_phi(iphi,n)
             !cosangle = cos(m*theta-n*phi) = cos(m*theta) * cos(n*phi) + sin(m*theta) * sin(n*phi)
             cosangle = cos_m_theta(itheta,m) * cos_n_phi(iphi,n) + sin_m_theta(itheta,m) * sin_n_phi(iphi,n)
             print *,"cosangle:",cosangle

             d_cosangle_dtheta = -m*sinangle
             d_cosangle_dphi  =  n*nfp*sinangle
             d2_cosangle_dtheta2 = -m*m*cosangle
             d2_cosangle_dphi2  = -n*n*nfp*nfp*cosangle
             d2_cosangle_dtheta_dphi = m*n*nfp*cosangle
            
             R(itheta,iphi) = R(itheta,iphi) + amnc * cosangle * cos_theta(itheta)
             !X(itheta,iphi) = X(itheta,iphi) + amnc * cosangle * cos_theta(itheta) * cos_phi(iphi)
             !Y(itheta,iphi) = Y(itheta,iphi) + amnc * cosangle * cos_theta(itheta) * sin_phi(iphi)
             Z(itheta,iphi) = Z(itheta,iphi) + amnc * cosangle * sin_theta(itheta)
            
             dXdtheta(itheta,iphi) = dXdtheta(itheta,iphi) + amnc * (cosangle * (-sin_theta(itheta)) + d_cosangle_dtheta * cos_theta(itheta)) * cos_phi(iphi)
             dYdtheta(itheta,iphi) = dYdtheta(itheta,iphi) + amnc * (cosangle * (-sin_theta(itheta)) + d_cosangle_dtheta * cos_theta(itheta)) * sin_phi(iphi)
             dZdtheta(itheta,iphi) = dZdtheta(itheta,iphi) + amnc * (cosangle *   cos_theta(itheta)  + d_cosangle_dtheta * sin_theta(itheta))
            
             dXdphi(itheta,iphi) = dXdphi(itheta,iphi) + amnc * (cosangle * (-sin_phi(iphi)) + d_cosangle_dphi * cos_phi(iphi)) * cos_theta(itheta)
             dYdphi(itheta,iphi) = dYdphi(itheta,iphi) + amnc * (cosangle *   cos_phi(iphi)  + d_cosangle_dphi * sin_phi(iphi)) * cos_theta(itheta)
             dZdphi(itheta,iphi) = dZdphi(itheta,iphi) + amnc * d_cosangle_dphi * sin_theta(itheta)
            
             d2Xdtheta2(itheta,iphi) = d2Xdtheta2(itheta,iphi) + amnc * (cosangle * (-cos_theta(itheta)) + 2*d_cosangle_dtheta * (-sin_theta(itheta)) + d2_cosangle_dtheta2 * cos_theta(itheta)) * cos_phi(iphi)
             d2Ydtheta2(itheta,iphi) = d2Ydtheta2(itheta,iphi) + amnc * (cosangle * (-cos_theta(itheta)) + 2*d_cosangle_dtheta * (-sin_theta(itheta)) + d2_cosangle_dtheta2 * cos_theta(itheta)) * sin_phi(iphi)
             d2Zdtheta2(itheta,iphi) = d2Zdtheta2(itheta,iphi) + amnc * (cosangle * (-sin_theta(itheta)) + 2*d_cosangle_dtheta *   cos_theta(itheta)  + d2_cosangle_dtheta2 * sin_theta(itheta))
             
             d2Xdthetadphi(itheta,iphi) = d2Xdthetadphi(itheta,iphi) + amnc * ((cosangle * (-sin_theta(itheta)) + d_cosangle_dtheta * cos_theta(itheta)) * (-sin_phi(iphi)) + (d_cosangle_dphi * (-sin_theta(itheta)) + d2_cosangle_dtheta_dphi * cos_theta(itheta)) * cos_phi(iphi))
             d2Ydthetadphi(itheta,iphi) = d2Ydthetadphi(itheta,iphi) + amnc * ((cosangle * (-sin_theta(itheta)) + d_cosangle_dtheta * cos_theta(itheta)) *   cos_phi(iphi)  + (d_cosangle_dphi * (-sin_theta(itheta)) + d2_cosangle_dtheta_dphi * cos_theta(itheta)) * sin_phi(iphi))
             d2Zdthetadphi(itheta,iphi) = d2Zdthetadphi(itheta,iphi) + amnc * (d_cosangle_dphi *   cos_theta(itheta)  + d2_cosangle_dtheta_dphi * sin_theta(itheta))
             
             d2Xdphi2(itheta,iphi) = d2Xdphi2(itheta,iphi) + amnc * (cosangle * (-cos_phi(iphi)) + 2*d_cosangle_dphi * (-sin_phi(iphi)) + d2_cosangle_dphi2 * cos_phi(iphi)) * cos_theta(itheta)
             d2Ydphi2(itheta,iphi) = d2Ydphi2(itheta,iphi) + amnc * (cosangle * (-sin_phi(iphi)) + 2*d_cosangle_dphi *   cos_phi(iphi)  + d2_cosangle_dphi2 * sin_phi(iphi)) * cos_theta(itheta)
             d2Zdphi2(itheta,iphi) = d2Zdphi2(itheta,iphi) + amnc * d2_cosangle_dphi2 * sin_theta(itheta)
            
          end do
       end do
    end do

    print *,"End of big imn loop"
    print *,"R:",R
    print *,"Z:",Z

    if (.not. computing_Jacobian) then
       !X_base = X
       !Y_base = Y
       Z_base = Z
       R_base = R
       dXdtheta_base = dXdtheta
       dYdtheta_base = dYdtheta
       dZdtheta_base = dZdtheta
       dXdphi_base = dXdphi
       dYdphi_base = dYdphi
       dZdphi_base = dZdphi
       d2Xdtheta2_base = d2Xdtheta2
       d2Ydtheta2_base = d2Ydtheta2
       d2Zdtheta2_base = d2Zdtheta2
       d2Xdphi2_base = d2Xdphi2
       d2Ydphi2_base = d2Ydphi2
       d2Zdphi2_base = d2Zdphi2
       d2Xdthetadphi_base = d2Xdthetadphi
       d2Ydthetadphi_base = d2Ydthetadphi
       d2Zdthetadphi_base = d2Zdthetadphi
    end if

    do iphi = 1, N_phi
       phi_copied = phi(iphi)
       call qfm_surfaces_compute_B(N_theta, R(:,iphi), phi_copied, Z(:,iphi), B_R(:,iphi), B_phi(:,iphi), B_Z(:,iphi))
       B_X(:,iphi) = B_R(:,iphi) * cos_phi(iphi) - B_phi(:,iphi) * sin_phi(iphi)
       B_Y(:,iphi) = B_R(:,iphi) * sin_phi(iphi) + B_phi(:,iphi) * cos_phi(iphi)
    end do

    print *,"B_R:",B_R
    print *,"B_phi:",B_phi
    print *,"B_Z:",B_Z
        
    NX = normal_vector_sign * (dYdtheta * dZdphi - dYdphi * dZdtheta)
    NY = normal_vector_sign * (dZdtheta * dXdphi - dZdphi * dXdtheta)
    NZ = normal_vector_sign * (dXdtheta * dYdphi - dXdphi * dYdtheta)
    do iphi = 1, N_phi
       NR(:,iphi) = NX(:,iphi) * cos_phi(iphi) + NY(:,iphi) * sin_phi(iphi)
    end do
        
    norm_normal = sqrt(NX*NX + NY*NY + NZ*NZ)
    normal_X = NX / norm_normal
    normal_Y = NY / norm_normal
    normal_Z = NZ / norm_normal
        
    Bnormal = B_X * normal_X + B_Y * normal_Y + B_Z * normal_Z
        
    fundamental_form_E = dXdtheta * dXdtheta + dYdtheta * dYdtheta + dZdtheta * dZdtheta
    fundamental_form_F = dXdtheta * dXdphi  + dYdtheta * dYdphi  + dZdtheta * dZdphi
    fundamental_form_G = dXdphi  * dXdphi  + dYdphi  * dYdphi  + dZdphi  * dZdphi
        
    fundamental_form_L = d2Xdtheta2     * normal_X + d2Ydtheta2     * normal_Y + d2Zdtheta2     * normal_Z
    fundamental_form_M = d2Xdthetadphi * normal_X + d2Ydthetadphi * normal_Y + d2Zdthetadphi * normal_Z
    fundamental_form_N = d2Xdphi2      * normal_X + d2Ydphi2      * normal_Y + d2Zdphi2      * normal_Z
        
    mean_curvature = (fundamental_form_E * fundamental_form_N + fundamental_form_G * fundamental_form_L - 2  * fundamental_form_F * fundamental_form_M) &
         / (2*(fundamental_form_E * fundamental_form_G - fundamental_form_F * fundamental_form_F))
        
    ! For the next few lines, see note 20180712-01.
    B_dot_e_theta = B_X * dXdtheta + B_Y * dYdtheta + B_Z * dZdtheta
    B_dot_e_phi  = B_X * dXdphi  + B_Y * dYdphi  + B_Z * dZdphi
        
    N_dot_B_cross_e_phi  = normal_vector_sign * (B_dot_e_theta * fundamental_form_G - B_dot_e_phi  * fundamental_form_F)
    N_dot_e_theta_cross_B = normal_vector_sign * (B_dot_e_phi  * fundamental_form_E - B_dot_e_theta * fundamental_form_F)

    ! Use BLAS here for speed?
    d_Bnormal_d_theta = matmul(ddtheta,Bnormal)
    d_Bnormal_d_phi  = transpose(matmul(ddphi,transpose(Bnormal)))
        
    shape_gradient = mean_curvature * Bnormal * Bnormal &
         + normal_vector_sign * (N_dot_B_cross_e_phi * d_Bnormal_d_theta + N_dot_e_theta_cross_B * d_Bnormal_d_phi) / (norm_normal * norm_normal)

    residual = 0
        
    integrand = (shape_gradient + state_vector(vector_size)) * norm_normal
    do imn = 1,mnmax
       ! Need to add explicit theta + phi loops here
       do iphi = 1, N_phi
          do itheta = 1, N_theta
             !!sinangle = sin(m*theta-n*phi) = sin(m*theta) * cos(n*phi) - cos(m*theta) * sin(n*phi)
             !sinangle = sin_m_theta(itheta,m+1) * cos_n_phi(iphi,n+1) - cos_m_theta(itheta,m+1) * sin_n_phi(iphi,n+1)
             !!cosangle = cos(m*theta-n*phi) = cos(m*theta) * cos(n*phi) + sin(m*theta) * sin(n*phi)
             cosangle = cos_m_theta(itheta,xm(imn)) * cos_n_phi(iphi,xn(imn)) + sin_m_theta(itheta,xm(imn)) * sin_n_phi(iphi,xn(imn))
             !angle = xm(imn) * theta2D - xn(imn) * phi2D;
             !residual(imn) = dtheta * dphi * nfp * sum(sum(integrand .* cos(angle)));
             residual(imn) = residual(imn) + integrand(itheta,iphi) * cosangle
          end do
       end do
    end do
    residual = residual * dtheta * dphi * nfp
        
    ! Compute plasma volume using \int (1/2) R^2 dZ dphi.
    residual(vector_size) = 0.5d+0 * dtheta * dphi * nfp * sum(R * R * dZdtheta) - volume_target

  end subroutine qfm_surfaces_residual

end subroutine qfm_surfaces_single_solve
