subroutine qfm_surfaces_init_solve()

  use qfm_surfaces_variables
  use qfm_surfaces_solver_variables

  implicit none

  mpol = mpols(j_resolution)
  ntor = ntors(j_resolution)
  N_theta = max(min_N_theta,mpol*4)
  N_phi = max(min_N_phi,ntor*4)
  if (verbose) print "(a,i5,a,i5,a,i6,a,i6,a)","  Using resolution mpol=",mpol,", ntor=",ntor,", N_theta=",N_theta,", N_phi=",N_phi,"."
  max_mpol_used = max(max_mpol_used,mpol)
  max_ntor_used = max(max_ntor_used,ntor)

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
  call qfm_surfaces_differentiation_matrix(N_phi, 0.0d+0, 2*pi/nfp, 20, 1, phi, temp_1D, ddphi, temp_2D)
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
  allocate(IPIV(vector_size))
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
!!$        print *,"last_xm:",last_xm
!!$        print *,"     xm:",xm
!!$        print *,"last_xn:",last_xn
!!$        print *,"     xn:",xn
!!$        print *,"last_state_vector:",last_state_vector
!!$        print *,"     state_vector:",state_vector

     deallocate(last_xm, last_xn, last_state_vector)
  end if


end subroutine qfm_surfaces_init_solve
