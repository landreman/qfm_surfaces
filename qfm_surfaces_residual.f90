subroutine qfm_surfaces_residual(Jacobian_column)

  use qfm_surfaces_variables
  use qfm_surfaces_solver_variables

  implicit none

  integer, intent(in) :: Jacobian_column
  integer :: mnmax_to_use
  logical :: computing_Jacobian
  integer :: itheta, iphi
  real(dp) :: this_amnc, d_cosangle_dtheta, d_cosangle_dphi, d2_cosangle_dtheta2, d2_cosangle_dphi2, d2_cosangle_dtheta_dphi
  real :: segment_start_time, segment_end_time

  call cpu_time(segment_start_time)

  !print "(a,i3)","Entering qfm_surfaces_residual. Jacobian_column=",Jacobian_column
  !print *,"state_vector:",state_vector
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
        dZdphi(itheta,:) = d_Z0_d_phi
        d2Xdphi2(itheta,:) = d2_R0_d_phi2 * cos_phi + 2 * d_R0_d_phi * (-sin_phi) + R0 * (-cos_phi)
        d2Ydphi2(itheta,:) = d2_R0_d_phi2 * sin_phi + 2 * d_R0_d_phi * cos_phi + R0 * (-sin_phi)
        d2Zdphi2(itheta,:) = d2_Z0_d_phi2
     end do
     mnmax_to_use = mnmax
  end if

  !print *,"R0:",R0
  !print *,"Z0:",Z0
  !print *,"cos_theta:",cos_theta

  call cpu_time(segment_end_time)
  init_time = init_time + segment_end_time - segment_start_time
  call cpu_time(segment_start_time)

  do imn = 1, mnmax_to_use
     if (computing_Jacobian) then
        m = xm(Jacobian_column)
        n = xn(Jacobian_column)
        this_amnc = epsilon
     else
        m = xm(imn)
        n = xn(imn)
        this_amnc = state_vector(imn)
     end if
     !print *,"imn=",imn,"amnc=",this_amnc
     do iphi = 1, N_phi
        do itheta = 1, N_theta
           !sinangle = sin(m*theta-n*phi) = sin(m*theta) * cos(n*phi) - cos(m*theta) * sin(n*phi)
           sinangle = sin_m_theta(itheta,m) * cos_n_phi(iphi,n) - cos_m_theta(itheta,m) * sin_n_phi(iphi,n)
           !cosangle = cos(m*theta-n*phi) = cos(m*theta) * cos(n*phi) + sin(m*theta) * sin(n*phi)
           cosangle = cos_m_theta(itheta,m) * cos_n_phi(iphi,n) + sin_m_theta(itheta,m) * sin_n_phi(iphi,n)

           d_cosangle_dtheta = -m*sinangle
           d_cosangle_dphi  =  n*nfp*sinangle
           d2_cosangle_dtheta2 = -m*m*cosangle
           d2_cosangle_dphi2  = -n*n*nfp*nfp*cosangle
           d2_cosangle_dtheta_dphi = m*n*nfp*cosangle

           R(itheta,iphi) = R(itheta,iphi) + this_amnc * cosangle * cos_theta(itheta)
           !X(itheta,iphi) = X(itheta,iphi) + this_amnc * cosangle * cos_theta(itheta) * cos_phi(iphi)
           !Y(itheta,iphi) = Y(itheta,iphi) + this_amnc * cosangle * cos_theta(itheta) * sin_phi(iphi)
           Z(itheta,iphi) = Z(itheta,iphi) + this_amnc * cosangle * sin_theta(itheta)

           dXdtheta(itheta,iphi) = dXdtheta(itheta,iphi) + this_amnc * (cosangle * (-sin_theta(itheta)) + d_cosangle_dtheta * cos_theta(itheta)) * cos_phi(iphi)
           dYdtheta(itheta,iphi) = dYdtheta(itheta,iphi) + this_amnc * (cosangle * (-sin_theta(itheta)) + d_cosangle_dtheta * cos_theta(itheta)) * sin_phi(iphi)
           dZdtheta(itheta,iphi) = dZdtheta(itheta,iphi) + this_amnc * (cosangle *   cos_theta(itheta)  + d_cosangle_dtheta * sin_theta(itheta))

           dXdphi(itheta,iphi) = dXdphi(itheta,iphi) + this_amnc * (cosangle * (-sin_phi(iphi)) + d_cosangle_dphi * cos_phi(iphi)) * cos_theta(itheta)
           dYdphi(itheta,iphi) = dYdphi(itheta,iphi) + this_amnc * (cosangle *   cos_phi(iphi)  + d_cosangle_dphi * sin_phi(iphi)) * cos_theta(itheta)
           dZdphi(itheta,iphi) = dZdphi(itheta,iphi) + this_amnc * d_cosangle_dphi * sin_theta(itheta)

           d2Xdtheta2(itheta,iphi) = d2Xdtheta2(itheta,iphi) + this_amnc * (cosangle * (-cos_theta(itheta)) + 2*d_cosangle_dtheta * (-sin_theta(itheta)) + d2_cosangle_dtheta2 * cos_theta(itheta)) * cos_phi(iphi)
           d2Ydtheta2(itheta,iphi) = d2Ydtheta2(itheta,iphi) + this_amnc * (cosangle * (-cos_theta(itheta)) + 2*d_cosangle_dtheta * (-sin_theta(itheta)) + d2_cosangle_dtheta2 * cos_theta(itheta)) * sin_phi(iphi)
           d2Zdtheta2(itheta,iphi) = d2Zdtheta2(itheta,iphi) + this_amnc * (cosangle * (-sin_theta(itheta)) + 2*d_cosangle_dtheta *   cos_theta(itheta)  + d2_cosangle_dtheta2 * sin_theta(itheta))

           d2Xdthetadphi(itheta,iphi) = d2Xdthetadphi(itheta,iphi) + this_amnc * ((cosangle * (-sin_theta(itheta)) + d_cosangle_dtheta * cos_theta(itheta)) * (-sin_phi(iphi)) + (d_cosangle_dphi * (-sin_theta(itheta)) + d2_cosangle_dtheta_dphi * cos_theta(itheta)) * cos_phi(iphi))
           d2Ydthetadphi(itheta,iphi) = d2Ydthetadphi(itheta,iphi) + this_amnc * ((cosangle * (-sin_theta(itheta)) + d_cosangle_dtheta * cos_theta(itheta)) *   cos_phi(iphi)  + (d_cosangle_dphi * (-sin_theta(itheta)) + d2_cosangle_dtheta_dphi * cos_theta(itheta)) * sin_phi(iphi))
           d2Zdthetadphi(itheta,iphi) = d2Zdthetadphi(itheta,iphi) + this_amnc * (d_cosangle_dphi *   cos_theta(itheta)  + d2_cosangle_dtheta_dphi * sin_theta(itheta))

           d2Xdphi2(itheta,iphi) = d2Xdphi2(itheta,iphi) + this_amnc * (cosangle * (-cos_phi(iphi)) + 2*d_cosangle_dphi * (-sin_phi(iphi)) + d2_cosangle_dphi2 * cos_phi(iphi)) * cos_theta(itheta)
           d2Ydphi2(itheta,iphi) = d2Ydphi2(itheta,iphi) + this_amnc * (cosangle * (-sin_phi(iphi)) + 2*d_cosangle_dphi *   cos_phi(iphi)  + d2_cosangle_dphi2 * sin_phi(iphi)) * cos_theta(itheta)
           d2Zdphi2(itheta,iphi) = d2Zdphi2(itheta,iphi) + this_amnc * d2_cosangle_dphi2 * sin_theta(itheta)

        end do
     end do
  end do

  call cpu_time(segment_end_time)
  compute_derivatives_time = compute_derivatives_time + segment_end_time - segment_start_time
  call cpu_time(segment_start_time)

  !print *,"End of big imn loop"
  !print *,"R:",R
  !print *,"Z:",Z

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

  !print *,"B_R:",B_R
  !print *,"B_phi:",B_phi
  !print *,"B_Z:",B_Z

  NX = normal_vector_sign * (dYdtheta * dZdphi - dYdphi * dZdtheta)
  NY = normal_vector_sign * (dZdtheta * dXdphi - dZdphi * dXdtheta)
  NZ = normal_vector_sign * (dXdtheta * dYdphi - dXdphi * dYdtheta)
  do iphi = 1, N_phi
     NR(:,iphi) = NX(:,iphi) * cos_phi(iphi) + NY(:,iphi) * sin_phi(iphi)
  end do
  !print *,"NX:",NX
  !print *,"NY:",NY
  !print *,"NZ:",NZ
  !print *,"dXdtheta:",dXdtheta
  !print *,"dZdtheta:",dZdtheta
  !print *,"dXdphi:",dXdphi
  !print *,"dZdphi:",dZdphi

  norm_normal = sqrt(NX*NX + NY*NY + NZ*NZ)
  normal_X = NX / norm_normal
  normal_Y = NY / norm_normal
  normal_Z = NZ / norm_normal
  !print *,"norm_normal:",norm_normal
  !print *,"normal_X:",normal_X

  Bnormal = B_X * normal_X + B_Y * normal_Y + B_Z * normal_Z
  !print *,"Bnormal:",Bnormal

  fundamental_form_E = dXdtheta * dXdtheta + dYdtheta * dYdtheta + dZdtheta * dZdtheta
  fundamental_form_F = dXdtheta * dXdphi  + dYdtheta * dYdphi  + dZdtheta * dZdphi
  fundamental_form_G = dXdphi  * dXdphi  + dYdphi  * dYdphi  + dZdphi  * dZdphi

  fundamental_form_L = d2Xdtheta2     * normal_X + d2Ydtheta2     * normal_Y + d2Zdtheta2     * normal_Z
  fundamental_form_M = d2Xdthetadphi * normal_X + d2Ydthetadphi * normal_Y + d2Zdthetadphi * normal_Z
  fundamental_form_N = d2Xdphi2      * normal_X + d2Ydphi2      * normal_Y + d2Zdphi2      * normal_Z

  mean_curvature = (fundamental_form_E * fundamental_form_N + fundamental_form_G * fundamental_form_L - 2  * fundamental_form_F * fundamental_form_M) &
       / (2*(fundamental_form_E * fundamental_form_G - fundamental_form_F * fundamental_form_F))

  !print *,"mean_curvature:",mean_curvature

  ! For the next few lines, see note 20180712-01.
  B_dot_e_theta = B_X * dXdtheta + B_Y * dYdtheta + B_Z * dZdtheta
  B_dot_e_phi  = B_X * dXdphi  + B_Y * dYdphi  + B_Z * dZdphi

  N_dot_B_cross_e_phi  = normal_vector_sign * (B_dot_e_theta * fundamental_form_G - B_dot_e_phi  * fundamental_form_F)
  N_dot_e_theta_cross_B = normal_vector_sign * (B_dot_e_phi  * fundamental_form_E - B_dot_e_theta * fundamental_form_F)

  call cpu_time(segment_end_time)
  multiply_time = multiply_time + segment_end_time - segment_start_time
  call cpu_time(segment_start_time)

  ! Use BLAS here for speed?
  d_Bnormal_d_theta = matmul(ddtheta,Bnormal)
  d_Bnormal_d_phi  = transpose(matmul(ddphi,transpose(Bnormal)))

  call cpu_time(segment_end_time)
  matmul_time = matmul_time + segment_end_time - segment_start_time
  call cpu_time(segment_start_time)

  shape_gradient = mean_curvature * Bnormal * Bnormal &
       + normal_vector_sign * (N_dot_B_cross_e_phi * d_Bnormal_d_theta + N_dot_e_theta_cross_B * d_Bnormal_d_phi) / (norm_normal * norm_normal)

  residual = 0
  integrand = (shape_gradient + state_vector(vector_size)) * norm_normal

  call cpu_time(segment_end_time)
  multiply_time = multiply_time + segment_end_time - segment_start_time
  call cpu_time(segment_start_time)

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

  call cpu_time(segment_end_time)
  transform_time = transform_time + segment_end_time - segment_start_time
  call cpu_time(segment_start_time)

  ! Compute plasma volume using \int (1/2) R^2 dZ dphi.
  residual(vector_size) = 0.5d+0 * dtheta * dphi * nfp * sum(R * R * dZdtheta) - volume_target

  call cpu_time(segment_end_time)
  multiply_time = multiply_time + segment_end_time - segment_start_time

end subroutine qfm_surfaces_residual
