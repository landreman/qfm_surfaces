subroutine qfm_surfaces_flux

  use qfm_surfaces_variables
  use qfm_surfaces_solver_variables

  implicit none

  integer :: mnmax_to_use
  integer :: itheta, iphi
  real(dp) :: this_amnc, d_cosangle_dtheta, d_cosangle_dphi

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
  ! Add contribution from the axis
  do itheta = 1, N_theta
     !X(itheta,:) = R0 * cos_phi
     !Y(itheta,:) = R0 * sin_phi
     Z(itheta,:) = Z0
     R(itheta,:) = R0
     dXdphi(itheta,:) = d_R0_d_phi * cos_phi + R0 * (-sin_phi)
     dYdphi(itheta,:) = d_R0_d_phi * sin_phi + R0 * cos_phi
     dZdphi(itheta,:) = d_Z0_d_phi
  end do
  mnmax_to_use = mnmax

  !print *,"R0:",R0
  !print *,"Z0:",Z0
  !print *,"cos_theta:",cos_theta

  do imn = 1, mnmax_to_use
     m = xm(imn)
     n = xn(imn)
     this_amnc = state_vector(imn)
     !print *,"imn=",imn,"amnc=",this_amnc
     do iphi = 1, N_phi
        do itheta = 1, N_theta
           !sinangle = sin(m*theta-n*phi) = sin(m*theta) * cos(n*phi) - cos(m*theta) * sin(n*phi)
           sinangle = sin_m_theta(itheta,m) * cos_n_phi(iphi,n) - cos_m_theta(itheta,m) * sin_n_phi(iphi,n)
           !cosangle = cos(m*theta-n*phi) = cos(m*theta) * cos(n*phi) + sin(m*theta) * sin(n*phi)
           cosangle = cos_m_theta(itheta,m) * cos_n_phi(iphi,n) + sin_m_theta(itheta,m) * sin_n_phi(iphi,n)

           d_cosangle_dtheta = -m*sinangle
           d_cosangle_dphi  =  n*nfp*sinangle

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

        end do
     end do
  end do

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

  norm_normal = sqrt(NX*NX + NY*NY + NZ*NZ)
  normal_X = NX / norm_normal
  normal_Y = NY / norm_normal
  normal_Z = NZ / norm_normal
  !print *,"norm_normal:",norm_normal
  !print *,"normal_X:",normal_X

  Bnormal = B_X * normal_X + B_Y * normal_Y + B_Z * normal_Z
  !print *,"Bnormal:",Bnormal

  this_quadratic_flux = dtheta * dphi * nfp * sum(Bnormal * Bnormal * norm_normal)
  this_area = dtheta * dphi * nfp * sum(norm_normal)

end subroutine qfm_surfaces_flux
