subroutine qfm_surfaces_deallocate

  use qfm_surfaces_variables
  use qfm_surfaces_solver_variables

  implicit none


  ! Done with this resolution. Deallocate everything.
  deallocate(phi, ddphi, theta, ddtheta)
  deallocate(R0, d_R0_d_phi, d2_R0_d_phi2, Z0, d_Z0_d_phi, d2_Z0_d_phi2)
  deallocate(xm, xn, state_vector, state_vector0, residual, residual0, Jacobian, step_direction, IPIV)
  deallocate(sin_m_theta, cos_m_theta, sin_n_phi, cos_n_phi, sin_phi, cos_phi, sin_theta, cos_theta)
  deallocate(Z,R)
  deallocate(dXdtheta,dYdtheta,dZdtheta,dXdphi,dYdphi,dZdphi)
  deallocate(d2Xdtheta2,d2Ydtheta2,d2Zdtheta2,d2Xdphi2,d2Ydphi2,d2Zdphi2)
  deallocate(d2Xdthetadphi,d2Ydthetadphi,d2Zdthetadphi)
  deallocate(Z_base,R_base)
  deallocate(dXdtheta_base,dYdtheta_base,dZdtheta_base,dXdphi_base,dYdphi_base,dZdphi_base)
  deallocate(d2Xdtheta2_base,d2Ydtheta2_base,d2Zdtheta2_base,d2Xdphi2_base,d2Ydphi2_base,d2Zdphi2_base)
  deallocate(d2Xdthetadphi_base,d2Ydthetadphi_base,d2Zdthetadphi_base)
  deallocate(B_R, B_phi, B_Z, B_X, B_Y, NX, NY, NZ, NR, normal_X, normal_Y, normal_Z, Bnormal, norm_normal)
  deallocate(fundamental_form_E, fundamental_form_F, fundamental_form_G, fundamental_form_L, fundamental_form_M, fundamental_form_N)
  deallocate(mean_curvature, shape_gradient, B_dot_e_theta, B_dot_e_phi, N_dot_B_cross_e_phi, N_dot_e_theta_cross_B)
  deallocate(d_Bnormal_d_theta, d_Bnormal_d_phi, integrand, phi_copied)


end subroutine qfm_surfaces_deallocate
