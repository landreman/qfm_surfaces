module qfm_surfaces_solver_variables

  use stel_kinds

  implicit none

  integer :: j_resolution, m, n, j, k, imn
  integer :: N_theta, N_phi, mpol, ntor
  real(dp), dimension(:), allocatable :: theta, phi, temp_1D
  real(dp), dimension(:,:), allocatable :: ddtheta, ddphi, temp_2D, Jacobian
  real(dp), dimension(:), allocatable :: R0, d_R0_d_phi, d2_R0_d_phi2
  real(dp), dimension(:), allocatable :: Z0, d_Z0_d_phi, d2_Z0_d_phi2
  real(dp), dimension(:), allocatable :: state_vector, state_vector0, last_state_vector, residual, residual0, step_direction
  real(dp) :: angle, sinangle, cosangle, volume_target, dtheta, dphi
  integer :: vector_size, last_vector_size=0, index, jm, jn
  integer, dimension(:), allocatable :: last_xm, last_xn
  logical :: found_match
  real(dp), dimension(:,:), allocatable :: sin_m_theta, cos_m_theta, sin_n_phi, cos_n_phi
  real(dp), dimension(:), allocatable :: sin_phi, cos_phi, sin_theta, cos_theta
  real(dp), dimension(:,:), allocatable :: Z, R, dXdtheta, dYdtheta, dZdtheta, dXdphi, dYdphi, dZdphi
  real(dp), dimension(:,:), allocatable :: d2Xdtheta2, d2Ydtheta2, d2Zdtheta2, d2Xdphi2, d2Ydphi2, d2Zdphi2
  real(dp), dimension(:,:), allocatable :: d2Xdthetadphi, d2Ydthetadphi, d2Zdthetadphi
  real(dp), dimension(:,:), allocatable :: Z_base, R_base, dXdtheta_base, dYdtheta_base, dZdtheta_base, dXdphi_base, dYdphi_base, dZdphi_base
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
  real :: start_time, end_time
  real :: init_time, compute_derivatives_time, multiply_time, matmul_time, transform_time
  real(dp) :: this_area, this_quadratic_flux

  ! Variables needed by LAPACK:                                                                                            
  integer :: INFO
  integer, dimension(:), allocatable :: IPIV

end module qfm_surfaces_solver_variables
