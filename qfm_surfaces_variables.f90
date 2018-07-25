module qfm_surfaces_variables

  use stel_kinds

  implicit none

  real(dp), parameter :: pi = 3.14159265358979d+0

  character(len=*), parameter :: &
       resolution_option_fixed = "fixed", &
       resolution_option_adaptive = "adaptive"
  character(len=50) :: resolution_option = resolution_option_fixed
  ! "fixed"    = Run using the specified N_phi.
  ! "adaptive" = Keep doubling N_phi (approximately, so N_phi remains odd) until iota_tolerance is achieved, or N_phi > max_N_phi.

  character(len=*), parameter :: &
       general_option_single = "single", &
       general_option_scan = "scan"
  character(len=50) :: general_option = general_option_scan

  character(len=*), parameter :: &
       constraint_option_no_Z_component = "no_Z_component", &
       constraint_option_sigma_initial = "sigma_initial"
  character(len=50) :: constraint_option = constraint_option_sigma_initial
  ! "no_Z_component" = Force the theta = 0 line to have no Z component at phi = 0.
  ! "sigma_initial"  = Force sigma = sigma_initial at phi = 0.

  character(len=200) :: output_filename

  integer :: N_procs, mpi_rank
  logical :: proc0

  integer :: nfp = 3
  integer :: N_phi_axis = 5
  real(dp) :: amplitude_m1 = 0.0
  real(dp) :: amplitude_m2 = 0.25d+0
  real(dp) :: total_time
  logical :: stellarator_symmetry = .true.
  integer :: N_iterations_axis = 20
  integer :: N_line_search_axis = 5
  real(dp) :: Newton_tolerance_axis = 1.0d-12
  integer :: N_iterations = 20
  integer :: N_line_search = 5
  real(dp) :: Newton_tolerance = 1.0d-12
  integer :: N_volumes = 1
  real(dp) :: max_volume = 1
  integer :: max_mpol = 12
  integer :: max_ntor = 9
  integer :: min_N_theta = 32
  integer :: min_N_phi = 16
  real(dp) :: min_accurate_quadratic_flux = 1.0e-9
  real(dp), dimension(:), allocatable :: volumes, quadratic_flux
  integer :: N_resolutions
  integer, dimension(:), allocatable :: mpols, ntors
  integer :: nmax_axis, max_mpol_used, max_ntor_used, mnmax
  real(dp), dimension(:), allocatable :: R0c, R0s, Z0c, Z0s, lambda
  real(dp), dimension(:,:,:), allocatable :: amnc_big
  real(dp), dimension(:,:), allocatable :: amnc
  integer, dimension(:), allocatable :: xm, xn

  namelist / qfm_surfaces / nfp, N_phi_axis, amplitude_m1, amplitude_m2, stellarator_symmetry, N_volumes, max_volume, max_mpol, max_ntor, min_accurate_quadratic_flux, &
       min_N_theta, min_N_phi, N_iterations, N_line_search, Newton_tolerance

end module qfm_surfaces_variables

