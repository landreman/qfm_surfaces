subroutine qfm_surfaces_write_output

  use qfm_surfaces_variables
  use ezcdf

  implicit none

  integer :: ierr, ncid

  ! Same convention as in VMEC:
  ! Prefix vn_ indicates the variable name used in the .nc file.

  ! Scalars:
  character(len=*), parameter :: &
       vn_nfp = "nfp", &
       vn_N_volumes = "N_volumes", &
       vn_Newton_tolerance = "Newton_tolerance", &
       vn_N_iterations = "N_iterations", &
       vn_N_line_search = "N_line_search", &
       vn_Newton_tolerance_axis = "Newton_tolerance_axis", &
       vn_N_iterations_axis = "N_iterations_axis", &
       vn_N_line_search_axis = "N_line_search_axis", &
       vn_max_mpol = "max_mpol", &
       vn_max_ntor = "max_ntor"

  ! Arrays with dimension 1
  character(len=*), parameter :: &
       vn_volumes = "volumes", &
       vn_areas = "areas", &
       vn_quadratic_flux = "quadratic_flux", &
       vn_lambda = "lambda", &
       vn_R0c = "R0c", &
       vn_R0s = "R0s", &
       vn_Z0c = "Z0c", &
       vn_Z0s = "Z0s", &
       vn_xm = "xm", &
       vn_xn = "xn"
       

  ! Arrays with dimension 2
  character(len=*), parameter :: &
       vn_N_scan_array  = "N_scan_array", &
       vn_scan_R0c  = "scan_R0c", &
       vn_scan_R0s  = "scan_R0s", &
       vn_scan_Z0c  = "scan_Z0c", &
       vn_scan_Z0s  = "scan_Z0s", &
       vn_amnc = "amnc", &
       vn_quadratic_flux_convergence = "quadratic_flux_convergence"

!!$  ! Arrays with dimension 3
!!$  character(len=*), parameter :: &
!!$       vn_r_plasma  = "r_plasma"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now create variables that name the dimensions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Arrays with dimension 1:
  character(len=*), parameter, dimension(1) :: &
       N_volumes_dim = (/'N_volumes'/), &
       nmax_axis_plus_1_dim = (/'nmax_axis_plus_1'/), &
       mnmax_dim = (/'mnmax'/)

  ! Arrays with dimension 2:
  ! The form of the array declarations here is inspired by
  ! http://stackoverflow.com/questions/21552430/gfortran-does-not-allow-character-arrays-with-varying-component-lengths
  character(len=*), parameter, dimension(2) :: &
       mnmax_N_volumes_dim = (/ character(len=50) :: 'mnmax','N_volumes'/), &
       N_scan_max_axis_nmax_plus_1_dim = (/ character(len=50) :: 'N_scan','max_axis_nmax_plus_1'/), &
       N_resolutions_N_volumes_dim = (/ character(len=50) :: "N_resolutions","N_volumes" /)

!!$  ! Arrays with dimension 3:
!!$  character(len=*), parameter, dimension(3) :: &
!!$       xyz_ntheta_nzetal_plasma_dim = (/ character(len=50) :: 'xyz','ntheta_plasma','nzetal_plasma'/), &

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Only proc 0 writes.
  if (.not. proc0) return

  call cdf_open(ncid,output_filename,'w',ierr)
  IF (ierr .ne. 0) then
     print *,"Error opening output file ",output_filename
     stop
  end IF

  ! Scalars

  call cdf_define(ncid, vn_nfp, nfp)
  call cdf_setatt(ncid, vn_nfp, 'Number of field periods, i.e. the number of identical toroidal segments, 5 for W7-X, 4 for HSX, etc. ' // &
       'Equivalent to the VMEC variable of the same name.')

  call cdf_define(ncid, vn_N_volumes, N_volumes)
  !call cdf_setatt(ncid, vn_resolution_option, 'Method used to define the geometry of the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_Newton_tolerance, Newton_tolerance)
  call cdf_define(ncid, vn_N_iterations, N_iterations)
  call cdf_define(ncid, vn_N_line_search, N_line_search)
  call cdf_define(ncid, vn_Newton_tolerance_axis, Newton_tolerance_axis)
  call cdf_define(ncid, vn_N_iterations_axis, N_iterations_axis)
  call cdf_define(ncid, vn_N_line_search_axis, N_line_search_axis)
  call cdf_define(ncid, vn_max_mpol, max_mpol)
  call cdf_define(ncid, vn_max_ntor, max_ntor)
  

  ! Arrays with dimension 1

  call cdf_define(ncid, vn_volumes, volumes, dimname=N_volumes_dim)
  call cdf_define(ncid, vn_areas, areas, dimname=N_volumes_dim)
  call cdf_define(ncid, vn_quadratic_flux, quadratic_flux, dimname=N_volumes_dim)
  call cdf_define(ncid, vn_lambda, lambda, dimname=N_volumes_dim)
  call cdf_define(ncid, vn_R0c, R0c, dimname=nmax_axis_plus_1_dim)
  call cdf_define(ncid, vn_R0s, R0s, dimname=nmax_axis_plus_1_dim)
  call cdf_define(ncid, vn_Z0c, Z0c, dimname=nmax_axis_plus_1_dim)
  call cdf_define(ncid, vn_Z0s, Z0s, dimname=nmax_axis_plus_1_dim)
  call cdf_define(ncid, vn_xm, xm, dimname=mnmax_dim)
  call cdf_define(ncid, vn_xn, xn, dimname=mnmax_dim)


  ! Arrays with dimension 2

  call cdf_define(ncid, vn_amnc,  amnc, dimname=mnmax_N_volumes_dim)
  call cdf_define(ncid, vn_quadratic_flux_convergence,  quadratic_flux_convergence, dimname=N_resolutions_N_volumes_dim)
!!$  call cdf_define(ncid, vn_scan_R0c,  scan_R0c, dimname=N_scan_max_axis_nmax_plus_1_dim)
!!$  call cdf_define(ncid, vn_scan_R0s,  scan_R0s, dimname=N_scan_max_axis_nmax_plus_1_dim)
!!$  call cdf_define(ncid, vn_scan_Z0c,  scan_Z0c, dimname=N_scan_max_axis_nmax_plus_1_dim)
!!$  call cdf_define(ncid, vn_scan_Z0s,  scan_Z0s, dimname=N_scan_max_axis_nmax_plus_1_dim)

  ! Arrays with dimension 3

  !call cdf_define(ncid, vn_r_plasma,  r_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Done with cdf_define calls. Now write the data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  ! Scalars

  call cdf_write(ncid, vn_nfp, nfp)
  call cdf_write(ncid, vn_N_volumes, N_volumes)
  call cdf_write(ncid, vn_Newton_tolerance, Newton_tolerance)
  call cdf_write(ncid, vn_N_iterations, N_iterations)
  call cdf_write(ncid, vn_N_line_search, N_line_search)
  call cdf_write(ncid, vn_Newton_tolerance_axis, Newton_tolerance_axis)
  call cdf_write(ncid, vn_N_iterations_axis, N_iterations_axis)
  call cdf_write(ncid, vn_N_line_search_axis, N_line_search_axis)
  call cdf_write(ncid, vn_max_mpol, max_mpol)
  call cdf_write(ncid, vn_max_ntor, max_ntor)

  ! Arrays with dimension 1

  call cdf_write(ncid, vn_volumes, volumes)
  call cdf_write(ncid, vn_areas, areas)
  call cdf_write(ncid, vn_quadratic_flux, quadratic_flux)
  call cdf_write(ncid, vn_lambda, lambda)
  call cdf_write(ncid, vn_R0c, R0c)
  call cdf_write(ncid, vn_R0s, R0s)
  call cdf_write(ncid, vn_Z0c, Z0c)
  call cdf_write(ncid, vn_Z0s, Z0s)
  call cdf_Write(ncid, vn_xm, xm)
  call cdf_Write(ncid, vn_xn, xn)

  ! Arrays with dimension 2

  call cdf_write(ncid, vn_amnc,  amnc)
  call cdf_write(ncid, vn_quadratic_flux_convergence,  quadratic_flux_convergence)

  ! Arrays with dimension 3

  !call cdf_write(ncid, vn_r_plasma, r_plasma)

  ! Finish up:
  call cdf_close(ncid)

end subroutine qfm_surfaces_write_output
