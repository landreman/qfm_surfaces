subroutine qfm_surfaces_volume_scan

  use qfm_surfaces_variables

  implicit none

  integer :: j_volume, j

  N_resolutions = max(max_mpol, max_ntor)
  allocate(mpols(N_resolutions))
  allocate(ntors(N_resolutions))
  mpols = [( 1 + int((max_mpol-1)*(j-1.0d+0)/(N_resolutions-1)), j=1, N_resolutions )]
  ntors = [( 1 + int((max_ntor-1)*(j-1.0d+0)/(N_resolutions-1)), j=1, N_resolutions )]

  print "(a,*(i4))"," mpols:",mpols
  print "(a,*(i4))"," ntors:",ntors

  allocate(volumes(N_volumes))
  allocate(quadratic_flux(N_volumes))

  volumes = [( max_volume * sqrt((j_volume * 1.0d+0)/N_volumes), j_volume = 1, N_volumes )]
  print *,"Quadratic-flux-minimizing surfaces with the following volumes will be computed:"
  print "(*(es10.3))",volumes

  do j_volume = 1, N_volumes
     call qfm_surfaces_single_solve(j_volume)
  end do

end subroutine qfm_surfaces_volume_scan
