subroutine qfm_surfaces_volume_scan

  use qfm_surfaces_variables

  implicit none

  integer :: j_volume, j, index, jm, jn, mpol, ntor

  allocate(amnc_big(0:max_mpol, -max_ntor:max_ntor, N_volumes))
  amnc_big = 0

  N_resolutions = max(max_mpol, max_ntor)
  allocate(mpols(N_resolutions))
  allocate(ntors(N_resolutions))
  mpols = [( 1 + int((max_mpol-1)*(j-1.0d+0)/(N_resolutions-1)), j=1, N_resolutions )]
  ntors = [( 1 + int((max_ntor-1)*(j-1.0d+0)/(N_resolutions-1)), j=1, N_resolutions )]
  max_mpol_used = 0
  max_ntor_used = 0

  print "(a,*(i4))"," mpols:",mpols
  print "(a,*(i4))"," ntors:",ntors

  allocate(volumes(N_volumes))
  allocate(areas(N_volumes))
  allocate(quadratic_flux(N_volumes))
  allocate(quadratic_flux_convergence(N_resolutions,N_volumes))
  allocate(shape_convergence(N_resolutions,N_shape_convergence_locations,N_volumes))
  allocate(lambda(N_volumes))

  volumes = [( max_volume * ((j_volume * 1.0d+0)/N_volumes)**2, j_volume = 1, N_volumes )]
  print *,"Quadratic-flux-minimizing surfaces with the following volumes will be computed:"
  print "(*(es10.3))",volumes

  do j_volume = 1, N_volumes
     call qfm_surfaces_single_volume(j_volume)
  end do

  ! Trim down amnc_big to the arrays saved to output
  mpol = max_mpol_used
  ntor = max_ntor_used
  mnmax = (ntor*2 + 1) * mpol + ntor + 1
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
  allocate(amnc(mnmax,N_volumes))
  do index = 1, mnmax
     amnc(index,:) = amnc_big(xm(index),xn(index),:)
  end do

end subroutine qfm_surfaces_volume_scan
