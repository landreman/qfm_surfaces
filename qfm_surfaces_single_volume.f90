subroutine qfm_surfaces_single_volume(j_volume)

  use qfm_surfaces_variables
  use qfm_surfaces_solver_variables

  implicit none

  integer, intent(in) :: j_volume

  volume_target = volumes(j_volume)
  print "(a,i5,a,i5,a)"," Solving for volume",j_volume," of",N_volumes,"."
  if (verbose) print *,"volume_target=",volume_target

  do j_resolution = 1, N_resolutions
     call qfm_surfaces_init_solve()

     call qfm_surfaces_Newton()

     allocate(last_state_vector(vector_size))
     allocate(last_xm(mnmax))
     allocate(last_xn(mnmax))
     last_vector_size = vector_size
     last_state_vector = state_vector
     last_xm = xm
     last_xn = xn

     do imn = 1, mnmax
        amnc_big(xm(imn), xn(imn), j_volume) = state_vector(imn)
     end do
     lambda(j_volume) = state_vector(vector_size)

     call qfm_surfaces_flux()
     areas(j_volume) = this_area
     quadratic_flux(j_volume) = this_quadratic_flux
     quadratic_flux_convergence(j_resolution,j_volume) = this_quadratic_flux

     call qfm_surfaces_deallocate()

  end do ! Loop over resolution
  deallocate(last_xm, last_xn, last_state_vector)


end subroutine qfm_surfaces_single_volume
