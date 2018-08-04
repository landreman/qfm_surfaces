subroutine qfm_surfaces_Jacobian()

  use qfm_surfaces_variables
  use qfm_surfaces_solver_variables

  implicit none

  call cpu_time(start_time)
  residual0 = residual
  state_vector0 = state_vector

  init_time = 0
  compute_derivatives_time = 0
  multiply_time = 0
  matmul_time = 0
  transform_time = 0

  do j = 1, vector_size
     if (j==vector_size) state_vector(vector_size) = state_vector(vector_size) + epsilon
     ! Aside from the last column, the perturbations to the state vector are handled inside qfm_surfaces_residual for the sake of speed.
     call qfm_surfaces_residual(j)
     Jacobian(:,j) = (residual - residual0) / epsilon
  end do

  residual = residual0
  state_vector = state_vector0
  call cpu_time(end_time)
  if (trim(verbose_option)==verbose_option_detailed) then
     print *,"Time to compute finite-difference Jacobian:",end_time-start_time
     print "(5(a,es9.2))","  init:",init_time,"  compute_derivatives:",compute_derivatives_time,"  multiply:",multiply_time,"  matmul:",matmul_time,"  transform:",transform_time
  end if

end subroutine qfm_surfaces_Jacobian

