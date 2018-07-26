subroutine qfm_surfaces_Jacobian()

  use qfm_surfaces_variables
  use qfm_surfaces_solver_variables

  implicit none

  call cpu_time(start_time)
  residual0 = residual
  state_vector0 = state_vector

  do j = 1, vector_size
     if (j==vector_size) state_vector(vector_size) = state_vector(vector_size) + epsilon
     ! Aside from the last column, the perturbations to the state vector are handled inside qfm_surfaces_residual for the sake of speed.
     call qfm_surfaces_residual(j)
     Jacobian(:,j) = (residual - residual0) / epsilon
  end do

  residual = residual0
  state_vector = state_vector0
  call cpu_time(end_time)
  print *,"Time to compute finite-difference Jacobian:",end_time-start_time

end subroutine qfm_surfaces_Jacobian

