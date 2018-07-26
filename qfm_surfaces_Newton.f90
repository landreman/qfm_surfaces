subroutine qfm_surfaces_Newton()

  use qfm_surfaces_variables
  use qfm_surfaces_solver_variables

  implicit none


  ! As a prelude to the Newton iteration, get the initial residual:
  call qfm_surfaces_residual(0)
  initial_residual_norm = sqrt(sum(residual * residual))
  residual_norm = initial_residual_norm
  print "(a,es10.3)","                 Initial residual L2 norm:",residual_norm
  !print *,"residual:",residual
  !return

  ! Here is the main Newton iteration:
  Newton: do iteration = 1, N_iterations
     last_residual_norm = residual_norm
     !if (residual_norm / initial_residual_norm < Newton_tolerance) then
     if (residual_norm < Newton_tolerance) then
        exit Newton
     end if

     call qfm_surfaces_Jacobian()

     state_vector0 = state_vector
     if (verbose) print "(a,i3)","  Newton iteration ",iteration
     ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
     ! step_direction = - matrix \ residual
     step_direction = -residual ! Note that LAPACK will over-write step_direction and with the solution, and over-write Jacobian with the LU factorization.
     call cpu_time(start_time)
     call DGESV(vector_size, 1, Jacobian, vector_size, IPIV, step_direction, vector_size, INFO)
     call cpu_time(end_time)
     print *,"Time for DGESV:",end_time-start_time
     if (INFO /= 0) then
        print *, "Error in LAPACK call DGESV: info = ", INFO
        stop
     end if

     step_scale = 1
     line_search: do j_line_search = 1, N_line_search
        state_vector = state_vector0 + step_scale * step_direction

        call cpu_time(start_time)
        call qfm_surfaces_residual(0)
        call cpu_time(end_time)
        print *,"Time to compute residual:",end_time-start_time
        residual_norm = sqrt(sum(residual * residual))
        !if (verbose) print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Relative residual L2 norm:",residual_norm / initial_residual_norm,"  iota:",iota
        if (verbose) print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Residual L2 norm:",residual_norm
        if (residual_norm < last_residual_norm) exit line_search

        step_scale = step_scale / 2
     end do line_search

     if (residual_norm > last_residual_norm) then
        print *,"Line search failed to reduce residual."
        exit Newton
     end if
  end do Newton
  ! End of Newton solve.

end subroutine qfm_surfaces_Newton
