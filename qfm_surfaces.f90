! Main program

program qfm_surfaces

  use stel_kinds
  use qfm_surfaces_variables, only: proc0, mpi_rank, N_procs, total_time

  implicit none

  include 'mpif.h'

  integer :: ierr
!integer :: tic, toc, countrate, ierr
  real :: start_time, end_time
  real(dp) :: B_R, B_phi, B_Z

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
  proc0 = (mpi_rank==0)

  if (proc0) then
     print "(a)"," -------------------------------------------------------------"
     print *,"Computing quadratic-flux-minimizing surfaces."
  end if
  !call system_clock(tic,countrate)
  call cpu_time(start_time)

  call qfm_surfaces_read_input()
  call qfm_surfaces_validate_input()

!!$  call qfm_surfaces_compute_B(1.1d+0, 0.2d+0, -0.05d+0, B_R, B_phi, B_Z)
!!$  print *,"B_R:",B_R
!!$  print *,"B_phi:",B_phi
!!$  print *,"B_Z:",B_Z

  call qfm_surfaces_compute_axis()

  call qfm_surfaces_volume_scan()

  !call system_clock(toc)
  !total_time = real(toc-tic)/countrate
  call cpu_time(end_time)
  total_time = end_time - start_time

  call qfm_surfaces_write_output()

  if (proc0) then
     print "(a)"," -------------------------------------------------------------"
     print "(a,es10.3,a)"," QFM-surfaces solver is complete. Total time=",total_time," sec."
  end if

  call mpi_finalize(ierr)

end program qfm_surfaces
