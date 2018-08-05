subroutine qfm_surfaces_volume_scan

  use qfm_surfaces_variables

  implicit none

  include 'mpif.h'

  integer :: j_volume, j, index, jm, jn, mpol, ntor
  integer :: ierr, tag, scan_index_min, scan_index_max, N_volumes_local
  integer :: mpi_status(MPI_STATUS_SIZE), send_data(1), recv_data(1), dummy(1)
  integer, parameter :: buffer_length = 100
  character(len=buffer_length) :: proc_assignments_string
  real :: start_time, end_time


  if (N_volumes >= N_procs) then
     scan_index_min = 1 + (mpi_rank * N_volumes) / N_procs
     scan_index_max = ((mpi_rank+1) * N_volumes) / N_procs
  else
     ! There are more procs than solves to do
     scan_index_min = min(mpi_rank+1, N_volumes)
     if (mpi_rank < N_volumes) then
        scan_index_max = scan_index_min
     else
        scan_index_max = scan_index_min - 1
     end if
  end if
  N_volumes_local = scan_index_max - scan_index_min + 1

  write (proc_assignments_string,fmt="(a,i5,a,i5,a,i9,a,i9)") "Proc ",mpi_rank," of",N_procs," will handle radii",scan_index_min," to",scan_index_max

  ! Print the processor/radius assignments in a coordinated manner.
  if (proc0) then
     print *,trim(proc_assignments_string)
     do j = 1,N_procs - 1
        ! Use mpi rank of the sender as the tag
        tag = j
        call MPI_RECV(proc_assignments_string,buffer_length,MPI_CHAR,j,tag,MPI_COMM_WORLD,mpi_status,ierr)
        print *,trim(proc_assignments_string)
     end do
  else
     tag = mpi_rank
     call MPI_SEND(proc_assignments_string,buffer_length,MPI_CHAR,0,tag,MPI_COMM_WORLD,ierr)
  end if

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  allocate(amnc_big(0:max_mpol, -max_ntor:max_ntor, N_volumes))
  amnc_big = 0

  N_resolutions = max(max_mpol, max_ntor)
  allocate(mpols(N_resolutions))
  allocate(ntors(N_resolutions))
  mpols = [( 1 + int((max_mpol-1)*(j-1.0d+0)/(N_resolutions-1)), j=1, N_resolutions )]
  ntors = [( 1 + int((max_ntor-1)*(j-1.0d+0)/(N_resolutions-1)), j=1, N_resolutions )]
  max_mpol_used = 0
  max_ntor_used = 0

  if (proc0) then
     print "(a,*(i4))"," mpols:",mpols
     print "(a,*(i4))"," ntors:",ntors
  end if

  allocate(volumes(N_volumes))
  allocate(areas(N_volumes))
  allocate(quadratic_flux(N_volumes))
  allocate(quadratic_flux_convergence(N_resolutions,N_volumes))
  allocate(shape_convergence(N_resolutions,N_shape_convergence_locations,N_volumes))
  allocate(lambda(N_volumes))

  if (min_volume <= 0) then
     volumes = [( max_volume * ((j_volume * 1.0d+0)/N_volumes)**2, j_volume = 1, N_volumes )]
  else
     volumes = [( (sqrt(min_volume) + (sqrt(max_volume) - sqrt(min_volume)) * (j_volume - 1.0d+0) / (N_volumes - 1.0d+0) ) ** 2, j_volume = 1, N_volumes )]
  end if
  if (proc0) print *,"Quadratic-flux-minimizing surfaces with the following volumes will be computed:"
  if (proc0) print "(*(es10.3))",volumes

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call cpu_time(start_time)

  !do j_volume = 1, N_volumes
  do j_volume = scan_index_min, scan_index_max
     call qfm_surfaces_single_volume(j_volume)
  end do

  call cpu_time(end_time)
  print "(a,i4,a,es10.3,a)"," proc",mpi_rank," finished, took",end_time-start_time," seconds."

  ! Send results to proc 0
  if (proc0) then
     do j = 1, N_procs-1
        ! Use mpi_rank as the tag
        call mpi_recv(dummy,1,MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        scan_index_min = dummy(1)
        call mpi_recv(dummy,1,MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        scan_index_max = dummy(1)
        N_volumes_local = scan_index_max - scan_index_min + 1

        ! Receive 1D arrays
        call mpi_recv(areas(scan_index_min:scan_index_max),N_volumes_local,MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(lambda(scan_index_min:scan_index_max),N_volumes_local,MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(quadratic_flux(scan_index_min:scan_index_max),N_volumes_local,MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)

        ! Receive multi-D arrays
        call mpi_recv(quadratic_flux_convergence(:,scan_index_min:scan_index_max),N_resolutions*N_volumes_local,MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(amnc_big(:,:,scan_index_min:scan_index_max),(max_mpol+1)*(max_ntor*2+1)*N_volumes_local,MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
     end do
  else
     send_data = scan_index_min
     ! Use mpi_rank as the tag
     call mpi_send(send_data,1,MPI_INT,0,mpi_rank,MPI_COMM_WORLD,ierr)

     send_data = scan_index_max
     call mpi_send(send_data,1,MPI_INT,0,mpi_rank,MPI_COMM_WORLD,ierr)

     ! Send 1D arrays
     call mpi_send(areas(scan_index_min:scan_index_max),N_volumes_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(lambda(scan_index_min:scan_index_max),N_volumes_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(quadratic_flux(scan_index_min:scan_index_max),N_volumes_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)

     ! Send multi-D arrays
     call mpi_send(quadratic_flux_convergence(:,scan_index_min:scan_index_max),N_resolutions*N_volumes_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(amnc_big(:,:,scan_index_min:scan_index_max),(max_mpol+1)*(max_ntor*2+1)*N_volumes_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)

  end if

  send_data = max_mpol_used
  call mpi_reduce(send_data,recv_data,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  mpol = recv_data(1)

  send_data = max_ntor_used
  call mpi_reduce(send_data,recv_data,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  ntor = recv_data(1)

  ! Trim down amnc_big to the arrays saved to output
  !mpol = max_mpol_used
  !ntor = max_ntor_used
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
