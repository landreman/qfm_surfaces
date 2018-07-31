subroutine qfm_surfaces_read_input

  use qfm_surfaces_variables

  implicit none

  integer :: numargs
  character(len=200) :: input_filename
  integer :: fileUnit, didFileAccessWork, i
  integer, parameter :: uninitialized = -9999

  ! getcarg is in LIBSTELL
  call getcarg(1, input_filename, numargs)

  if (numargs<1) then
     stop "One argument is required: the input namelist file, which must be named qfm_surfaces_in.XXXXX"
  end if
  if (numargs>1) then
     print *,"WARNING: Arguments after the first will be ignored."
  end if
  if (input_filename(1:16) .ne. "qfm_surfaces_in.") then
     stop "Input file must be named qfm_surfaces_in.XXX for some extension XXX"
  end if

  output_filename = "qfm_surfaces_out" // trim(input_filename(16:)) // ".nc"

  fileUnit=11
  open(unit=fileUnit, file=input_filename, action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error opening input file ", trim(input_filename)
     stop
  else
     read(fileUnit, nml=qfm_surfaces, iostat=didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Error!  I was able to open the file ", trim(input_filename), &
               " but not read data from the qfm_surfaces namelist in it."
        if (didFileAccessWork==-1) then
           print *,"Make sure there is a carriage return after the / at the end of the namelist!"
        end if
        stop
     end if
     if (proc0) print *,"Successfully read parameters from qfm_surfaces namelist in ", trim(input_filename), "."
  end if
  close(unit = fileUnit)


end subroutine qfm_surfaces_read_input
