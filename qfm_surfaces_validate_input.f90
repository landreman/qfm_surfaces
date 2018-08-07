subroutine qfm_surfaces_validate_input

  use qfm_surfaces_variables

  implicit none

  if (nfp<1) stop "Error! nfp must be positive."
  !if (nfp .ne. 3) stop "nfp must be 3 for now."

  select case (trim(resolution_option))
  case (resolution_option_fixed)
  case (resolution_option_adaptive)
  case default
     print *,"Error! Invalid resolution_option:",resolution_option
  end select

  select case (trim(general_option))
  case (general_option_single)
  case (general_option_scan)
  case default
     print *,"Error! Invalid general_option:",general_option
  end select

  select case (trim(constraint_option))
  case (constraint_option_no_Z_component)
  case (constraint_option_sigma_initial)
  case default
     print *,"Error! Invalid constraint_option:",constraint_option
  end select

  select case (trim(verbose_option))
  case (verbose_option_detailed)
     verbose = .true.
  case (verbose_option_all)
     verbose = .true.
  case (verbose_option_proc0)
     verbose = proc0
  case (verbose_option_summary)
     verbose = .false.
  case default
     print *,"Error! Invalid verbose_option:",verbose_option
  end select

end subroutine qfm_surfaces_validate_input
