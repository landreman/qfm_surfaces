subroutine qfm_surfaces_compute_B(N, R, phi, Z, B_R, B_phi, B_Z)

  use qfm_surfaces_variables

  implicit none

  integer, intent(in) :: N
  real(dp), dimension(N), intent(in) :: R, phi, Z
  real(dp), dimension(N), intent(out) :: B_R, B_phi, B_Z
  real(dp), dimension(N) :: sinnphi, cosnphi

!!$  print *,"size(phi):",size(phi)
!!$
!!$  N = size(phi)
!!$  if (N==0) N=1
!!$  print *,"N:",N
!!$  allocate(sinnphi(N))
!!$  allocate(cosnphi(N))

  sinnphi = sin(nfp * phi)
  cosnphi = cos(nfp * phi)

  B_R = amplitude_m1 * ((3*(-1/R**4 + R*R) * Z / 2) * cosnphi + (((1/R**4) + R*R)/2) * sinnphi) &
       + amplitude_m2 * (((5 - 2 * R*R - 5 * R**8 - 24 * Z*Z + R**6 * (2 + 24 * Z*Z))/(32 * R**4)) * sinnphi &
       + (((1 + R**6) * Z)/(2 * R**4)) * cosnphi)

  B_phi = (1 + amplitude_m1 * nfp * (((1/R**3 + R**3) * Z / 2) * (-sinnphi) + ((-(1/R**3) + R**3)/6) * cosnphi) &
       + amplitude_m2 * nfp * (((-5 + 6 *R*R + 2 *R**6 - 3 *R**8)/(96 *R**3) + (1/R**3 + R**3)/4 * Z*Z) * cosnphi &
       + ((-(1/R**3) + R**3) * Z)/6 * (-sinnphi))) / R

  B_Z = amplitude_m1 * (((1/R**3 + R**3) / 2) * cosnphi + 0 * sinnphi) &
       + amplitude_m2 * (((1/R**3 + R**3) * Z)/2 * sinnphi &
       + ((-(1/R**3) + R**3)/6) * cosnphi)

!  deallocate(sinnphi,cosnphi)

end subroutine qfm_surfaces_compute_B
