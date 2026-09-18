Program timmes_reference_driver
  Include 'implno.dek'
  Include 'vector_eos.dek'

  Integer, Parameter :: nstate = 5
  Integer :: i
  Real(8) :: abar(nstate), den(nstate), temp(nstate), zbar(nstate)

  temp = (/ 3.2d7, 1.7d8, 8.0d8, 4.0d9, 9.0d9 /)
  den = (/ 2.5d2, 3.0d6, 2.0d9, 7.0d7, 1.0d10 /)
  abar = (/ 1.0d0/(0.70d0 + 0.28d0/4.0d0 + 0.02d0/12.0d0), &
    & 4.0d0, 12.0d0, 56.0d0, &
    & 1.0d0/(0.50d0/12.0d0 + 0.50d0/56.0d0) /)
  zbar = (/ abar(1)*(0.70d0 + 2.0d0*0.28d0/4.0d0 + 6.0d0*0.02d0/12.0d0), &
    & 2.0d0, 6.0d0, 26.0d0, &
    & abar(5)*(6.0d0*0.50d0/12.0d0 + 26.0d0*0.50d0/56.0d0) /)

  jlo_eos = 1
  jhi_eos = nstate
  Do i = 1, nstate
    temp_row(i) = temp(i)
    den_row(i) = den(i)
    abar_row(i) = abar(i)
    zbar_row(i) = zbar(i)
  EndDo
  Call eosfxt

  Write(*,'(a)') 'state temp_K rho_g_cm3 abar zbar cv_erg_g_K eta detadt_per_K'
  Do i = 1, nstate
    Write(*,'(i2,7(1x,es25.16e3))') i,temp(i),den(i),abar(i),zbar(i), &
      & cv_row(i),etaele_row(i),detat_row(i)
  EndDo
End Program timmes_reference_driver
