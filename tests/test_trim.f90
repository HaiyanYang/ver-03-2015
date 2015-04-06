PROGRAM test_trim
  CHARACTER(len=10) :: s = "GFORTRAN  "
  WRITE(*,*) LEN(s), LEN(TRIM(s)), LEN(s)
END PROGRAM