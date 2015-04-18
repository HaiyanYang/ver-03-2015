include 'globals/parameter_module.f90'
include 'materials/cohesive_material_module.f90'

program test_cohesive_material
! Purpose:
! to perform unit testing on cohesive_material_module
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    06/04/15  B. Y. Chen            Original code
!    11/04/15  B. Y. Chen            added testing of cohesive_ig_point
!
!
use parameter_module, NST => NST_COHESIVE
use cohesive_material_module

implicit none

! declare variables:
! this
! modulus
! strength
! toughness
! sdv
! dee
! traction, separation
! istat
! emsg
! d_max

type(cohesive_material)       :: this
type(cohesive_sdv)            :: sdv
real(DP)                      :: dee(NST,NST)
real(DP)                      :: traction(NST), separation(NST)
integer                       :: istat
character(MSGLENGTH)          :: emsg
real(DP)                      :: d_max

type(cohesive_ig_point) :: ig_point
real(DP)           :: igx(NDIM), igu(NDIM), igtract(NST), igsepar(NST)
type(cohesive_sdv) :: igsdv1, igsdv2

character(len=20) :: display_fmt
character(len=10) :: cndim, cnst

logical :: nofailure


! initialize local variables
! all derived types have been initialized in definition
dee         = ZERO
traction    = ZERO
separation  = ZERO
istat       = STAT_SUCCESS
emsg        = ''
d_max       = ZERO

display_fmt = ''
cndim = ''
cnst  = ''

nofailure   = .false.

igx     = ZERO
igu     = ZERO
igtract = ZERO
igsepar = ZERO

! define separation
separation(1) = -0.006_DP
separation(2) = 0.02_DP
separation(3) = 0.02_DP

! define d_max
d_max = ONE

! store ndim and nst as strings in cndim and cnst
write(cndim,'(i5)') NDIM
write(cnst,'(i5)') NST

! call all public procedures, test their correctness

call empty (this)

call display (this)

call set (this, modulus=cohesive_modulus(Dnn=1000000._DP, Dtt=1000000._DP, Dll=1000000._DP), &
        & strength=cohesive_strength(tau_nc=60._DP, tau_tc=90._DP, tau_lc=90._DP), &
        & toughness=cohesive_toughness(Gnc=0.2_DP, Gtc=1._DP, Glc=1._DP, alpha=1._DP), &
        & istat=istat, emsg=emsg)

if(istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if

call display (this)

!nofailure = .true.
nofailure = .false.

if (nofailure) then

  call ddsdde (this, dee=dee, traction=traction, separation=separation)

else

  call ddsdde (this, dee=dee, traction=traction, sdv=sdv, separation=separation, &
  & istat=istat, emsg=emsg, d_max=d_max)

end if

if(istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if

! check to see if outputs are correct
! variables to check:
! traction
! sdv
! istat
! emsg

! set display format, note that for scientific real, ESw.d, w>=d+7
display_fmt = '(1X, A, ES10.3)'
write(*,'(A)') ''
write(*,'(A)') 'display the traction components:'
write(*,display_fmt) 'tau_n: ', traction(1)
write(*,display_fmt) 'tau_t: ', traction(2)
write(*,display_fmt) 'tau_l: ', traction(3)
write(*,'(A)') ''

call display (sdv)

write(*,'(A)') ''
write(*,'(A)') 'display the status and message:'
write(*,*) istat
write(*,*) emsg
write(*,'(A)') ''




write(*,'(A)') ''
write(*,'(A)') 'test cohesive_ig_point:'
write(*,'(A)') ''
write(*,'(A)') 'display the cohesive_ig_point before any update:'
call display(ig_point)

igx = ONE
igu = ONE

call update(ig_point, x=igx, u=igu)
write(*,'(A)') 'display the cohesive_ig_point after updates on x and u:'
call display(ig_point)

igtract = ONE
igsepar = ONE
!igsdv1=cohesive_sdv(dm=0.5_DP, u0=0.5_DP, uf=1.5_DP, fstat=COH_MAT_ONSET)
!igsdv2=cohesive_sdv(dm=0.9_DP, u0=0.5_DP, uf=1.5_DP, fstat=COH_MAT_FAILED)

call update(ig_point, traction=igtract, separation=igsepar, &
& converged_sdv=igsdv1, iterating_sdv=igsdv2)
write(*,'(A)') 'display the cohesive_ig_point after all updates:'
call display(ig_point)

write(*,'(A)') ''
write(*,'(A)') 'check extracted values from ig_point'

igx      = ZERO
igu      = ZERO
igtract = ZERO
igsepar = ZERO
!igsdv1 = cohesive_sdv(ZERO, ZERO, ZERO, INTACT)
!igsdv2 = cohesive_sdv(ZERO, ZERO, ZERO, INTACT)

call extract(ig_point, x=igx, u=igu, traction=igtract, separation=igsepar, &
& converged_sdv=igsdv1, iterating_sdv=igsdv2)

write(*,'(A)') ''
write(*,'(A)') 'display extracted values from ig_point'
display_fmt = '(1X, A,'//trim(adjustl(cndim))//'ES10.3)'
write(*,display_fmt) '- x          :', igx
write(*,display_fmt) '- u          :', igu
display_fmt = '(1X, A,'//trim(adjustl(cnst))//'ES10.3)'
write(*,display_fmt) '- traction   :', igtract
write(*,display_fmt) '- separation :', igsepar
call display(igsdv1)
call display(igsdv2)

call empty(ig_point)
write(*,'(A)') 'display the cohesive_ig_point after being emptied:'
call display(ig_point)



end program test_cohesive_material
