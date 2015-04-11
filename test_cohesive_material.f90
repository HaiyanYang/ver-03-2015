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
use parameter_module
use cohesive_material_module

implicit none

integer, parameter :: NDIM=3, NST = 3

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
type(cohesive_modulus)        :: modulus
type(cohesive_strength)       :: strength
type(cohesive_toughness)      :: toughness
type(cohesive_sdv)            :: sdv
real(DP)                      :: dee(NST,NST)
real(DP)                      :: traction(NST), separation(NST)
integer                       :: istat
character(MSGLENGTH)          :: emsg
real(DP)                      :: d_max

type(cohesive_ig_point)         :: ig_point
real(DP), allocatable           :: igx(:), igu(:), igstress(:), igstrain(:)
type(cohesive_sdv), allocatable :: igsdv(:)

character(len=20) :: display_fmt

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
nofailure   = .false.

allocate(igx(NDIM),igu(NDIM),igstress(NST),igstrain(NST),igsdv(2))
igx      = ZERO
igu      = ZERO
igstress = ZERO
igstrain = ZERO

! define values of input modulus
modulus%Dnn  = 1000000._DP
modulus%Dtt  = 1000000._DP
modulus%Dll  = 1000000._DP

! define values of input strengths
strength%tau_nc = 60._DP
strength%tau_tc = 90._DP
strength%tau_lc = 90._DP

! define values of input toughness
toughness%Gnc   = 0.2_DP
toughness%Gtc   = 1._DP
toughness%Glc   = 1._DP
toughness%alpha = 1._DP

! define separation
separation(1) = -0.006_DP
separation(2) = 0.02_DP
separation(3) = 0.02_DP

! define d_max
d_max = ONE


! call all public procedures, test their correctness

call empty (this)

call display (this)

call set (this, modulus=modulus, strength=strength, toughness=toughness, &
& istat=istat, emsg=emsg)

if(istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if

call display (this)

!nofailure = .true.
nofailure = .false.

if (nofailure) then

  call ddsdde (this, dee=dee, traction=traction, separation=separation, &
  & istat=istat, emsg=emsg)

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

call update(ig_point, istat=istat, emsg=emsg, x=igx, u=igu)
if(istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if
write(*,'(A)') 'display the cohesive_ig_point after updates on x and u:'
call display(ig_point)

call update(ig_point, istat=istat, emsg=emsg, stress=igstress, strain=igstrain, sdv=igsdv)
if(istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if
write(*,'(A)') 'display the cohesive_ig_point after all updates:'
call display(ig_point)

call empty(ig_point)
write(*,'(A)') 'display the cohesive_ig_point after being emptied:'
call display(ig_point)



end program test_cohesive_material
