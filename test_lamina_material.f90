include 'globals/parameter_module.f90'
include 'materials/lamina_material_module.f90'

program test_lamina_material
! Purpose:
! to perform unit testing on lamina_material_module
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    03/04/15  B. Y. Chen            Original code
!    11/04/15  B. Y. Chen            added testing of lamina_ig_point
!
!
use parameter_module, NST => NST_STANDARD
use lamina_material_module

implicit none

! declare variables:
! this
! modulus
! strength
! fibreToughness
! sdv
! dee
! stress, strain
! clength
! istat
! emsg
! d_max

type(lamina_material)         :: this
type(lamina_modulus)          :: modulus
type(lamina_strength)         :: strength
type(lamina_fibreToughness)   :: fibreToughness
type(lamina_sdv)              :: sdv
real(DP)                      :: dee(NST,NST)
real(DP)                      :: stress(NST), strain(NST)
real(DP)                      :: clength
integer                       :: istat
character(MSGLENGTH)          :: emsg
real(DP)                      :: d_max

type(lamina_ig_point) :: ig_point
real(DP)              :: igx(NDIM), igu(NDIM), igstress(NST), igstrain(NST)
type(lamina_sdv)      :: igsdv1, igsdv2

character(len=20) :: display_fmt
character(len=10) :: cndim, cnst

logical :: nofailure


! initialize local variables
! all derived types have been initialized in definition
dee     = ZERO
stress  = ZERO
strain  = ZERO
clength = ZERO
istat   = STAT_SUCCESS
emsg    = ''
d_max   = ZERO

display_fmt = ''
cndim = ''
cnst  = ''

nofailure   = .false.

igx      = ZERO
igu      = ZERO
igstress = ZERO
igstrain = ZERO

! define values of input modulus
modulus%E1   = 161000._DP
modulus%E2   = 11400._DP
modulus%G12  = 5170._DP
modulus%G23  = 4980._DP
modulus%nu12 = ZERO
modulus%nu23 = ZERO

! define values of input strengths
strength%Xt  = 2820._DP
strength%Xc  = 1140._DP
strength%Yt  = 60._DP
strength%Yc  = 110._DP
strength%St  = 90._DP
strength%Sl  = 90._DP

! define values of input fibre  toughness
fibreToughness%GfcT  = 112.7_DP
fibreToughness%GfcC  = 112.7_DP

! define strain
strain(1) = 0.6_DP

! define clength
clength = 0.5_DP

! define d_max
d_max = ONE

! store ndim and nst as strings in cndim and cnst
write(cndim,'(i5)') NDIM
write(cnst,'(i5)') NST


! call all public procedures, test their correctness

call empty (this)

call display (this)

call set (this, modulus, strength, fibreToughness, &
& istat=istat, emsg=emsg)

if(istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if

call display (this)

!nofailure = .true.
nofailure = .false.

if (nofailure) then
  ! ddsdde_intact
  call ddsdde (this, dee=dee, stress=stress, strain=strain)

else

  call ddsdde (this, dee=dee, stress=stress, sdv=sdv, strain=strain, &
  & clength=clength, istat=istat, emsg=emsg, d_max=d_max)

end if


if(istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if

! check to see if outputs are correct
! variables to check:
! stress
! sdv
! istat
! emsg

! set display format, note that for scientific real, ESw.d, w>=d+7
display_fmt = '(1X, A, ES10.3)'
write(*,'(A)') ''
write(*,'(A)') 'display the stress components:'
write(*,display_fmt) 'sigma_1: ', stress(1)
write(*,display_fmt) 'sigma_2: ', stress(2)
write(*,display_fmt) 'sigma_3: ', stress(3)
write(*,display_fmt) 'tau_12 : ', stress(4)
write(*,display_fmt) 'tau_13 : ', stress(5)
write(*,display_fmt) 'tau_23 : ', stress(6)
write(*,'(A)') ''

call display (sdv)

write(*,'(A)') ''
write(*,'(A)') 'display the status and message:'
write(*,*) istat
write(*,*) emsg
write(*,'(A)') ''

write(*,'(A)') ''
write(*,'(A)') 'test lamina_ig_point:'
write(*,'(A)') ''
write(*,'(A)') 'display the lamina_ig_point before any update:'
call display(ig_point)

igx      = ONE
igu      = ONE

call update(ig_point, x=igx, u=igu)
write(*,'(A)') 'display the lamina_ig_point after updates on x and u:'
call display(ig_point)

igstress = ONE
igstrain = ONE
igsdv1=lamina_sdv(df=0.5_DP, u0=0.5_DP, uf=1.5_DP, &
     & fstat=FIBRE_ONSET, ffstat=FIBRE_ONSET, mfstat=MATRIX_ONSET)
igsdv2=lamina_sdv(df=0.9_DP, u0=0.5_DP, uf=1.5_DP, &
     & fstat=FIBRE_FAILED, ffstat=FIBRE_FAILED, mfstat=MATRIX_ONSET)

call update(ig_point, stress=igstress, strain=igstrain, &
& converged_sdv=igsdv1, iterating_sdv=igsdv2)
write(*,'(A)') 'display the lamina_ig_point after all updates:'
call display(ig_point)

write(*,'(A)') ''
write(*,'(A)') 'check extracted values from ig_point'

igx      = ZERO
igu      = ZERO
igstress = ZERO
igstrain = ZERO
igsdv1 = lamina_sdv(ZERO, ZERO, ZERO, INTACT, INTACT, INTACT)
igsdv2 = lamina_sdv(ZERO, ZERO, ZERO, INTACT, INTACT, INTACT)

call extract(ig_point, x=igx, u=igu, stress=igstress, strain=igstrain, &
& converged_sdv=igsdv1, iterating_sdv=igsdv2)

write(*,'(A)') ''
write(*,'(A)') 'display extracted values from ig_point'
display_fmt = '(1X, A,'//trim(adjustl(cndim))//'ES10.3)'
write(*,display_fmt) '- x      :', igx
write(*,display_fmt) '- u      :', igu
display_fmt = '(1X, A,'//trim(adjustl(cnst))//'ES10.3)'
write(*,display_fmt) '- stress :', igstress
write(*,display_fmt) '- strain :', igstrain
call display(igsdv1)
call display(igsdv2)

call empty(ig_point)
write(*,'(A)') 'display the lamina_ig_point after being emptied:'
call display(ig_point)

end program test_lamina_material
