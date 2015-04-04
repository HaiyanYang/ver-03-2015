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
!
!
use parameter_module
use lamina_material_module, only : lamina_modulus, lamina_strength, &
& lamina_matrixToughness, lamina_fibreToughness, lamina_material,   &
& lamina_sdv, empty, update, extract, display, ddsdde

implicit none

integer, parameter :: NST = 6

! declare variables:
! this
! modulus
! strength
! fibreToughness
! matrixToughness
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
type(lamina_matrixToughness)  :: matrixToughness
type(lamina_fibreToughness)   :: fibreToughness
type(lamina_sdv)              :: sdv
real(DP)                      :: dee(NST,NST)
real(DP)                      :: stress(NST), strain(NST)
real(DP)                      :: clength
integer                       :: istat
character(MSGLENGTH)          :: emsg
real(DP)                      :: d_max

character(len=20) :: display_fmt


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

! define values of input modulus
modulus%E1   = 161000._DP
modulus%E2   = 11400._DP
modulus%G12  = 5170._DP
modulus%G23  = 4980._DP
modulus%nu12 = 0.32_DP
modulus%nu23 = 0.43_DP

! define values of input strengths
strength%Xt  = 2820._DP
strength%Xc  = 1140._DP
strength%Yt  = 60._DP
strength%Yc  = 110._DP
strength%St  = 90._DP
strength%Sl  = 90._DP

! define values of input matrix toughness
matrixToughness%GIc  = 0.293_DP
matrixToughness%GIIc = 0.631_DP
matrixToughness%eta  = 1._DP

! define values of input fibre  toughness
fibreToughness%GfcT  = 112.7_DP
fibreToughness%GfcC  = 112.7_DP

! define strain
strain(1) = 0.02_DP

! define clength
clength = 0.5_DP

! define d_max
d_max = ONE


! call all public procedures, test their correctness

call empty(this)

call display (this)

call update (this, modulus, strength, fibreToughness,&
  & matrixToughness)

call display (this)
!~
!~call extract (this, modulus, strength, matrixToughness,&
!~  & fibreToughness)
!~
call ddsdde (this, dee, stress, sdv, strain, clength, &
& istat, emsg, d_max)

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

end program test_lamina_material
