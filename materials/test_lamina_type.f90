program test_lamina_type
! Purpose:
! to perform unit testing on lamina_type_module
!    
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    03/04/15  B. Y. Chen            Original code
!
! Pseudo code:
! declare local variables
! 
!
!
!
!
!
use parameter_module,   only : 
use lamina_type_module, only : lamina_modulus, lamina_strength, &
& lamina_matrixToughness, lamina_fibreToughness, lamina_type,   &
& empty, update, extract, ddsdde

implicit none

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
! maxdm

type(lamina_type)             :: this
type(lamina_modulus)          :: modulus
type(lamina_strength)         :: strength
type(lamina_matrixToughness)  :: matrixToughness
type(lamina_fibreToughness)   :: fibreToughness
type(lamina_sdv)              :: sdv
real(DP)                      :: dee(:,:)
real(DP)                      :: stress(:), strain(:)
real(DP)                      :: clength
integer                       :: istat
character(MSGLENGTH)          :: emsg
real(DP)                      :: maxdm

! initialize local variables
! all derived types have been initialized in definition
dee     = ZERO
stress  = ZERO
strain  = ZERO
clength = ZERO
istat   = STAT_SUCCESS
emsg    = ''
maxdm   = ZERO

! define values of input modulus
modulus%E1   = 
modulus%E2   =
modulus%G12  =
modulus%G23  =
modulus%nu12 =
modulus%nu23 =

! define values of input strengths
strength%Xt  =
strength%Xc  =
strength%Yt  =
strength%Yc  =
strength%St  =
strength%Sl  =

! define values of input matrix toughness
matrixToughness%GIc  =
matrixToughness%GIIc =
matrixToughness%eta  =

! define values of input fibre  toughness
fibreToughness%GfcT  =
fibreToughness%GfcC  =

! define strain
strain(1) =

! define clength
clength =

! define maxdm
maxdm =


! call all public procedures, test their correctness

call empty(this)

call update (this, modulus, strength, fibreToughness,&
  & matrixToughness)
  
call extract (this, modulus, strength, matrixToughness,&
  & fibreToughness)

call ddsdde (this, dee, stress, sdv, strain, clength, &
& istat, emsg, maxdm)

! check to see if outputs are correct
! variables to check:
! dee
! stress
! sdv
! istat
! emsg

end program test_lamina_type