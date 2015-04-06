include 'globals/parameter_module.f90'
include 'materials/global_matkey_module.f90'

program test_global_matkey
! Purpose:
! to perform unit testing on global_matkey_module
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    06/04/15  B. Y. Chen            Original code
!
!
use parameter_module
use global_matkey_module, only : global_matkey, empty, set, extract, display

implicit none

integer, parameter :: NST = 3

type(global_matkey)          :: this
character(len=MATNAMELENGTH) :: matname
character(len=MATTYPELENGTH) :: mattype
integer                      :: type_array_index
integer                      :: istat
character(len=MSGLENGTH)     :: emsg

! initialize local variables
! all derived types have been initialized in definition
matname          = ''
mattype          = ''
type_array_index = 0
istat            = STAT_SUCCESS
emsg             = ''

! define inputs
matname = 'user material 1'
mattype = 'cohesive'
type_array_index = 10

! call all public procedures, test their correctness

call empty (this)

call display (this)

call set (this, matname=matname, mattype=mattype, &
& type_array_index=type_array_index, istat=istat, emsg=emsg)

if(istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if

call display (this)

matname          = ''
mattype          = ''
type_array_index = 0

call extract (this, matname=matname, mattype=mattype, &
& type_array_index=type_array_index)

write(*,'(1X,a)')   ''
write(*,'(1X,a)')   trim(matname)
write(*,'(1X,a)')   trim(mattype)
write(*,'(1X,I10)') type_array_index

end program test_global_matkey
