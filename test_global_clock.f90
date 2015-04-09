include 'globals/global_clock_module.f90'

program test_global_clock
! Purpose:
! to perform unit testing on global_clock_module
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    09/04/15  B. Y. Chen            Original code
!
!
use global_clock_module ! use everything

implicit none

type(program_clock) :: local_clock

! initialize local variables
! all derived types have been initialized in definition

! define inputs
local_clock%step_number = 2
local_clock%increment_number = 10

! call all public procedures, test their correctness
call empty (global_clock)

call set (global_clock, curr_step=1, curr_inc=10)

if (clock_in_sync(global_clock, local_clock)) then
  write(*,'(1X,a)')   ''
  write(*,'(1X,a)')   'local clock is in sync with global clock'
  write(*,'(1X,a)')   ''
else
  write(*,'(1X,a)')   ''
  write(*,'(1X,a)')   'local clock is out of sync with global clock'
  write(*,'(1X,a)')   ''
  local_clock = global_clock
  if (clock_in_sync(global_clock, local_clock)) &
  & write(*,'(1X,a)')   'local clock is in sync with global clock'
end if

end program test_global_clock
