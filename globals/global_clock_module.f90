module global_clock_module
!
!  Purpose:
!   the global progress of analysis; this module is continuously updated  
!   at the start of each increment       
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    09/04/15  B. Y. Chen            Original code
!

implicit none
private

! define a program_clock object
type, public :: program_clock
  integer  :: step_number      = 0
  integer  :: increment_number = 0
end type program_clock

! define the global program_clock, publicly available and saved
type(program_clock), public, save :: global_clock


! related procedures
interface empty
  module procedure empty_global_clock
end interface empty

interface set
  module procedure set_global_clock
end interface set

public :: empty, set, clock_in_sync




contains




pure subroutine empty_global_clock (global_clock)

  type(program_clock), intent(inout) :: global_clock

  global_clock%step_number      = 0
  global_clock%increment_number = 0

end subroutine empty_global_clock



pure subroutine set_global_clock (global_clock, curr_step, curr_inc)

  type(program_clock), intent(inout) :: global_clock
  integer,  optional,  intent(in)    :: curr_step, curr_inc
  
  if (present(curr_step)) global_clock%step_number      = curr_step
  if (present(curr_inc))  global_clock%increment_number = curr_inc

end subroutine set_global_clock



pure logical function clock_in_sync (global_clock, local_clock)
! Purpose:
! to see if the local clock is in sync with the global clock
  type(program_clock), intent(in) :: global_clock
  type(program_clock), intent(in) :: local_clock
  
  if (local_clock%step_number == global_clock%step_number .and. &
  &   local_clock%increment_number == global_clock%increment_number) then
    clock_in_sync = .true.
  else
    clock_in_sync = .false.
  end if
  
end function clock_in_sync




end module global_clock_module