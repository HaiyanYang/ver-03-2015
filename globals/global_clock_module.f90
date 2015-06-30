module global_clock_module
!
!  Purpose:
!   the global progress of analysis; the global clock object is reset  
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
  private
  integer  :: step_number      = 0
  integer  :: increment_number = 0
end type program_clock

! define the global variable GLOBAL_CLOCK, it is saved
type(program_clock), public, save :: GLOBAL_CLOCK


! related procedures
interface empty
  module procedure empty_program_clock
end interface empty

interface set
  module procedure set_program_clock
end interface set

public :: empty, set, clock_in_sync




contains




pure subroutine empty_program_clock (GLOBAL_CLOCK)

  type(program_clock), intent(inout) :: GLOBAL_CLOCK

  GLOBAL_CLOCK%step_number      = 0
  GLOBAL_CLOCK%increment_number = 0

end subroutine empty_program_clock



pure subroutine set_program_clock (GLOBAL_CLOCK, curr_step, curr_inc, &
& istat, emsg)
! Purpose:
! this module sets the current global clock values
use parameter_module, only : MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

  type(program_clock), intent(inout) :: GLOBAL_CLOCK
  integer,             intent(in)    :: curr_step, curr_inc
  integer,                  optional, intent(out) :: istat
  character(len=MSGLENGTH), optional, intent(out) :: emsg
  
  if (present(istat) .and. present(emsg)) then
    istat = STAT_SUCCESS
    emsg  = ''
    ! check validity of inputs
    if (curr_step < 1 .or. curr_inc < 1) then
      istat = STAT_FAILURE
      emsg  = 'current step or increment no. is < 1, set, global_clock_module'
      return
    end if
  end if
  
  ! update if inputs are valid
  GLOBAL_CLOCK%step_number      = curr_step
  GLOBAL_CLOCK%increment_number = curr_inc

end subroutine set_program_clock



pure logical function clock_in_sync (GLOBAL_CLOCK, local_clock)
! Purpose:
! to see if the local clock is in sync with the global clock
  type(program_clock), intent(in) :: GLOBAL_CLOCK
  type(program_clock), intent(in) :: local_clock
  
  if (local_clock%step_number == GLOBAL_CLOCK%step_number .and. &
  &   local_clock%increment_number == GLOBAL_CLOCK%increment_number) then
    clock_in_sync = .true.
  else
    clock_in_sync = .false.
  end if
  
end function clock_in_sync




end module global_clock_module