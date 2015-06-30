module global_matkey_module
!
!  Purpose:
!    define an object which functions as a key to access a particular
!    material definition in the global material library; it comes
!    with the associated procedures to empty, set, and extract
!    its components
!    
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    06/04/15  B. Y. Chen            Original code
!

! list out ALL the parameter used, no matter how long:
! MATNAMELENGTH : character length of a material name
! MATTYPELENGTH : character length of a material type name
! MSGLENGTH     : character length of error message
! STAT_SUCCESS  : integer value of a successful status
! STAT_FAILURE  : integer value of a failure    status
use parameter_module, only : MATNAMELENGTH, MATTYPELENGTH, &
                           & MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

implicit none
private

! define a key to access a material definition in the global material library
type, public :: global_matkey
  private
  ! type components:
  ! matname           : user-given name from input
  ! mattype           : defined material types: lamina or cohesive
  ! type_list_index   : index in the array of the material type
  character(len=MATNAMELENGTH) :: matname          = ''  
  character(len=MATTYPELENGTH) :: mattype          = ''  
  integer                      :: type_list_index  = 0
end type

interface empty
  module procedure empty_global_matkey
end interface

interface set
  module procedure set_global_matkey
end interface

interface extract
  module procedure extract_global_matkey
end interface

interface display
  module procedure display_global_matkey
end interface

public :: empty, set, extract, display



contains



  pure subroutine empty_global_matkey (this)
  ! Purpose:
  ! this subroutine is used in the preprocessing to format the 
  ! global material library lib_mat
    type(global_matkey), intent(inout):: this
    
    this%matname = ''
    this%mattype = ''
    this%type_list_index = 0

  end subroutine empty_global_matkey



  pure subroutine set_global_matkey (this, matname, mattype, type_list_index, &
  & istat, emsg)
  ! Purpose:
  ! this subroutine is used in the preprocessing to fill in the
  ! material information in the global material library lib_mat
  ! property check subroutine, error status and message are needed to 
  ! flag an error when the input properties are unacceptable
  
    ! - istat       : status variable of this procedure   to output
    ! - emsg        : error message                       to output
    type(global_matkey),      intent(inout) :: this
    character(len=*),         intent(in)    :: matname 
    character(len=*),         intent(in)    :: mattype 
    integer,                  intent(in)    :: type_list_index
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    
    ! initialize intent(out) & local variables
    istat = STAT_SUCCESS  ! default
    emsg  = ''
    
    this%matname = matname
    this%mattype = mattype
    this%type_list_index = type_list_index
    
    ! check this_mat properties
    call check_prop (this, istat, emsg)
    if (istat == STAT_FAILURE) return

  end subroutine set_global_matkey
 
 
  
  pure subroutine check_prop (this, istat, emsg)
  ! Purpose:
  ! to check the validity of the properties of the pass arg. this
    
    ! - istat       : status variable of this procedure   to output
    ! - emsg        : error message                       to output
    type(global_matkey),      intent(in)  :: this
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    
    ! initialize intent(out) & local variables
    istat = STAT_SUCCESS  ! default
    emsg  = ''
    
    ! check the material type name, see if it is defined
    select case (trim(this%mattype))
    
      case ('lamina','cohesive')
        continue
        
      case default
      ! unsupported material type name
      ! flag error status and message and exit
        istat = STAT_FAILURE
        emsg  = 'unsupported material type name, global_matkey_module; &
        & supported material types: lamina, cohesive'
        return
        
    end select
    
    
    ! check the material type_list_index, see if it is legal
    if (.not. (this%type_list_index > 0)) then
    ! type list index must be larger than 0
    ! flag error status and message and exit
      istat = STAT_FAILURE
      emsg  = 'material type list index must be > 0, global_matkey_module'
      return
    end if
  
  end subroutine check_prop



  pure subroutine extract_global_matkey (this, matname, mattype, &
  & type_list_index)
  ! Purpose:
  ! this subroutine is used anywhere to extract material information
    type(global_matkey),                    intent(in)  :: this
    character(len=MATNAMELENGTH), optional, intent(out) :: matname 
    character(len=MATTYPELENGTH), optional, intent(out) :: mattype 
    integer,                      optional, intent(out) :: type_list_index
    
    if (present(matname))          matname = this%matname
    if (present(mattype))          mattype = this%mattype
    if (present(type_list_index))  type_list_index = this%type_list_index

  end subroutine extract_global_matkey
  

  
  subroutine display_global_matkey (this)
  ! Purpose:
  ! to display this object's components on cmd window
  ! this is useful for debugging
    type(global_matkey), intent(in):: this
    
    ! local variable to set the output format
    character(len=20) :: display_fmt
    
    ! initialize local variable
    display_fmt = ''
    
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the components of the inquired global matkey:'
    write(*,'(1X, A)') 'the material name is: '//trim(this%matname)
    write(*,'(1X, A)') 'the material type is: '//trim(this%mattype)
    
    ! set display format for string and integer
    ! A for string, I for integer, 2 for width of the number
    display_fmt = '(1X, A, I10)' 
    
    write(*,display_fmt) "the material's index in its type list is: ",&
    &this%type_list_index
    write(*,'(1X, A)') ''

  end subroutine display_global_matkey


end module global_matkey_module