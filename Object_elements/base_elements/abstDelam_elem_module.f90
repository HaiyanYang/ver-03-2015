module abstDelam_elem_module
!
!  Purpose:
!    define an 'abstract element' object to interface with the
!    base coh6Delam and coh8Delam elements
!    with the associated procedures to empty, set, integrate and extract
!    its components
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    17/04/15  B. Y. Chen            Original code
!
use parameter_module, only : NDIM, NST => NST_COHESIVE, DP, ELTYPELENGTH, &
                      & MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

use coh6Delam_elem_module, only : coh6Delam_elem
use coh8Delam_elem_module, only : coh8Delam_elem

  implicit none
  private
  
  type,public :: abstDelam_elem
    ! encapsulate components of this type
    private 
    ! list of type components:
    ! eltype    : element type
    ! coh6Delam : allocated if eltype = 'coh6Delam'
    ! coh8Delam : allocated if eltype = 'coh8Delam'
    character(len=ELTYPELENGTH)       :: eltype = ''
    type(coh6Delam_elem), allocatable :: coh6Delam
    type(coh8Delam_elem), allocatable :: coh8Delam   
  end type

  interface empty
      module procedure empty_abstDelam_elem
  end interface

  interface set
      module procedure set_abstDelam_elem
  end interface
  
  interface integrate
      module procedure integrate_abstDelam_elem
  end interface
  
  interface extract
      module procedure extract_abstDelam_elem
  end interface




  public :: empty, set, integrate, extract

  
  
  
  contains
  
  
   
  
  pure subroutine empty_abstDelam_elem (elem)
  
    type(abstDelam_elem), intent(inout) :: elem
    
    elem%eltype=''

    if(allocated(elem%coh6Delam))      deallocate(elem%coh6Delam)
    if(allocated(elem%coh8Delam))      deallocate(elem%coh8Delam)
      
  end subroutine empty_abstDelam_elem
   
  
  
  pure subroutine set_abstDelam_elem (elem, eltype, connec, istat, emsg)
  ! Purpose:
  ! this is an interface used to set the contained base element according to
  ! eltype:
  ! - if eltype is coh8Delam, the input dummy args are passed to set_coh8Delam_elem
  ! to define the elem's coh8Delam component; elem's coh6Delam component is dealloc.
  ! - if eltype is coh6Delam, the input dummy args are passed to set_coh6Delam_elem
  ! to define the elem's coh6Delam component; elem's coh8Delam component is dealloc.
  ! ** note:
  ! - size of connec is checked here against eltype to ensure compatibility
  ! the values of connec and other dummy args are not checked here; they are
  ! left for the called eltype's set procedure for checking
  ! - local copies of elem's components are used for set operation;
  ! they are copied to actual elem's components only before successful return
  use coh6Delam_elem_module, only : coh6Delam_elem, set
  use coh8Delam_elem_module, only : coh8Delam_elem, set
  
      type(abstDelam_elem),     intent(inout) :: elem
      character(len=*),         intent(in)    :: eltype
      integer,                  intent(in)    :: connec(:)
      integer,                  intent(out)   :: istat
      character(len=MSGLENGTH), intent(out)   :: emsg
      
      ! local variables
      character(len=ELTYPELENGTH)       :: eltype_lcl
      type(coh6Delam_elem), allocatable :: coh6Delam_lcl
      type(coh8Delam_elem), allocatable :: coh8Delam_lcl
                 
      ! initialize intent out and local variables (non derived type)
      istat = STAT_SUCCESS
      emsg  = ''
      eltype_lcl = ''

      eltype_lcl = adjustl(eltype)  
      
      select case (trim(eltype_lcl))
        
        case ('coh6Delam')
            ! check no. of nodes, exit program if incorrect
            if ( size(connec) /= 6 ) then
              istat = STAT_FAILURE
              emsg  = 'size of connec is not 6 for coh6Delam base element, &
              & set, abstDelam_elem_module'
              return
            end if
            
            ! allocate local coh6Delam base elem
            allocate(coh6Delam_lcl)
            ! call the set procedure of the base element
            call set (coh6Delam_lcl, connec, istat, emsg)
            ! check istat, if istat is failure, clean up and exit program
            if (istat == STAT_FAILURE) then
              deallocate(coh6Delam_lcl)
              return
            end if
            
            ! update intent inout arg. if no error has been encountered

            ! set elem type
            elem%eltype = eltype_lcl
            ! allocate the appropriate base element
            if (.not. allocated(elem%coh6Delam)) allocate(elem%coh6Delam)
            ! deallocate the other base element
            if (allocated(elem%coh8Delam)) deallocate(elem%coh8Delam)
            ! copy definition from local variable
            elem%coh6Delam = coh6Delam_lcl
        
        case ('coh8Delam')
            ! check no. of nodes, exit program if incorrect
            if ( size(connec) /= 8 ) then
              istat = STAT_FAILURE
              emsg  = 'size of connec is not 8 for coh8Delam base element, &
              & set, abstDelam_elem_module'
              return
            end if
            
            ! allocate local coh8Delam base elem
            allocate(coh8Delam_lcl)
            ! call the set procedure of the base element
            call set (coh8Delam_lcl, connec, istat, emsg)
            ! check istat, if istat is failure, clean up and exit program
            if (istat == STAT_FAILURE) then
              deallocate(coh8Delam_lcl)
              return
            end if
            
            ! update intent inout arg. if no error has been encountered

            ! set elem type
            elem%eltype = eltype_lcl
            ! allocate the appropriate base element
            if (.not. allocated(elem%coh8Delam)) allocate(elem%coh8Delam)
            ! deallocate the other base element
            if (allocated(elem%coh6Delam)) deallocate(elem%coh6Delam)
            ! copy definition from local variable
            elem%coh8Delam = coh8Delam_lcl
        
        case default
            ! this should not be reached, flag an error and return
            istat = STAT_FAILURE
            emsg  = 'unsupported eltype in set, baseCoh element module'
            return
      
      end select

  end subroutine set_abstDelam_elem
  
  
  
  pure subroutine extract_abstDelam_elem (elem, eltype, fstat, connec, &
  & ig_points, ig_angles, traction, separation, dm)
  ! extra modules needed to declare the type of some dummy args
  use cohesive_material_module, only : cohesive_ig_point
  use coh6Delam_elem_module,    only : extract
  use coh8Delam_elem_module,    only : extract

    type(abstDelam_elem),                           intent(in)  :: elem
    character(len=ELTYPELENGTH),          optional, intent(out) :: eltype
    integer,                              optional, intent(out) :: fstat
    integer,                 allocatable, optional, intent(out) :: connec(:)
    type(cohesive_ig_point), allocatable, optional, intent(out) :: ig_points(:)
    real(DP),                allocatable, optional, intent(out) :: ig_angles(:)
    real(DP),                             optional, intent(out) :: traction(NST)
    real(DP),                             optional, intent(out) :: separation(NST)
    real(DP),                             optional, intent(out) :: dm

    if (present(eltype)) eltype = elem%eltype

    ! based on eltype, call the respective extract procedure to extract the
    ! requested components

    select case (trim(elem%eltype))

      case ('coh6Delam')
        if (present(fstat))       call extract (elem%coh6Delam, fstat=fstat)
        if (present(connec))      call extract (elem%coh6Delam, connec=connec)
        if (present(ig_points))   call extract (elem%coh6Delam, ig_points=ig_points)
        if (present(ig_angles))   call extract (elem%coh6Delam, ig_angles=ig_angles)
        if (present(traction))    call extract (elem%coh6Delam, traction=traction)
        if (present(separation))  call extract (elem%coh6Delam, separation=separation)
        if (present(dm))          call extract (elem%coh6Delam, dm=dm)

      case ('coh8Delam')
        if (present(fstat))       call extract (elem%coh8Delam, fstat=fstat)
        if (present(connec))      call extract (elem%coh8Delam, connec=connec)
        if (present(ig_points))   call extract (elem%coh8Delam, ig_points=ig_points)
        if (present(ig_angles))   call extract (elem%coh8Delam, ig_angles=ig_angles)
        if (present(traction))    call extract (elem%coh8Delam, traction=traction)
        if (present(separation))  call extract (elem%coh8Delam, separation=separation)
        if (present(dm))          call extract (elem%coh8Delam, dm=dm)

      case default
      ! this should not be reached

    end select


  end subroutine extract_abstDelam_elem



  pure subroutine integrate_abstDelam_elem (elem, nodes, material, theta1, theta2,&
  & Kmatrix, Fvector, istat, emsg, nofailure)
  ! extra modules needed to declare the type of some dummy args
  use fnode_module,             only : fnode
  use cohesive_material_module, only : cohesive_material
  use coh6Delam_elem_module,    only : integrate
  use coh8Delam_elem_module,    only : integrate

      type(abstDelam_elem),     intent(inout) :: elem
      type(fnode),              intent(in)    :: nodes(:)
      type(cohesive_material),  intent(in)    :: material
      real(DP),                 intent(in)    :: theta1, theta2
      real(DP), allocatable,    intent(out)   :: Kmatrix(:,:), Fvector(:)
      integer,                  intent(out)   :: istat
      character(len=MSGLENGTH), intent(out)   :: emsg
      logical,        optional, intent(in)    :: nofailure

      ! local variables
      logical :: nofail

      ! initialize intent out and local variable
      istat  = STAT_SUCCESS
      emsg   = ''
      nofail = .false.

      if (present(nofailure)) nofail = nofailure

      ! call the respective integrate procedure
      ! note that in case of error, the appropriate actions on K and F should
      ! have already been done in the called procedure. nothing need to be done
      ! here.
      select case(elem%eltype)

        case('coh6Delam')

            ! check no. of nodes, exit program if incorrect
            if ( size(nodes) /= 6 ) then
              istat = STAT_FAILURE
              emsg  = 'size of nodes is not 6 for coh6Delam base element, &
              & integrate, abstDelam_elem_module'
              return
            end if
            
            call integrate(elem%coh6Delam, nodes, material, theta1, theta2, &
            & Kmatrix, Fvector, istat, emsg, nofail)

        case('coh8Delam')

            ! check no. of nodes, exit program if incorrect
            if ( size(nodes) /= 8 ) then
              istat = STAT_FAILURE
              emsg  = 'size of nodes is not 8 for coh8Delam base element, &
              & integrate, abstDelam_elem_module'
              return
            end if
            
            call integrate(elem%coh8Delam, nodes, material, theta1, theta2, &
            & Kmatrix, Fvector, istat, emsg, nofail)

        case default
            ! this should not be reached
            istat = STAT_FAILURE
            emsg  = 'unexpected elem type, integrate, abstDelam_elem_module'
            return
            
      end select

  end subroutine integrate_abstDelam_elem


  

end module abstDelam_elem_module
