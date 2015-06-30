module baseCoh_element_module
!
!  Purpose:
!    define an 'abstract element' object to interface with the
!    base delam6 and delam8 elements
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

use delam6_element_module, only : delam6_element
use delam8_element_module, only : delam8_element

  implicit none
  private
  
  type,public :: baseCoh_element
    ! encapsulate components of this type
    private 
    ! list of type components:
    ! eltype  : element type
    ! delam6  : allocated if eltype = 'delam6'
    ! delam8  : allocated if eltype = 'delam8'
    character(len=ELTYPELENGTH)       :: eltype = ''
    type(delam6_element), allocatable :: delam6
    type(delam8_element), allocatable :: delam8   
  end type

  interface empty
      module procedure empty_baseCoh_element
  end interface

  interface set
      module procedure set_baseCoh_element
  end interface
  
  interface integrate
      module procedure integrate_baseCoh_element
  end interface
  
  interface extract
      module procedure extract_baseCoh_element
  end interface




  public :: empty, set, integrate, extract

  
  
  
  contains
  
  
   
  
  pure subroutine empty_baseCoh_element (elem)
  
    type(baseCoh_element), intent(inout) :: elem
    
    elem%eltype=''

    if(allocated(elem%delam6))      deallocate(elem%delam6)
    if(allocated(elem%delam8))      deallocate(elem%delam8)
      
  end subroutine empty_baseCoh_element
   
  
  
  pure subroutine set_baseCoh_element (elem, eltype, connec, istat, emsg)
  ! Purpose:
  ! this is an interface used to set the contained base element according to
  ! eltype:
  ! - if eltype is delam8, the input dummy args are passed to set_delam8_element
  ! to define the elem's delam8 component; elem's delam6 component is dealloc.
  ! - if eltype is delam6, the input dummy args are passed to set_delam6_element
  ! to define the elem's delam6 component; elem's delam8 component is dealloc.
  ! ** note:
  ! - size of connec is checked here against eltype to ensure compatibility
  ! the values of connec and other dummy args are not checked here; they are
  ! left for the called eltype's set procedure for checking
  ! - local copies of elem's components are used for set operation;
  ! they are copied to actual elem's components only before successful return
  use delam6_element_module, only : delam6_element, set
  use delam8_element_module, only : delam8_element, set
  
      type(baseCoh_element),    intent(inout) :: elem
      character(len=*),         intent(in)    :: eltype
      integer,                  intent(in)    :: connec(:)
      integer,                  intent(out)   :: istat
      character(len=MSGLENGTH), intent(out)   :: emsg
      
      ! local variables
      character(len=ELTYPELENGTH)       :: eltype_lcl
      type(delam6_element), allocatable :: delam6_lcl
      type(delam8_element), allocatable :: delam8_lcl
                 
      ! initialize intent out and local variables (non derived type)
      istat = STAT_SUCCESS
      emsg  = ''
      eltype_lcl = ''

      eltype_lcl = adjustl(eltype)  
      
      select case (trim(eltype_lcl))
        
        case ('delam6')
            ! check no. of nodes, exit program if incorrect
            if ( size(connec) /= 6 ) then
              istat = STAT_FAILURE
              emsg  = 'size of connec is not 6 for delam6 base element, &
              & set, baseCoh_element_module'
              return
            end if
            
            ! allocate local delam6 base elem
            allocate(delam6_lcl)
            ! call the set procedure of the base element
            call set (delam6_lcl, connec, istat, emsg)
            ! check istat, if istat is failure, clean up and exit program
            if (istat == STAT_FAILURE) then
              deallocate(delam6_lcl)
              return
            end if
            
            ! update intent inout arg. if no error has been encountered

            ! set elem type
            elem%eltype = eltype_lcl
            ! allocate the appropriate base element
            if (.not. allocated(elem%delam6)) allocate(elem%delam6)
            ! deallocate the other base element
            if (allocated(elem%delam8)) deallocate(elem%delam8)
            ! copy definition from local variable
            elem%delam6 = delam6_lcl
        
        case ('delam8')
            ! check no. of nodes, exit program if incorrect
            if ( size(connec) /= 8 ) then
              istat = STAT_FAILURE
              emsg  = 'size of connec is not 8 for delam8 base element, &
              & set, baseCoh_element_module'
              return
            end if
            
            ! allocate local delam8 base elem
            allocate(delam8_lcl)
            ! call the set procedure of the base element
            call set (delam8_lcl, connec, istat, emsg)
            ! check istat, if istat is failure, clean up and exit program
            if (istat == STAT_FAILURE) then
              deallocate(delam8_lcl)
              return
            end if
            
            ! update intent inout arg. if no error has been encountered

            ! set elem type
            elem%eltype = eltype_lcl
            ! allocate the appropriate base element
            if (.not. allocated(elem%delam8)) allocate(elem%delam8)
            ! deallocate the other base element
            if (allocated(elem%delam6)) deallocate(elem%delam6)
            ! copy definition from local variable
            elem%delam8 = delam8_lcl
        
        case default
            ! this should not be reached, flag an error and return
            istat = STAT_FAILURE
            emsg  = 'unsupported eltype in set, baseCoh element module'
            return
      
      end select

  end subroutine set_baseCoh_element
  
  
  
  pure subroutine extract_baseCoh_element (elem, eltype, fstat, connec, &
  & ig_points, ig_angles, traction, separation, dm)
  ! extra modules needed to declare the type of some dummy args
  use cohesive_material_module, only :cohesive_ig_point
  use delam6_element_module,    only : extract
  use delam8_element_module,    only : extract

    type(baseCoh_element),                          intent(in)  :: elem
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

      case ('delam6')
        if (present(fstat))       call extract (elem%delam6, fstat=fstat)
        if (present(connec))      call extract (elem%delam6, connec=connec)
        if (present(ig_points))   call extract (elem%delam6, ig_points=ig_points)
        if (present(ig_angles))   call extract (elem%delam6, ig_angles=ig_angles)
        if (present(traction))    call extract (elem%delam6, traction=traction)
        if (present(separation))  call extract (elem%delam6, separation=separation)
        if (present(dm))          call extract (elem%delam6, dm=dm)

      case ('delam8')
        if (present(fstat))       call extract (elem%delam8, fstat=fstat)
        if (present(connec))      call extract (elem%delam8, connec=connec)
        if (present(ig_points))   call extract (elem%delam8, ig_points=ig_points)
        if (present(ig_angles))   call extract (elem%delam8, ig_angles=ig_angles)
        if (present(traction))    call extract (elem%delam8, traction=traction)
        if (present(separation))  call extract (elem%delam8, separation=separation)
        if (present(dm))          call extract (elem%delam8, dm=dm)

      case default
      ! this should not be reached

    end select


  end subroutine extract_baseCoh_element



  pure subroutine integrate_baseCoh_element (elem, nodes, material, theta1, theta2,&
  & Kmatrix, Fvector, istat, emsg, nofailure)
  ! extra modules needed to declare the type of some dummy args
  use fnode_module,             only : fnode
  use cohesive_material_module, only : cohesive_material
  use delam6_element_module,    only : integrate
  use delam8_element_module,    only : integrate

      type(baseCoh_element),    intent(inout) :: elem
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

        case('delam6')

            ! check no. of nodes, exit program if incorrect
            if ( size(nodes) /= 6 ) then
              istat = STAT_FAILURE
              emsg  = 'size of nodes is not 6 for delam6 base element, &
              & integrate, baseCoh_element_module'
              return
            end if
            
            call integrate(elem%delam6, nodes, material, theta1, theta2, &
            & Kmatrix, Fvector, istat, emsg, nofail)

        case('delam8')

            ! check no. of nodes, exit program if incorrect
            if ( size(nodes) /= 8 ) then
              istat = STAT_FAILURE
              emsg  = 'size of nodes is not 8 for delam8 base element, &
              & integrate, baseCoh_element_module'
              return
            end if
            
            call integrate(elem%delam8, nodes, material, theta1, theta2, &
            & Kmatrix, Fvector, istat, emsg, nofail)

        case default
            ! this should not be reached
            istat = STAT_FAILURE
            emsg  = 'unexpected elem type, integrate, baseCoh_element_module'
            return
            
      end select

  end subroutine integrate_baseCoh_element


  

end module baseCoh_element_module
