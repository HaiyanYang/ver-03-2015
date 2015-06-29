module baseCoh_element_module
!
!  Purpose:
!    define an 'abstract element' object to interface with the
!    base coh3d6 and coh3d8 elements
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

use coh3d6_element_module, only : coh3d6_element
use coh3d8_element_module, only : coh3d8_element

  implicit none
  private
  
  type,public :: baseCoh_element
    ! encapsulate components of this type
    private 
    ! list of type components:
    ! eltype  : element type
    ! coh3d6  : allocated if eltype = 'coh3d6'
    ! coh3d8  : allocated if eltype = 'coh3d8'
    character(len=ELTYPELENGTH)       :: eltype = ''
    type(coh3d6_element), allocatable :: coh3d6
    type(coh3d8_element), allocatable :: coh3d8   
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

    if(allocated(elem%coh3d6))      deallocate(elem%coh3d6)
    if(allocated(elem%coh3d8))      deallocate(elem%coh3d8)
      
  end subroutine empty_baseCoh_element
   
  
  
  pure subroutine set_baseCoh_element (elem, eltype, connec, istat, emsg)
  ! Purpose:
  ! this is an interface used to set the contained base element according to
  ! eltype:
  ! - if eltype is coh3d8, the input dummy args are passed to set_coh3d8_element
  ! to define the elem's coh3d8 component; elem's coh3d6 component is dealloc.
  ! - if eltype is coh3d6, the input dummy args are passed to set_coh3d6_element
  ! to define the elem's coh3d6 component; elem's coh3d8 component is dealloc.
  ! ** note:
  ! - size of connec is checked here against eltype to ensure compatibility
  ! the values of connec and other dummy args are not checked here; they are
  ! left for the called eltype's set procedure for checking
  ! - local copies of elem's components are used for set operation;
  ! they are copied to actual elem's components only before successful return
  use coh3d6_element_module, only : coh3d6_element, set
  use coh3d8_element_module, only : coh3d8_element, set
  
      type(baseCoh_element),    intent(inout) :: elem
      character(len=*),         intent(in)    :: eltype
      integer,                  intent(in)    :: connec(:)
      integer,                  intent(out)   :: istat
      character(len=MSGLENGTH), intent(out)   :: emsg
      
      ! local variables
      character(len=ELTYPELENGTH)       :: eltype_lcl
      type(coh3d6_element), allocatable :: coh3d6_lcl
      type(coh3d8_element), allocatable :: coh3d8_lcl
                 
      ! initialize intent out and local variables (non derived type)
      istat = STAT_SUCCESS
      emsg  = ''
      eltype_lcl = ''

      eltype_lcl = adjustl(eltype)  
      
      select case (trim(eltype_lcl))
        
        case ('coh3d6')
            ! check no. of nodes, exit program if incorrect
            if ( size(connec) /= 6 ) then
              istat = STAT_FAILURE
              emsg  = 'size of connec is not 6 for coh3d6 base element, &
              & set, baseCoh_element_module'
              return
            end if
            
            ! allocate local coh3d6 base elem
            allocate(coh3d6_lcl)
            ! call the set procedure of the base element
            call set (coh3d6_lcl, connec, istat, emsg)
            ! check istat, if istat is failure, clean up and exit program
            if (istat == STAT_FAILURE) then
              deallocate(coh3d6_lcl)
              return
            end if
            
            ! update intent inout arg. if no error has been encountered

            ! set elem type
            elem%eltype = eltype_lcl
            ! allocate the appropriate base element
            if (.not. allocated(elem%coh3d6)) allocate(elem%coh3d6)
            ! deallocate the other base element
            if (allocated(elem%coh3d8)) deallocate(elem%coh3d8)
            ! copy definition from local variable
            elem%coh3d6 = coh3d6_lcl
        
        case ('coh3d8')
            ! check no. of nodes, exit program if incorrect
            if ( size(connec) /= 8 ) then
              istat = STAT_FAILURE
              emsg  = 'size of connec is not 8 for coh3d8 base element, &
              & set, baseCoh_element_module'
              return
            end if
            
            ! allocate local coh3d8 base elem
            allocate(coh3d8_lcl)
            ! call the set procedure of the base element
            call set (coh3d8_lcl, connec, istat, emsg)
            ! check istat, if istat is failure, clean up and exit program
            if (istat == STAT_FAILURE) then
              deallocate(coh3d8_lcl)
              return
            end if
            
            ! update intent inout arg. if no error has been encountered

            ! set elem type
            elem%eltype = eltype_lcl
            ! allocate the appropriate base element
            if (.not. allocated(elem%coh3d8)) allocate(elem%coh3d8)
            ! deallocate the other base element
            if (allocated(elem%coh3d6)) deallocate(elem%coh3d6)
            ! copy definition from local variable
            elem%coh3d8 = coh3d8_lcl
        
        case default
            ! this should not be reached, flag an error and return
            istat = STAT_FAILURE
            emsg  = 'unsupported eltype in set, baseCoh element module'
            return
      
      end select

  end subroutine set_baseCoh_element
  
  
  
  pure subroutine extract_baseCoh_element (elem, eltype, fstat, connec, &
  & ig_points, traction, separation, dm)
  ! extra modules needed to declare the type of some dummy args
  use cohesive_material_module, only :cohesive_ig_point
  use coh3d6_element_module,    only : extract
  use coh3d8_element_module,    only : extract

    type(baseCoh_element),                          intent(in)  :: elem
    character(len=ELTYPELENGTH),          optional, intent(out) :: eltype
    integer,                              optional, intent(out) :: fstat
    integer,                 allocatable, optional, intent(out) :: connec(:)
    type(cohesive_ig_point), allocatable, optional, intent(out) :: ig_points(:)
    real(DP),                             optional, intent(out) :: traction(NST)
    real(DP),                             optional, intent(out) :: separation(NST)
    real(DP),                             optional, intent(out) :: dm

    if (present(eltype)) eltype = elem%eltype

    ! based on eltype, call the respective extract procedure to extract the
    ! requested components

    select case (trim(elem%eltype))

      case ('coh3d6')
        if (present(fstat))       call extract (elem%coh3d6, fstat=fstat)
        if (present(connec))      call extract (elem%coh3d6, connec=connec)
        if (present(ig_points))   call extract (elem%coh3d6, ig_points=ig_points)
        if (present(traction))    call extract (elem%coh3d6, traction=traction)
        if (present(separation))  call extract (elem%coh3d6, separation=separation)
        if (present(dm))          call extract (elem%coh3d6, dm=dm)

      case ('coh3d8')
        if (present(fstat))       call extract (elem%coh3d8, fstat=fstat)
        if (present(connec))      call extract (elem%coh3d8, connec=connec)
        if (present(ig_points))   call extract (elem%coh3d8, ig_points=ig_points)
        if (present(traction))    call extract (elem%coh3d8, traction=traction)
        if (present(separation))  call extract (elem%coh3d8, separation=separation)
        if (present(dm))          call extract (elem%coh3d8, dm=dm)

      case default
      ! this should not be reached

    end select


  end subroutine extract_baseCoh_element



  pure subroutine integrate_baseCoh_element (elem, nodes, material, &
  & Kmatrix, Fvector, istat, emsg, nofailure)
  ! extra modules needed to declare the type of some dummy args
  use xnode_module,             only : xnode
  use cohesive_material_module, only : cohesive_material
  use coh3d6_element_module,    only : integrate
  use coh3d8_element_module,    only : integrate

      type(baseCoh_element),    intent(inout) :: elem
      type(xnode),              intent(in)    :: nodes(:)
      type(cohesive_material),  intent(in)    :: material
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

        case('coh3d6')

            ! check no. of nodes, exit program if incorrect
            if ( size(nodes) /= 6 ) then
              istat = STAT_FAILURE
              emsg  = 'size of nodes is not 6 for coh3d6 base element, &
              & integrate, baseCoh_element_module'
              return
            end if
            
            call integrate(elem%coh3d6, nodes, material, Kmatrix, Fvector, &
            & istat, emsg, nofail)

        case('coh3d8')

            ! check no. of nodes, exit program if incorrect
            if ( size(nodes) /= 8 ) then
              istat = STAT_FAILURE
              emsg  = 'size of nodes is not 8 for coh3d8 base element, &
              & integrate, baseCoh_element_module'
              return
            end if
            
            call integrate(elem%coh3d8, nodes, material, Kmatrix, Fvector, &
            & istat, emsg, nofail)

        case default
            ! this should not be reached
            istat = STAT_FAILURE
            emsg  = 'unexpected elem type, integrate, baseCoh_element_module'
            return
            
      end select

  end subroutine integrate_baseCoh_element


  

end module baseCoh_element_module
