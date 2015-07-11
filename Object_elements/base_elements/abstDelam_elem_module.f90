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
  
  interface integrate
      module procedure integrate_abstDelam_elem
  end interface
  
  interface extract
      module procedure extract_abstDelam_elem
  end interface




  public :: integrate, extract

  
  
  
  contains
  
 
  
  pure subroutine extract_abstDelam_elem (elem, eltype, fstat, &
  & ig_points, ig_angles, traction, separation, dm)
  ! extra modules needed to declare the type of some dummy args
  use cohesive_material_module, only : cohesive_ig_point
  use coh6Delam_elem_module,    only : extract
  use coh8Delam_elem_module,    only : extract

    type(abstDelam_elem),                           intent(in)  :: elem
    character(len=ELTYPELENGTH),          optional, intent(out) :: eltype
    integer,                              optional, intent(out) :: fstat
    type(cohesive_ig_point), allocatable, optional, intent(out) :: ig_points(:)
    real(DP),                allocatable, optional, intent(out) :: ig_angles(:)
    real(DP),                             optional, intent(out) :: traction(NST)
    real(DP),                             optional, intent(out) :: separation(NST)
    real(DP),                             optional, intent(out) :: dm

    if (present(eltype)) eltype = elem%eltype

    ! based on eltype, call the respective extract procedure to extract the
    ! requested components

    select case (trim(adjustl(elem%eltype)))

      case ('coh6Delam')
        if (present(fstat))       call extract (elem%coh6Delam, fstat=fstat)
        if (present(ig_points))   call extract (elem%coh6Delam, ig_points=ig_points)
        if (present(ig_angles))   call extract (elem%coh6Delam, ig_angles=ig_angles)
        if (present(traction))    call extract (elem%coh6Delam, traction=traction)
        if (present(separation))  call extract (elem%coh6Delam, separation=separation)
        if (present(dm))          call extract (elem%coh6Delam, dm=dm)

      case ('coh8Delam')
        if (present(fstat))       call extract (elem%coh8Delam, fstat=fstat)
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
      select case(size(nodes))

        case(6)
            
            elem%eltype = 'coh6Delam'
            
            if (.not. allocated(elem%coh6Delam)) allocate(elem%coh6Delam)
            if (allocated(elem%coh8Delam))     deallocate(elem%coh8Delam)
            
            call integrate(elem%coh6Delam, nodes, material, theta1, theta2, &
            & Kmatrix, Fvector, istat, emsg, nofail)

        case(8)
        
            elem%eltype = 'coh8Delam'
            
            if (.not. allocated(elem%coh8Delam)) allocate(elem%coh8Delam)
            if (allocated(elem%coh6Delam))     deallocate(elem%coh6Delam)
            
            call integrate(elem%coh8Delam, nodes, material, theta1, theta2, &
            & Kmatrix, Fvector, istat, emsg, nofail)

        case default
            ! this should not be reached
            istat = STAT_FAILURE
            emsg  = 'unexpected no. of nodes, integrate, abstDelam_elem_module'
            return
            
      end select

  end subroutine integrate_abstDelam_elem


  

end module abstDelam_elem_module
