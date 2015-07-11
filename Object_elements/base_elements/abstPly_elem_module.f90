module abstPly_elem_module
!
!  Purpose:
!    define an 'abstract element' object to interface with the 
!    base brick and wedge elements
!    with the associated procedures to empty, set, integrate and extract
!    its components
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    17/04/15  B. Y. Chen            Original code
!
use parameter_module, only : NST => NST_STANDARD, DP, ELTYPELENGTH, &
                      & MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

use wedgePly_elem_module, only : wedgePly_elem
use brickPly_elem_module, only : brickPly_elem

  implicit none
  private

  type, public :: abstPly_elem
    ! encapsulate components of this type
    private
    ! list of type components:
    ! eltype : element type
    ! wedge  : allocated if eltype = 'wedge'
    ! brick  : allocated if eltype = 'brick'
    character(len=ELTYPELENGTH)      :: eltype = ''
    type(wedgePly_elem), allocatable :: wedge      ! wedge sub elements
    type(brickPly_elem), allocatable :: brick      ! brick sub elements
  end type

  interface integrate
      module procedure integrate_abstPly_elem
  end interface

  interface extract
      module procedure extract_abstPly_elem
  end interface




  public :: integrate, extract




  contains
  


  pure subroutine extract_abstPly_elem (elem, eltype, fstat, ig_points, stress, strain, df)
  ! extra modules needed to declare the type of some dummy args
  use lamina_material_module, only : lamina_ig_point
  use wedgePly_elem_module,   only : extract
  use brickPly_elem_module,   only : extract

    type(abstPly_elem),                           intent(in)  :: elem
    character(len=ELTYPELENGTH),        optional, intent(out) :: eltype
    integer,                            optional, intent(out) :: fstat
    type(lamina_ig_point), allocatable, optional, intent(out) :: ig_points(:)
    real(DP),                           optional, intent(out) :: stress(NST)
    real(DP),                           optional, intent(out) :: strain(NST)
    real(DP),                           optional, intent(out) :: df

    if (present(eltype)) eltype = elem%eltype

    ! based on eltype, call the respective extract procedure to extract the
    ! requested components

    select case (trim(adjustl(elem%eltype)))

      case ('wedge')
        if (present(fstat))       call extract (elem%wedge, fstat=fstat)
        if (present(ig_points))   call extract (elem%wedge, ig_points=ig_points)
        if (present(stress))      call extract (elem%wedge, stress=stress)
        if (present(strain))      call extract (elem%wedge, strain=strain)
        if (present(df))          call extract (elem%wedge, df=df)

      case ('brick')
        if (present(fstat))       call extract (elem%brick, fstat=fstat)
        if (present(ig_points))   call extract (elem%brick, ig_points=ig_points)
        if (present(stress))      call extract (elem%brick, stress=stress)
        if (present(strain))      call extract (elem%brick, strain=strain)
        if (present(df))          call extract (elem%brick, df=df)

      case default
      ! this should not be reached

    end select


  end subroutine extract_abstPly_elem



  pure subroutine integrate_abstPly_elem (elem, nodes, ply_angle, material, Kmatrix, &
  & Fvector, istat, emsg, nofailure)
  ! extra modules needed to declare the type of some dummy args
  use fnode_module,           only : fnode
  use lamina_material_module, only : lamina_material
  use wedgePly_elem_module,   only : integrate
  use brickPly_elem_module,   only : integrate

      type(abstPly_elem),       intent(inout) :: elem
      type(fnode),              intent(in)    :: nodes(:)
      real(DP),                 intent(in)    :: ply_angle
      type(lamina_material),    intent(in)    :: material
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
            
            elem%eltype = 'wedge'
            
            if (.not. allocated(elem%wedge)) allocate(elem%wedge)
            if (allocated(elem%brick))     deallocate(elem%brick)
            
            call integrate(elem%wedge, nodes, ply_angle, material, Kmatrix, Fvector, &
            & istat, emsg, nofail)

        case(8)
        
            elem%eltype = 'brick'
            
            if (.not. allocated(elem%brick)) allocate(elem%brick)
            if (allocated(elem%wedge))     deallocate(elem%wedge)
            
            call integrate(elem%brick, nodes, ply_angle, material, Kmatrix, Fvector, &
            & istat, emsg, nofail)

        case default
            istat = STAT_FAILURE
            emsg  = 'unexpected no. of nodes, integrate, abstPly_elem_module'
            return
            
      end select

  end subroutine integrate_abstPly_elem




end module abstPly_elem_module
