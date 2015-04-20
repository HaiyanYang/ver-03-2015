module basePly_element_module
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

use wedge_element_module  ! use everything available
use brick_element_module  ! use everything available

  implicit none
  private

  type, public :: basePly_element
    ! encapsulate components of this type
    private
    ! list of type components:
    ! eltype : element type
    ! wedge  : allocated if eltype = 'wedge'
    ! brick  : allocated if eltype = 'brick'
    character(len=ELTYPELENGTH)      :: eltype = ''
    type(wedge_element), allocatable :: wedge      ! wedge sub elements
    type(brick_element), allocatable :: brick      ! brick sub elements
  end type

  interface empty
      module procedure empty_basePly_element
  end interface

  interface set
      module procedure set_basePly_element
  end interface

  interface integrate
      module procedure integrate_basePly_element
  end interface

  interface extract
      module procedure extract_basePly_element
  end interface




  public :: empty, set, integrate, extract




  contains




  pure subroutine empty_basePly_element (elem)

      type(basePly_element), intent(inout) :: elem

      elem%eltype = ''

      if(allocated(elem%wedge)) deallocate(elem%wedge)
      if(allocated(elem%brick)) deallocate(elem%brick)

  end subroutine empty_basePly_element



  pure subroutine set_basePly_element (elem, eltype, connec, ID_matlist, &
  & ply_angle, istat, emsg)
  ! Purpose:
  ! this is an interface used to set the contained base element according to
  ! eltype:
  ! - if eltype is brick, the input dummy args are passed to set_brick_element
  ! to define the elem's brick component; elem's wedge component is dealloc.
  ! - if eltype is wedge, the input dummy args are passed to set_wedge_element
  ! to define the elem's wedge component; elem's brick component is dealloc.
  ! ** note:
  ! - size of connec is checked here against eltype to ensure compatibility
  ! the values of connec and other dummy args are not checked here; they are
  ! left for the called eltype's set procedure for checking
  ! - local copies of elem's components are used for set operation;
  ! they are copied to actual elem's components only before successful return

    type(basePly_element),       intent(inout) :: elem
    character(len=*),            intent(in)    :: eltype
    integer,                     intent(in)    :: connec(:)
    integer,                     intent(in)    :: ID_matlist
    real(DP),                    intent(in)    :: ply_angle
    integer,                     intent(out)   :: istat
    character(len=MSGLENGTH),    intent(out)   :: emsg

    ! local var., local copies of intent inout arg. components
    character(len=ELTYPELENGTH) :: eltype_lcl
    type(wedge_element), allocatable :: wedge_lcl
    type(brick_element), allocatable :: brick_lcl

    ! initialize intent out and local variables (non derived type)
    istat = STAT_SUCCESS
    emsg  = ''
    eltype_lcl = ''

    eltype_lcl = adjustl(eltype)

    select case (trim(eltype_lcl))

      case ('wedge')
          ! check no. of nodes, exit program if incorrect
          if ( size(connec) /= 6 ) then
            istat = STAT_FAILURE
            emsg  = 'size of connec is not 6 for wedge base element, &
            & set, basePly_element_module'
            return
          end if

          ! allocate local wedge base elem
          allocate(wedge_lcl)
          ! call the set procedure of the base element
          call set (wedge_lcl, connec, ID_matlist, ply_angle, istat, emsg)
          ! check istat, if istat is failure, clean up and exit program
          if (istat == STAT_FAILURE) then
            deallocate(wedge_lcl)
            return
          end if

          ! update intent inout arg. if no error has been encountered

          ! set elem type
          elem%eltype = eltype_lcl
          ! allocate the appropriate base element
          if (.not. allocated(elem%wedge)) allocate(elem%wedge)
          ! deallocate the other base element
          if (allocated(elem%brick)) deallocate(elem%brick)
          ! copy definition from local variable
          elem%wedge = wedge_lcl


      case ('brick')
          ! check no. of nodes, exit program if incorrect
          if ( size(connec) /= 8 ) then
            istat = STAT_FAILURE
            emsg  = 'size of connec is not 8 for brick base element, &
            & set, basePly_element_module'
            return
          end if

          ! allocate local brick base elem
          allocate(brick_lcl)
          ! call the set procedure of the base element
          call set (brick_lcl, connec, ID_matlist, ply_angle, istat, emsg)
          ! check istat, if istat is failure, clean up and exit program
          if (istat == STAT_FAILURE) then
            deallocate(brick_lcl)
            return
          end if

          ! update intent inout arg. if no error has been encountered

          ! set elem type
          elem%eltype = eltype_lcl
          ! allocate the appropriate base element
          if (.not. allocated(elem%brick)) allocate(elem%brick)
          ! deallocate the other base element
          if (allocated(elem%wedge)) deallocate(elem%wedge)
          ! copy definition from local variable
          elem%brick = brick_lcl


      case default
          ! this should not be reached, flag an error and return
          istat = STAT_FAILURE
          emsg  = 'unsupported eltype in set, basePly element module'
          return


    end select

    if (allocated(wedge_lcl)) deallocate(wedge_lcl)
    if (allocated(brick_lcl)) deallocate(brick_lcl)

  end subroutine set_basePly_element



  pure subroutine extract_basePly_element (elem, eltype, fstat, connec, &
  & ID_matlist, ply_angle, local_clock, ig_points, stress, strain, df)
  ! extra modules needed to declare the type of some dummy args
  use global_clock_module
  use lamina_material_module

    type(basePly_element),                        intent(in)  :: elem
    character(len=ELTYPELENGTH),        optional, intent(out) :: eltype
    integer,                            optional, intent(out) :: fstat
    integer,               allocatable, optional, intent(out) :: connec(:)
    integer,                            optional, intent(out) :: ID_matlist
    real(DP),                           optional, intent(out) :: ply_angle
    type(program_clock),                optional, intent(out) :: local_clock
    type(lamina_ig_point), allocatable, optional, intent(out) :: ig_points(:)
    real(DP),                           optional, intent(out) :: stress(NST)
    real(DP),                           optional, intent(out) :: strain(NST)
    real(DP),                           optional, intent(out) :: df

    if (present(eltype)) eltype = elem%eltype

    ! based on eltype, call the respective extract procedure to extract the
    ! requested components

    select case (trim(elem%eltype))

      case ('wedge')
        if (present(fstat))       call extract (elem%wedge, fstat=fstat)
        if (present(connec))      call extract (elem%wedge, connec=connec)
        if (present(ID_matlist))  call extract (elem%wedge, ID_matlist=ID_matlist)
        if (present(ply_angle))   call extract (elem%wedge, ply_angle=ply_angle)
        if (present(local_clock)) call extract (elem%wedge, local_clock=local_clock)
        if (present(ig_points))   call extract (elem%wedge, ig_points=ig_points)
        if (present(stress))      call extract (elem%wedge, stress=stress)
        if (present(strain))      call extract (elem%wedge, strain=strain)
        if (present(df))          call extract (elem%wedge, df=df)

      case ('brick')
        if (present(fstat))       call extract (elem%brick, fstat=fstat)
        if (present(connec))      call extract (elem%brick, connec=connec)
        if (present(ID_matlist))  call extract (elem%brick, ID_matlist=ID_matlist)
        if (present(ply_angle))   call extract (elem%brick, ply_angle=ply_angle)
        if (present(local_clock)) call extract (elem%brick, local_clock=local_clock)
        if (present(ig_points))   call extract (elem%brick, ig_points=ig_points)
        if (present(stress))      call extract (elem%brick, stress=stress)
        if (present(strain))      call extract (elem%brick, strain=strain)
        if (present(df))          call extract (elem%brick, df=df)

      case default
      ! this should not be reached

    end select


  end subroutine extract_basePly_element



  pure subroutine integrate_basePly_element (elem, Kmatrix, Fvector, istat, &
  & emsg, nofailure, mnodes)
  use xnode_module

      type(basePly_element),    intent(inout) :: elem
      real(DP), allocatable,    intent(out)   :: Kmatrix(:,:), Fvector(:)
      integer,                  intent(out)   :: istat
      character(len=MSGLENGTH), intent(out)   :: emsg
      logical,        optional, intent(in)    :: nofailure
      type(xnode),    optional, intent(in)    :: mnodes(:)

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

          case('wedge')
            if (present(mnodes)) then
              call integrate(elem%wedge, Kmatrix, Fvector, istat, emsg,&
              & nofail, mnodes)
            else
              call integrate(elem%wedge, Kmatrix, Fvector, istat, emsg,&
              & nofail)
            end if

          case('brick')
            if (present(mnodes)) then
              call integrate(elem%brick, Kmatrix, Fvector, istat, emsg,&
              & nofail, mnodes)
            else
              call integrate(elem%brick, Kmatrix, Fvector, istat, emsg,&
              & nofail)
            end if

          case default
              ! this should not be reached
      end select

  end subroutine integrate_basePly_element




end module basePly_element_module
