module fCoh3d8_element_module
!
!  Purpose:
!    define a floating-node coh3d8 element object, with breakable top & bottom
!    surfaces.
!
!  topological definition of this element:
!
!  8__22___(31)___21___7
!  |\                  |\
!  | \23               | \20
!  |  \(32)            |  \(30)     Top    edges (anti-clock wise from front):
!  |   \24             |   \19      topE1, topE2, topE3, topE4
!  |    \              |    \
!  |     \5___17__(29)_|_18__\6
!  |______|____________|      |
! 4\   14 | (27)  13  3\      |
!   \     |             \     |
!  15\    |            12\    |
! (28)\   |           (26)\   |     Bottom edges (anti-clock wise from front):
!    16\  |              11\  |     botE1, botE2, botE3, botE4
!       \ |                 \ |
!        \|__________________\|
!         1    9  (25)  10    2
!
!  :::: bottom surface definition ::::
!  real     nodes on bot edges :  1,  2,  3,  4
!  floating nodes on bot edges :  9, 10, 11, 12, 13, 14, 15, 16
!  internal nodes on bot edges : 25, 26, 27, 28
!  nodes of botE1 ::    <end nodes: 1, 2>  <fl. nodes:  9, 10>  <in. node: 25>
!  nodes of botE2 ::    <end nodes: 2, 3>  <fl. nodes: 11, 12>  <in. node: 26>
!  nodes of botE3 ::    <end nodes: 3, 4>  <fl. nodes: 13, 14>  <in. node: 27>
!  nodes of botE4 ::    <end nodes: 4, 1>  <fl. nodes: 15, 16>  <in. node: 28>
!
!  :::: top surface definition ::::
!  real     nodes on top edges :  5,  6,  7,  8
!  floating nodes on top edges : 17, 18, 19, 20, 21, 22, 23, 24
!  internal nodes on top edges : 29, 30, 31, 32
!  nodes of topE1 ::    <end nodes: 5, 6>  <fl. nodes: 17, 18>  <in. node: 29>
!  nodes of topE2 ::    <end nodes: 6, 7>  <fl. nodes: 19, 20>  <in. node: 30>
!  nodes of topE3 ::    <end nodes: 7, 8>  <fl. nodes: 21, 22>  <in. node: 31>
!  nodes of topE4 ::    <end nodes: 8, 5>  <fl. nodes: 23, 24>  <in. node: 32>
!
!
!
!  topological definition of intact element, type coh3d8_element:
!
!  8___________________7
!  |\                  |\
!  | \                 | \
!  |  \                |  \
!  |   \               |   \
!  |    \              |    \
!  |     \5____________|_____\6
!  |______|____________|      |
! 4\      |           3\      |
!   \     |             \     |
!    \    |              \    |
!     \   |               \   |
!      \  |                \  |
!       \ |                 \ |
!        \|__________________\|
!         1                   2
!
!
!
!  topological definition of top sub element, type fCoh3d8_subelem:
!
!  8____22__TE3___21___7
!  |\                  |\
!  | \23               | \20
!  |  \TE4             |  \TE2      top subelem top edges in this elem top edges
!  |   \24             |   \19      (anti-clock wise from front):
!  |    \              |    \       topE1, topE2, topE3, topE4
!  |     \5___17___TE1_|_18__\6
!  |______|____________|      |
! 4\      |  27       3\      |
!   \     |             \     |
!    \    |              \    |
!   28\   |             26\   |
!      \  |                \  |
!       \ |                 \ |
!        \|__________________\|
!         1        25         2
!
!
!
!  topological definition of bot sub element (flip over), type fCoh3d8_subelem:
!
!  1_____9___BE1__10___2
!  |\                  |\
!  | \16               | \11
!  |  \BE4             |  \BE2      bot subelem top edges in this elem bot edges
!  |   \15             |   \12      (anti-clock wise from front):
!  |    \              |    \       botE3, botE2, botE1, botE4
!  |     \4___14___BE3_|_13__\3
!  |______|____________|      |
! 5\      |  29       6\      |
!   \     |             \     |
!    \    |              \    |
!   32\   |             30\   |
!      \  |                \  |
!       \ |                 \ |
!        \|__________________\|
!         8        31         7
!
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    23/04/15  B. Y. Chen            Original code
!

use parameter_module,       only : NDIM
use coh3d8_element_module,  only : coh3d8_element
use fCoh3d8_subelem_module, only : fCoh3d8_subelem


implicit none

private

! NODE AND EDGE PARAMETERS OF THIS ELEMENT
integer, parameter :: NNDRL       = 8,                      &
                   &  NEDGE       = 8,                      &
                   &  NNDFL       = 2 * NEDGE,              &
                   &  NNDIN       = 1 * NEDGE,              &
                   &  NNODE       = NNDRL + NNDFL + NNDIN,  &
                   &  NDOF        = NDIM * NNODE

! NODAL CONNEC OF INTACT  ELEM: REAL NODES ONLY
integer, parameter :: INTACT_ELEM_NODES(8) = [1,2,3,4,5,6,7,8]

! NODAL CONNEC OF TOP SUB ELEM: REAL NODES, FLOATING NODES AND INTERNAL NODES
integer, parameter :: TOP_SUBELEM_NODES(20) = &
& [1,2,3,4,5,6,7,8,17,18,19,20,21,22,23,24,25,26,27,28]

! NODAL CONNEC OF BOT SUB ELEM: REAL NODES, FLOATING NODES AND INTERNAL NODES
integer, parameter :: BOT_SUBELEM_NODES(20) = [8,7,6,5,4,3,2,1,  &
                   &  14,13,12,11,10, 9,16,15,   31,30,29,32]

type, public :: fCoh3d8_element
    private

    integer :: node_connec(NNODE) = 0
    logical :: top_subelem_set    = .false.
    logical :: bot_subelem_set    = .false.
    type(coh3d8_element),  allocatable   :: intact_elem
    type(fCoh3d8_subelem), allocatable   :: top_subelem
    type(fCoh3d8_subelem), allocatable   :: bot_subelem

end type fCoh3d8_element

interface empty
    module procedure empty_fCoh3d8_element
end interface

interface set
    module procedure set_fCoh3d8_element
end interface

interface update
    module procedure update_fCoh3d8_element
end interface

interface integrate
    module procedure integrate_fCoh3d8_element
end interface

interface extract
    module procedure extract_fCoh3d8_element
end interface


public :: empty, set, update, integrate, extract



contains



pure subroutine empty_fCoh3d8_element (elem)

  type(fCoh3d8_element), intent(inout) :: elem

  type(fCoh3d8_element) :: elem_lcl

  elem = elem_lcl

end subroutine empty_fCoh3d8_element



pure subroutine set_fCoh3d8_element (elem, node_connec, istat, emsg)
! Purpose:
! to set the element ready for first use
use parameter_module,       only : MSGLENGTH, STAT_FAILURE, STAT_SUCCESS
use coh3d8_element_module,  only : set

  type(fCoh3d8_element),    intent(inout) :: elem
  integer,                  intent(in)    :: node_connec(nnode)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local copy of elem
  type(fCoh3d8_element) :: elem_lcl
  ! global_connec of sub element
  integer, allocatable  :: global_connec(:)
  ! location for emsg
  character(len=MSGLENGTH) :: msgloc

  istat = STAT_SUCCESS
  emsg  = ''
  msgloc = ' set, fCoh3d8_element module'

  ! check validity of inputs
  if ( any(node_connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'node connec indices must be >=1'//trim(msgloc)
    return
  end if

  ! update to elem_lcl first
  elem_lcl%node_connec = node_connec

  ! allocate intact elem
  allocate(elem_lcl%intact_elem)
  allocate(global_connec(NNDRL))

  ! populate the global connec of intact element
  global_connec(:) = node_connec(INTACT_ELEM_NODES(:))

  ! set the intact element
  call set (elem_lcl%intact_elem, connec=global_connec, istat=istat, emsg=emsg)

  ! if an error is encountered in set, clean up and exit program
  if (istat == STAT_FAILURE) then
    if (allocated(global_connec)) deallocate(global_connec)
    emsg = emsg//trim(msgloc)
    return
  end if

  ! update to dummy arg. elem before successful return
  elem = elem_lcl

  if (allocated(global_connec)) deallocate(global_connec)

end subroutine set_fCoh3d8_element



pure subroutine update_fCoh3d8_element (elem, ply_edge_status, top_or_bottom, &
& istat, emsg)
! Purpose:
! this subroutine is used to update the edge_status array of the element
! passed in from the adjacent ply elements, and then alloc and set the subelem
use parameter_module, only : MSGLENGTH, STAT_FAILURE, STAT_SUCCESS,&
                      & INTACT, TRANSITION_EDGE, REFINEMENT_EDGE,  &
                      & CRACK_TIP_EDGE, WEAK_CRACK_EDGE,           &
                      & COH_CRACK_EDGE, STRONG_CRACK_EDGE
use fCoh3d8_subelem_module, only : set

  type(fCoh3d8_element),    intent(inout) :: elem
  integer,                  intent(in)    :: ply_edge_status(NEDGE/2)
  character(len=*),         intent(in)    :: top_or_bottom
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local copy of elem
  type(fCoh3d8_element)    :: el
  character(len=MSGLENGTH) :: msgloc
  integer                  :: n_crackedges

  ! if both top and bot sub elems have already been set, return directly
  if (elem%top_subelem_set .and. elem%bot_subelem_set) return

  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ' update, fCoh3d8_element module'
  n_crackedges = 0

  ! check edge status, see if there's any unexpected edge status value
  if ( any( .not. ( ply_edge_status == INTACT          .or.         &
  &                 ply_edge_status == TRANSITION_EDGE .or.         &
  &                 ply_edge_status == REFINEMENT_EDGE .or.         &
  &                 ply_edge_status == CRACK_TIP_EDGE  .or.         &
  &                 ply_edge_status == WEAK_CRACK_EDGE .or.         &
  &                 ply_edge_status == COH_CRACK_EDGE  .or.         &
  &                 ply_edge_status == STRONG_CRACK_EDGE )  )  ) then
    istat = STAT_FAILURE
    emsg  = 'ply edge status value is NOT recognized,'//trim(msgloc)
    return
  end if

  ! check the no. of broken edges; only accepts TWO cracked edges, as this is
  ! the final partition from the ply element
  n_crackedges = count (ply_edge_status >= COH_CRACK_EDGE)
  if (n_crackedges /= 2) then
    istat = STAT_FAILURE
    emsg  = 'no. of cracked edges must be TWO,'//trim(msgloc)
    return
  end if

  ! copy elem to its local copy
  el = elem

  ! update to elem component if checkings are passed
  select case (trim(adjustl(top_or_bottom)))

    case ('top')
        ! return if top subelem is already set
        if (el%top_subelem_set) return

        ! deallocate intact elem
        if (allocated(el%intact_elem)) deallocate(el%intact_elem)

        ! allocate top sub elem
        allocate(el%top_subelem)

        ! set top sub elem; note that the top sub elem's top edges are just
        ! this elem's top edges, without any change of order
        call set (el%top_subelem, node_connec=el%node_connec(TOP_SUBELEM_NODES),&
        & top_edge_status=ply_edge_status, istat=istat, emsg=emsg)
        if (istat == STAT_FAILURE) then
          emsg = emsg//trim(msgloc)
          return
        end if

        el%top_subelem_set = .true.

    case ('bottom')
        ! return if bot subelem is already set
        if (el%bot_subelem_set) return

        ! deallocate intact elem
        if (allocated(el%intact_elem)) deallocate(el%intact_elem)

        ! allocate top sub elem
        allocate(el%bot_subelem)

        ! set top sub elem; note that the bot sub elem's top edges are
        ! this elem's bot edges but with permutated order (see illustration at top
        ! of this module, in comment)
        call set (el%bot_subelem, node_connec=el%node_connec(BOT_SUBELEM_NODES),&
        & top_edge_status=ply_edge_status([3, 2, 1, 4]), istat=istat, &
        & emsg=emsg)
        if (istat == STAT_FAILURE) then
          emsg = emsg//trim(msgloc)
          return
        end if

        el%bot_subelem_set = .true.

    case default
        istat = STAT_FAILURE
        emsg  = "unsupported top_or_bottom value, input either 'top' or &
        &'bottom',"//trim(msgloc)
        return

  end select

  ! copy definition to input arg. before successful return
  elem = el

end subroutine update_fCoh3d8_element



pure subroutine extract_fCoh3d8_element(elem, top_subelem_set, bot_subelem_set,&
& intact_elem, top_subelem, bot_subelem)
use coh3d8_element_module,  only : coh3d8_element
use fCoh3d8_subelem_module, only : fCoh3d8_subelem

  type(fCoh3d8_element),                        intent(in)   :: elem
  logical,                            optional, intent(out)  :: top_subelem_set
  logical,                            optional, intent(out)  :: bot_subelem_set
  type(coh3d8_element),  allocatable, optional, intent(out)  :: intact_elem
  type(fCoh3d8_subelem), allocatable, optional, intent(out)  :: top_subelem
  type(fCoh3d8_subelem), allocatable, optional, intent(out)  :: bot_subelem

  if (present(top_subelem_set)) top_subelem_set = elem%top_subelem_set
  if (present(bot_subelem_set)) bot_subelem_set = elem%bot_subelem_set

  if (present(intact_elem)) then
    if (allocated(elem%intact_elem)) then
      allocate(intact_elem)
      intact_elem = elem%intact_elem
    end if
  end if

  if (present(top_subelem)) then
    if (allocated(elem%top_subelem)) then
      allocate(top_subelem)
      top_subelem = elem%top_subelem
    end if
  end if

  if (present(bot_subelem)) then
    if (allocated(elem%bot_subelem)) then
      allocate(bot_subelem)
      bot_subelem = elem%bot_subelem
    end if
  end if

end subroutine extract_fCoh3d8_element



pure subroutine integrate_fCoh3d8_element (elem, nodes, material, K_matrix, &
& F_vector, istat, emsg, nofailure)

use parameter_module, only : DP, MSGLENGTH, STAT_FAILURE, STAT_SUCCESS, ZERO, &
                      & NDIM, HALF
use xnode_module,             only : xnode
use cohesive_material_module, only : cohesive_material
use coh3d8_element_module,    only : integrate
use fCoh3d8_subelem_module,   only : integrate
use global_toolkit_module,    only : assembleKF

  type(fCoh3d8_element),    intent(inout) :: elem
  type(xnode),              intent(inout) :: nodes(NNODE)
  type(cohesive_material),  intent(in)    :: material
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure

  !:::: local variables ::::
  ! local copy of intent inout variables
  type(fCoh3d8_element) :: el
  type(xnode)           :: nds(NNODE), subelem_nds(20)
  ! sub elem K and F
  real(DP), allocatable :: Ki(:,:), Fi(:)
  ! local copy of optional input arg.
  logical :: nofail
  ! error msg location
  character(len=MSGLENGTH) :: msgloc

  ! initialize K & F and local variables
  allocate(K_matrix(NDOF,NDOF), F_vector(NDOF))
  K_matrix = ZERO
  F_vector = ZERO
  istat    = STAT_SUCCESS
  emsg     = ''
  nofail   = .false.
  msgloc   = ' integrate, fCoh3d8_element module'

  ! copy intent inout variables to their local copies
  el  = elem
  nds = nodes
  ! copy optional input to its local copy
  if(present(nofailure)) nofail = nofailure


  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! if intact_elem is present, then the top/bot edges are intact,
  ! just integrate the intact elem would do
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  if (allocated(el%intact_elem)) then

      ! integrate the intact elem
      call integrate (el%intact_elem, nodes=nds(INTACT_ELEM_NODES),          &
      & material=material, K_matrix=Ki, F_vector=Fi, istat=istat, emsg=emsg, &
      & nofailure=nofail)
      if (istat==STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call clean_up(Ki, Fi)
        return
      end if
      ! NOTE : vector subscript for nodes arg is allowed, as it is intent in

      ! assemble the sub elem K and F
      call assembleKF(K_matrix, F_vector, Ki, Fi, INTACT_ELEM_NODES, NDIM, &
      & istat, emsg)
      ! an error is encountered in assembly, zero K and F and exit
      if (istat==STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        K_matrix = ZERO
        F_vector = ZERO
        call clean_up(Ki, Fi)
        return
      end if

      ! update the intent inout dummy args before successful return
      ! note that the nodes are not changed by intact elem, so no need to update
      elem = el

      ! clean up local alloc. array before successful return
      call clean_up(Ki, Fi)

      return

  end if

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! if any subelem is set, integrate and assemble to elem K and F
  ! if both subelems are set, then the elem is integrated twice over its area
  ! half its K and F: K_el = half * (K_top+K_bot), F_el = half * (F_top+F_bot)
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  if (el%top_subelem_set) then
      
      ! MUST copy nds to subelem_nds first, before passing to integrate
      ! fortran does not allow vector subscript for intent out/inout argument
      subelem_nds = nds(TOP_SUBELEM_NODES)
      ! integrate the top subelem
      call integrate (el%top_subelem, nodes=subelem_nds, material=material, &
      & K_matrix=Ki, F_vector=Fi, istat=istat, emsg=emsg, nofailure=nofail)
      if (istat==STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call clean_up(Ki, Fi)
        return
      end if
      ! copy back to nds
      nds(TOP_SUBELEM_NODES) = subelem_nds

      ! assemble the top sub elem K and F
      call assembleKF(K_matrix, F_vector, Ki, Fi, TOP_SUBELEM_NODES, NDIM, &
      & istat, emsg)
      ! an error is encountered in assembly, zero K and F and exit
      if (istat==STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        K_matrix = ZERO
        F_vector = ZERO
        call clean_up(Ki, Fi)
        return
      end if

  end if


  if (el%bot_subelem_set) then

      ! MUST copy nds to subelem_nds first, before passing to integrate
      ! fortran does not allow vector subscript for intent out/inout argument
      subelem_nds = nds(BOT_SUBELEM_NODES)
      ! integrate the bot subelem
      call integrate (el%bot_subelem, nodes=subelem_nds, material=material, &
      & K_matrix=Ki, F_vector=Fi, istat=istat, emsg=emsg, nofailure=nofail)
      if (istat==STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call clean_up(Ki, Fi)
        return
      end if
      ! copy back to nds
      nds(BOT_SUBELEM_NODES) = subelem_nds

      ! assemble the bot sub elem K and F
      call assembleKF(K_matrix, F_vector, Ki, Fi, BOT_SUBELEM_NODES, NDIM, &
      & istat, emsg)
      ! an error is encountered in assembly, zero K and F and exit
      if (istat==STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        K_matrix = ZERO
        F_vector = ZERO
        call clean_up(Ki, Fi)
        return
      end if

  end if


  if (el%top_subelem_set .and. el%bot_subelem_set) then
      ! half the elem's K and F
      K_matrix = HALF * K_matrix
      F_vector = HALF * F_vector
  end if


  ! update the intent inout dummy args before successful return
  ! note that only the internal nodes are changed, and they are stored last
  elem = el
  nodes(NNODE-NNDIN+1 : NNODE) = nds(NNODE-NNDIN+1 : NNODE)

  ! clean up before successful return
  call clean_up(Ki, Fi)

  return


  contains


    pure subroutine clean_up (Ki, Fi)
      real(DP), allocatable, intent(inout) :: Ki(:,:), Fi(:)
      if(allocated(Ki))     deallocate(Ki)
      if(allocated(Fi))     deallocate(Fi)
    end subroutine clean_up


end subroutine integrate_fCoh3d8_element



end module fCoh3d8_element_module
