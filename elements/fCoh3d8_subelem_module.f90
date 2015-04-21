module fCoh3d8_subelem_module
!
!  Purpose:
!    define a floating-node coh3d8 sub-element object, with a breakable top
!    surface.
!
!  topological definition of this element:
!
!  8__14____(E3)__13___7
!  |\                  |\
!  | \15               | \12
!  |  \(E4)            |  \(E2)
!  |   \16             |   \11
!  |    \              |    \
!  |     \5___9__(E1)__|_10__\6
!  |______|____________|      |
! 4\      |  19       3\      |
!   \     |             \     |
!    \    |              \    |
!   20\   |             18\   |
!      \  |                \  |
!       \ |                 \ |
!        \|__________________\|
!         1        17         2
!
!  bottom surface rl. nodes (counter-clock vise): 1, 2, 3, 4
!  top    surface rl. nodes (counter-clock vise): 5, 6, 7, 8
!  top    surface     edges (counter-clock vise): E1, E2, E3, E4
!  bot    surface     edges (counter-clock vise): (same as above, omitted)
!
!  floating nodes on top edges : 9, 10, 11, 12, 13, 14, 15, 16
!  nodes of top E1 ::    <end nodes: 5, 6>   <fl. nodes:  9, 10>
!  nodes of top E2 ::    <end nodes: 6, 7>   <fl. nodes: 11, 12>
!  nodes of top E3 ::    <end nodes: 7, 8>   <fl. nodes: 13, 14>
!  nodes of top E4 ::    <end nodes: 8, 5>   <fl. nodes: 15, 16>
!
!  internal nodes on bottom edges : 17, 18, 19, 20
!  nodes of bot E1 ::    <end nodes: 1, 2>   <in. nodes: 17>
!  nodes of bot E2 ::    <end nodes: 2, 3>   <in. nodes: 18>
!  nodes of bot E3 ::    <end nodes: 3, 4>   <in. nodes: 19>
!  nodes of bot E4 ::    <end nodes: 4, 1>   <in. nodes: 20>
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    19/04/15  B. Y. Chen            Original code
!

! load the modules and entities used for derived type declaration only
use parameter_module, only : NDIM, DP, ZERO, INT_ALLOC_ARRAY
use baseCoh_element_module, only : baseCoh_element

implicit none

private

! list of private parameters in this module:
! NNDRL     : no. of real     nodes in this element
! NNDFL     : no. of floating nodes in this element
! NEDGE_TOP : no. of edges in this element on the top surface
! NNDIN     : no. of internal nodes in this element
! NNODE     : total no. of nodes in this element
! NDOF      : no. of degree-of-freedom in this element
! NODES_ON_TOP_EDGES(i,j) : this is a matrix to describe the nodal connectivity
!             of the top edges in this element. Each edge is topologically
!             composed of its two end nodes and two floating nodes. This
!             term stores the local nodal index of the i_th node on edge j
integer, parameter :: NNDRL=8, NEDGE_TOP=4, NEDGE_BOT=4, NNDFL=2*NEDGE_TOP, &
                    & NNDIN=1*NEDGE_BOT, NNODE=NNDRL+NNDFL+NNDIN, NDOF=NDIM*NNODE
integer, parameter :: NODES_ON_TOP_EDGES(4,NEDGE_TOP) = &
         & reshape([5,6,9,10, 6,7,11,12, 7,8,13,14, 8,5,15,16], [4,NEDGE_TOP])
integer, parameter :: NODES_ON_BOT_EDGES(3,NEDGE_BOT) = &
         & reshape([1,2,17,   2,3,18,    3,4,19,    4,1,20   ], [3,NEDGE_BOT])


type, public :: fCoh3d8_subelem
  private
  ! list of components of this type:
  ! pstat               : elem partition status
  ! node_connec         : nodes global connectivity
  ! edge_connec         : edges global connectivity
  ! lcl_ID_crack_edges  : local indices of cracked edges, this array is passed
  !                       from adj. ply elem. it is used to extract no. and ID
  !                       of cracked edges
  ! edge_lambda         : = (distance btw the crack pnt (if present) on an edge
  !                       and the endnode 1 of the edge) / (length of the edge)
  ! subelem             : list of subelems of type baseCoh
  ! subelem_lcl_connec  : lcl connec of subelems' nodes
  integer  :: pstat                         = 0
  integer  :: node_connec(NNODE)            = 0
  integer  :: edge_connec(NEDGE_TOP)        = 0
  integer  :: lcl_ID_crack_edges(NEDGE_TOP) = 0
  real(DP) :: edge_lambda(NEDGE_BOT)        = ZERO
  type(baseCoh_element), allocatable :: subelem(:)
  type(INT_ALLOC_ARRAY), allocatable :: subelem_lcl_connec(:)
end type fCoh3d8_subelem

interface empty
    module procedure empty_fCoh3d8_subelem
end interface

interface set
    module procedure set_fCoh3d8_subelem
end interface

interface update
    module procedure update_fCoh3d8_subelem
end interface

interface integrate
    module procedure integrate_fCoh3d8_subelem
end interface

interface extract
    module procedure extract_fCoh3d8_subelem
end interface


public :: empty, set, update, integrate, extract




contains




pure subroutine empty_fCoh3d8_subelem (elem)

  type(fCoh3d8_subelem), intent(inout) :: elem

  type(fCoh3d8_subelem) :: elem_lcl

  elem = elem_lcl

end subroutine empty_fCoh3d8_subelem



pure subroutine set_fCoh3d8_subelem (elem, node_connec, edge_connec, istat, emsg)
! Purpose:
! to set the element ready for first use
use parameter_module, only : MSGLENGTH, STAT_FAILURE, STAT_SUCCESS
use baseCoh_element_module, only : set

  type(fCoh3d8_subelem),    intent(inout) :: elem
  integer,                  intent(in)    :: node_connec(NNODE)
  integer,                  intent(in)    :: edge_connec(NEDGE_TOP)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local copy of elem
  type(fCoh3d8_subelem) :: elem_lcl
  ! global_connec of sub element
  integer, allocatable  :: global_connec(:)
  integer :: i

  istat = STAT_SUCCESS
  emsg  = ''
  i = 0

  ! check validity of inputs
  if ( any(node_connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'node connec indices must be >=1, set, &
    &fCoh3d8_subelem_module'
    return
  end if

  if ( any(edge_connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'edge connec indices must be >=1, set, &
    &fCoh3d8_subelem_module'
    return
  end if

  ! update to elem_lcl first
  elem_lcl%node_connec = node_connec
  elem_lcl%edge_connec = edge_connec

  ! allocate 1 coh3d8 subelem
  allocate(elem_lcl%subelem(1))
  allocate(elem_lcl%subelem_lcl_connec(1))
  allocate(elem_lcl%subelem_lcl_connec(1)%array(NNDRL))
  allocate(global_connec(NNDRL))

  ! sub elem 1 local connec
  elem_lcl%subelem_lcl_connec(1)%array = [(i, i=1,NNDRL)]

  ! sub elem 1 global connec
  global_connec(:) = elem_lcl%node_connec( elem_lcl%subelem_lcl_connec(1)%array(:) )

  ! set sub elem 1
  call set (elem_lcl%subelem(1), eltype='coh3d8', connec=global_connec, &
  & istat=istat, emsg=emsg)

  ! if an error is encountered in set, clean up and exit program
  if (istat == STAT_FAILURE) then
    if (allocated(global_connec)) deallocate(global_connec)
    return
  end if

  ! update to dummy arg. elem before successful return
  elem = elem_lcl

  if (allocated(global_connec)) deallocate(global_connec)

end subroutine set_fCoh3d8_subelem



pure subroutine update_fCoh3d8_subelem (elem, lcl_ID_crack_edges, istat, emsg)
! Purpose:
! this subroutine is used solely to update the lcl_ID_crack_edges array of the element
! the lcl_ID_crack_edges is passed in from the adjacent ply element
use parameter_module, only : MSGLENGTH, STAT_FAILURE, STAT_SUCCESS

  type(fCoh3d8_subelem),    intent(inout) :: elem
  integer,                  intent(in)    :: lcl_ID_crack_edges(NEDGE_TOP)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  istat = STAT_SUCCESS
  emsg  = ''

  ! check validity of inputs
  if ( any(lcl_ID_crack_edges < 0) .or. any(lcl_ID_crack_edges > NEDGE_TOP) ) then
    istat = STAT_FAILURE
    emsg  = "cracked edges' local indices must be within [0, 4], update, &
    &fCoh3d8_subelem_module"
    return
  end if

  elem%lcl_ID_crack_edges = lcl_ID_crack_edges

end subroutine update_fCoh3d8_subelem



pure subroutine extract_fCoh3d8_subelem (elem, pstat, node_connec, edge_connec,&
& lcl_ID_crack_edges, edge_lambda, subelem, subelem_lcl_connec)
! Purpose:
! to extract components of this element, generally used in output
! ** Note: **
! to output this element, subelem array needs to be extracted first, then
! extract the connec, fstat, dm, traction, separation, etc. from each sub elem
use parameter_module, only : DP, INT_ALLOC_ARRAY
use baseCoh_element_module, only : baseCoh_element

  type(fCoh3d8_subelem),                        intent(in)  :: elem
  integer,                            optional, intent(out) :: pstat
  integer,               allocatable, optional, intent(out) :: node_connec(:)
  integer,               allocatable, optional, intent(out) :: edge_connec(:)
  integer,               allocatable, optional, intent(out) :: lcl_ID_crack_edges(:)
  real(DP),              allocatable, optional, intent(out) :: edge_lambda(:)
  type(baseCoh_element), allocatable, optional, intent(out) :: subelem(:)
  type(INT_ALLOC_ARRAY), allocatable, optional, intent(out) :: subelem_lcl_connec(:)

  if(present(pstat)) pstat=elem%pstat

  if(present(node_connec)) then
      allocate(node_connec(NNODE))
      node_connec=elem%node_connec
  end if

  if(present(edge_connec)) then
      allocate(edge_connec(NEDGE_TOP))
      edge_connec=elem%edge_connec
  end if

  if(present(lcl_ID_crack_edges)) then
      allocate(lcl_ID_crack_edges(NEDGE_TOP))
      lcl_ID_crack_edges=elem%lcl_ID_crack_edges
  end if

  if(present(edge_lambda)) then
      allocate(edge_lambda(NEDGE_BOT))
      edge_lambda=elem%edge_lambda
  end if

  if(present(subelem)) then
    if(allocated(elem%subelem)) then
      allocate(subelem(size(elem%subelem)))
      subelem=elem%subelem
    end if
  end if

  if(present(subelem_lcl_connec)) then
    if(allocated(elem%subelem_lcl_connec)) then
      allocate(subelem_lcl_connec(size(elem%subelem_lcl_connec)))
      subelem_lcl_connec=elem%subelem_lcl_connec
    end if
  end if

end subroutine extract_fCoh3d8_subelem



pure subroutine integrate_fCoh3d8_subelem (elem, nodes, edge_status, material, &
& K_matrix, F_vector, istat, emsg, nofailure)
! Purpose:
! to integrate this fCoh3d8 subelem and update its internal nodes in nodes array
use parameter_module, only : DP, MSGLENGTH, STAT_FAILURE, STAT_SUCCESS, ZERO, &
                      & INTACT, PARTITIONED_FCOHSUB
use xnode_module, only : xnode
use cohesive_material_module, only : cohesive_material

  type(fCoh3d8_subelem),    intent(inout) :: elem
  type(xnode),              intent(inout) :: nodes(NNODE)
  integer,                  intent(in)    :: edge_status(NEDGE_TOP)
  type(cohesive_material),  intent(in)    :: material
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure

  ! local variables
  type(fCoh3d8_subelem) :: el
  type(xnode) :: nodes_lcl(NNODE)
  logical :: nofail

  ! initialize K & F and local variables
  allocate(K_matrix(NDOF,NDOF), F_vector(NDOF))
  K_matrix = ZERO
  F_vector = ZERO
  istat    = STAT_SUCCESS
  emsg     = ''
  nofail   = .false.

  ! check for input
  if (.not. allocated(elem%subelem)) then
    istat = STAT_FAILURE
    emsg  = 'sub element is NOT yet allocated, integrate, fCoh3d8_subelem module'
    return
  end if

  ! define local variables

  ! copy intent inout arg. to their local copies
  el = elem
  nodes_lcl = nodes

  ! by default, damage modelling is allowed, unless specified
  if(present(nofailure)) nofail = nofailure


  !**** MAIN CALCULATIONS ****
  !
  ! select what to do based on elem's pstat:
  ! if elem is at INTACT:
  !   - check edge status variables to see if the neighbour top ply has matrix
  !     crack; if so, update pstat to PARTITIONED_FCOHSUB and partition sub elems.
  !   - interpolate internal nodes and update to global node list
  !   - integrate and assemble sub elems
  ! if elem is at PARTITIONED_FCOHSUB:
  !   - no need to check edge status variables; no need to update partition.
  !   - interpolate internal nodes and update to global node list
  !   - integrate and assemble sub elems
  ! other status of pstat variable is NOT allowed
  !
  select case (el%pstat)

    case (INTACT)
    ! elem is not yet partitioned

        ! check for edge status, update pstat
        call edge_status_update (el, edge_status, istat, emsg)
        if (istat == STAT_FAILURE) then
          call clean_up(K_matrix, F_vector)
          return
        end if

        ! partition element into sub elems and update internal nodes when
        ! the elem reaches PARTITIONED_FCOHSUB status
        if (el%pstat == PARTITIONED_FCOHSUB) then

          ! partition elem into sub elems
          call partition_element (el, nodes_lcl, istat, emsg)
          if (istat == STAT_FAILURE) then
            call clean_up(K_matrix, F_vector)
            return
          end if

          ! update internal nodes
          call update_internal_nodes (el, nodes_lcl)

        end if

        ! integrate sub elems and assemble system matrices
        call integrate_assemble_subelem(el, nodes_lcl, material, K_matrix, &
        & F_vector, istat, emsg, nofail)
        if (istat == STAT_FAILURE) then
          call clean_up(K_matrix, F_vector)
          return
        end if



    case (PARTITIONED_FCOHSUB)
    ! elem is already partitioned

        ! update internal nodes
        call update_internal_nodes (el, nodes_lcl)

        ! integrate sub elems and assemble system matrices
        call integrate_assemble_subelem(el, nodes_lcl, material, K_matrix, &
        & F_vector, istat, emsg, nofail)
        if (istat == STAT_FAILURE) then
          call clean_up(K_matrix, F_vector)
          return
        end if

    case default
    ! this place should NOT be reached
        istat = STAT_FAILURE
        emsg  = 'unsupported pstat value in fCoh3d8 elem module'
        call clean_up(K_matrix, F_vector)
        return

  end select

  !**** END MAIN CALCULATIONS ****

  ! update to dummy arg. elem before successful return;
  ! for nodes, only need to update the internal nodes stored in the end
  elem  = el
  nodes(NNODE-NNDIN+1 : NNODE) = nodes_lcl(NNODE-NNDIN+1 : NNODE)


  contains


    pure subroutine clean_up (K_matrix, F_vector)
      real(DP), intent(inout) :: K_matrix(:,:), F_vector(:)
      K_matrix = ZERO
      F_vector = ZERO
    end subroutine clean_up

end subroutine integrate_fCoh3d8_subelem



pure subroutine edge_status_update (el, edge_status, istat, emsg)
! Purpose:
! update pass arg. el's pstat and its partition w.r.t its edge status variables
use parameter_module, only : MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                      & PARTITIONED_FCOHSUB, TRANSITION_EDGE, COH_CRACK_EDGE

  ! passed-in variables
  type(fCoh3d8_subelem),    intent(inout) :: el
  integer,                  intent(in)    :: edge_status(NEDGE_TOP)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variables
  ! n_crackedges : no. of cracked edges in this elem
  ! edge1, edge2 : indices of broken edges
  integer :: n_crackedges
  integer :: edge1, edge2
  integer :: i ! counters

  ! ----------------------------------------------------------------------------
  ! *** workings of edge_status, n_crackedges, lcl_ID_crack_edges ***
  !
  ! e.g.: element edge 1 and 3 are broken, then:
  !
  ! - n_crackedges=2
  ! - edge_status(1)>0; edge_status(2)=0; edge_status(3)>0; edge_status(4)=0
  ! - lcl_ID_crack_edges(1)=1; lcl_ID_crack_edges(2)=3; lcl_ID_crack_edges(3:)=0
  !
  ! ----------------------------------------------------------------------------

  ! initialize intent out and local variables
  istat = STAT_SUCCESS
  emsg  = ''
  n_crackedges = 0
  edge1 = 0
  edge2 = 0
  i = 0

  ! Note: no need to check for input validity or create local copy of intent
  ! inout dummy arg because this is an internal procedure

  ! no need to update partition if elem is already in final partition
  if(el%pstat == PARTITIONED_FCOHSUB) return

  ! find the no. of broken edges; local indices of broken edges are stored in
  ! lcl_ID_crack_edges array. this array is passed from adj. ply elem
  !************************* IMPORTANT NOTE: **** ***************************!
  ! use this array instead of edge_status array to extract no. of cracked
  ! edges and local indices of cracked edges.
  !**************************************************************************!
  n_crackedges = count (el%lcl_ID_crack_edges > 0)


  !**** MAIN CALCULATIONS ****

  ! update el%pstat w.r.t no. of failed edges and edge status
  ! only update el%pstat to the final partition status, PARTITIONED_FCOHSUB, when
  ! the top ply elem has reached its final partition status (MATRIX_CRACK_ELEM
  ! or FIBRE_FAILED_ELEM)
  select case (n_crackedges)

    case (0)
    ! adj. ply elem remains INTACT, do nothing
        continue

    case (1)
    ! adj. ply elem is in transition partition; edge status should be
    ! TRANSITION_EDGE
        if(edge_status(el%lcl_ID_crack_edges(1)) == TRANSITION_EDGE) then
          ! edge marks the refinement end and the trans elem start in the upper
          ! ply element. for this fCoh elem, no partition is done for this
          ! case
          continue
        else
          ! lcl_ID_crack_edges is not correct / edge_status not correct
          istat = STAT_FAILURE
          emsg  = 'wrong edge status for n_crackedges=1 in fCoh3d8'
          return
        end if

    case (2)
    ! adj. ply elem could be cracked(final), wake, tip, refinement elem
    ! only update el%pstat if adj. ply elem reaches final partition
        ! store local edge indices to edge 1 and edge2
        edge1 = el%lcl_ID_crack_edges(1)
        edge2 = el%lcl_ID_crack_edges(2)
        if (edge_status(edge1) >= COH_CRACK_EDGE .and. &
        &   edge_status(edge2) >= COH_CRACK_EDGE) then
        ! ply elem status is matrix_crack_elem or fibre_fail_elem, this is the
        ! final partition of the ply elem.
        ! update fCoh pstat
          el%pstat = PARTITIONED_FCOHSUB
        end if

    case default
        istat = STAT_FAILURE
        emsg  = 'unsupported n_crackedges value for edge and el stat update &
        &in fCoh3d8 edge stat partition!'
        return

  end select

  !**** END MAIN CALCULATIONS ****

end subroutine edge_status_update



pure subroutine partition_element (el, nodes, istat, emsg)
! Purpose:
! to partition the element into sub elements
! specifically, this subroutine will set the following el's components:
! - edge_lamda
! - subelem_lcl_connec
! - subelem
! according to the edge status of top edges of this element
use parameter_module, only : MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, INT_ALLOC_ARRAY, &
                      & DP
use xnode_module, only : xnode, extract
use baseCoh_element_module, only : set
use global_toolkit_module, only : distance

  ! passed-in variables
  type(fCoh3d8_subelem),    intent(inout) :: el
  type(xnode),              intent(in)    :: nodes(NNODE)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variables
  ! subelem_glb_connec  : global connec of sub elem nodes
  ! Icrackedge1, Icrackedge2  : indices of cracked edges
  ! e1, e2, e3, e4      : renumbered edge indices to facilitate elem partition
  ! nsub                : no. of sub elements
  ! jcracknode1/2       : fl. node indices on 2 cracked top edges
  ! jendnode1/2         : end node indices of an edge
  ! x1, x2, xc          : coords of crack edge end nodes (x1 & x2) and
  !                       crack tip (xc)

  type(INT_ALLOC_ARRAY), allocatable :: subelem_glb_connec(:)
  integer               :: n_crackedges, Icrackedge1, Icrackedge2
  integer               :: e1, e2, e3, e4
  integer               :: nsub
  integer               :: jcracknode1, jcracknode2, jendnode1, jendnode2
  real(DP), allocatable :: x1(:), x2(:), xc(:)
  ! counters
  integer :: i, j, l

  ! initialize intent out and local variables
  istat = STAT_SUCCESS
  emsg  = ''
  n_crackedges = 0
  Icrackedge1 = 0
  Icrackedge2 = 0
  e1 = 0
  e2 = 0
  e3 = 0
  e4 = 0
  nsub = 0
  jcracknode1 = 0
  jcracknode2 = 0
  jendnode1   = 0
  jendnode2   = 0
  i=0; j=0; l=0

  ! find the no. of broken edges; local indices of broken edges are stored in
  ! lcl_ID_crack_edges array. this array is passed from adj. ply elem
  !************************* IMPORTANT NOTE: **** ***************************!
  ! use this array instead of edge_status array to extract no. of cracked
  ! edges and local indices of cracked edges.
  !**************************************************************************!
  n_crackedges = count (el%lcl_ID_crack_edges > 0)

  select_ncrackedges: select case (n_crackedges)

    case (2) select_ncrackedges
    !- two edges cracked
        ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
        ! determine partition based on the indices of the two broken edges
        ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
        ! find the indices of the two broken edges, sorted in ascending order
        Icrackedge1 = min( el%lcl_ID_crack_edges(1), el%lcl_ID_crack_edges(2) )
        Icrackedge2 = max( el%lcl_ID_crack_edges(1), el%lcl_ID_crack_edges(2) )
        ! Icrackedge1 must be between 1 to 3, and Icrackedge2 between 2 to 4,
        ! with Icrackedge2 > Icrackedge1
        if (Icrackedge1<1 .or. Icrackedge1>3 .or. &
        &   Icrackedge2<2 .or. Icrackedge2>4 .or. &
        &   Icrackedge2 <= Icrackedge1) then
          istat = STAT_FAILURE
          emsg = 'wrong crack edge indices in fCoh3d8 update subelem_lcl_connec &
          & case n_crackedges=2'
          return
        end if
        ! find nsub and e1 - e4
        ! nsub: no. of sub domains
        ! e1 - e4: re-index edges to facilitate partitioning domain
        select case(Icrackedge1)
          case(1)
          ! 1st cracked edge is lcl edge 1, decide on nsub and e1-e4 based on
          ! the lcl ID of the 2nd cracked edge
              select case (Icrackedge2)
                case(2)
                  nsub=4
                  e1=1; e2=2; e3=3; e4=4
                case(3)
                  nsub=2
                  e1=1; e2=2; e3=3; e4=4
                case(4)
                  nsub=4
                  e1=4; e2=1; e3=2; e4=3
                case default
                  istat = STAT_FAILURE
                  emsg = 'wrong 2nd broken edge in update subelem_lcl_connec fCoh3d8'
                  return
              end select
          case(2)
              select case (Icrackedge2)
                case(3)
                  nsub=4
                  e1=2; e2=3; e3=4; e4=1
                case(4)
                  nsub=2
                  e1=2; e2=3; e3=4; e4=1
                case default
                  istat = STAT_FAILURE
                  emsg = 'wrong 2nd broken edge in update subelem_lcl_connec fCoh3d8'
                  return
              end select
          case(3)
              if(Icrackedge2==4) then
                  nsub=4
                  e1=3; e2=4; e3=1; e4=2
              else
                  istat = STAT_FAILURE
                  emsg = 'wrong 2nd broken edge in update subelem_lcl_connec fCoh3d8'
                  return
              end if
          case default
              istat = STAT_FAILURE
              emsg = 'wrong broken edge in update subelem_lcl_connec fCoh3d8'
              return
        end select
        !
        ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
        ! form sub elements based on nsub values
        ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
        select_nsub: select case(nsub)
          ! :::::::::::::::::::::::::::::::::::::::::!
          ! two quad subdomains, two coh3d8 sub elems
          ! :::::::::::::::::::::::::::::::::::::::::!
          case(2) select_nsub
              !:::::::::::::::::::::::::!
              !*** prepare arrays ***
              !:::::::::::::::::::::::::!
              ! allocate nsub no. of subelem, subelem_lcl_connec and
              ! subelem_glb_connec
              if(allocated(el%subelem)) &
              & deallocate(el%subelem)
              if(allocated(el%subelem_lcl_connec)) &
              & deallocate(el%subelem_lcl_connec)
              if(allocated(subelem_glb_connec)) &
              & deallocate(subelem_glb_connec)
              allocate(el%subelem(nsub))
              allocate(el%subelem_lcl_connec(nsub))
              allocate(subelem_glb_connec(nsub))
              ! allocate&initialize the internal arrays of these arrays
              do j=1, nsub
                ! allocate connec for the 8 nodes of coh3d8 sub elems
                allocate(el%subelem_lcl_connec(j)%array(8))
                allocate(subelem_glb_connec(j)%array(8))
                ! initialize these arrays
                el%subelem_lcl_connec(j)%array = 0
                subelem_glb_connec(j)%array    = 0
              end do
              !:::::::::::::::::::::::::!
              !*** calculate lambdas ***
              !:::::::::::::::::::::::::!
              ! find the rel. position of crack pnt on e1 w.r.t its endnode 1
              jcracknode1 = NODES_ON_TOP_EDGES(3,e1)
              jendnode1   = NODES_ON_TOP_EDGES(1,e1)
              jendnode2   = NODES_ON_TOP_EDGES(2,e1)
              call extract(nodes(jendnode1),   x=x1)
              call extract(nodes(jendnode2),   x=x2)
              call extract(nodes(jcracknode1), x=xc)
              el%edge_lambda(e1) = distance(x1,xc,NDIM) / distance(x1,x2,NDIM)
              ! find the rel. position of crack pnt on e3 w.r.t its endnode 1
              jcracknode2 = NODES_ON_TOP_EDGES(3,e3)
              jendnode1   = NODES_ON_TOP_EDGES(1,e3)
              jendnode2   = NODES_ON_TOP_EDGES(2,e3)
              call extract(nodes(jendnode1),   x=x1)
              call extract(nodes(jendnode2),   x=x2)
              call extract(nodes(jcracknode2), x=xc)
              el%edge_lambda(e3) = distance(x1,xc,NDIM) / distance(x1,x2,NDIM)
              !:::::::::::::::::::::::::!
              !*** define sub elm 1 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 4 nodes: e1 nodes 1 & 3, and e3 nodes 3 & 2
              ! top 4 nodes: e1 nodes 1 & 3, and e3 nodes 4 & 2
              el%subelem_lcl_connec(1)%array(1)=NODES_ON_BOT_EDGES(1,e1)
              el%subelem_lcl_connec(1)%array(2)=NODES_ON_BOT_EDGES(3,e1)
              el%subelem_lcl_connec(1)%array(3)=NODES_ON_BOT_EDGES(3,e3)
              el%subelem_lcl_connec(1)%array(4)=NODES_ON_BOT_EDGES(2,e3)
              el%subelem_lcl_connec(1)%array(5)=NODES_ON_TOP_EDGES(1,e1)
              el%subelem_lcl_connec(1)%array(6)=NODES_ON_TOP_EDGES(3,e1)
              el%subelem_lcl_connec(1)%array(7)=NODES_ON_TOP_EDGES(4,e3)
              el%subelem_lcl_connec(1)%array(8)=NODES_ON_TOP_EDGES(2,e3)
              ! define its global connec with global node list
              subelem_glb_connec(1)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(1)%array(:))
              ! set sub element 1
              call set(el%subelem(1), eltype='coh3d8', &
              & connec=subelem_glb_connec(1)%array, istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
              !:::::::::::::::::::::::::!
              !*** define sub elm 2 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 4 nodes: e1 nodes 3 & 2, and e3 nodes 1 & 3
              ! top 4 nodes: e1 nodes 4 & 2, and e3 nodes 1 & 3
              el%subelem_lcl_connec(2)%array(1)=NODES_ON_BOT_EDGES(3,e1)
              el%subelem_lcl_connec(2)%array(2)=NODES_ON_BOT_EDGES(2,e1)
              el%subelem_lcl_connec(2)%array(3)=NODES_ON_BOT_EDGES(1,e3)
              el%subelem_lcl_connec(2)%array(4)=NODES_ON_BOT_EDGES(3,e3)
              el%subelem_lcl_connec(2)%array(5)=NODES_ON_TOP_EDGES(4,e1)
              el%subelem_lcl_connec(2)%array(6)=NODES_ON_TOP_EDGES(2,e1)
              el%subelem_lcl_connec(2)%array(7)=NODES_ON_TOP_EDGES(1,e3)
              el%subelem_lcl_connec(2)%array(8)=NODES_ON_TOP_EDGES(3,e3)
              ! define its global connec with global node list
              subelem_glb_connec(2)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(2)%array(:))
              ! set sub element 2
              call set(el%subelem(2), eltype='coh3d8', &
              & connec=subelem_glb_connec(2)%array, istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
          ! :::::::::::::::::::::::::::::::::::::::::!
          ! four tri subdomains, four coh3d6 sub elems
          ! :::::::::::::::::::::::::::::::::::::::::!
          case(4) select_nsub
              !:::::::::::::::::::::::::!
              !*** prepare arrays ***
              !:::::::::::::::::::::::::!
              ! allocate nsub no. of subelem, subelem_lcl_connec and subelem_T
              if(allocated(el%subelem)) &
              & deallocate(el%subelem)
              if(allocated(el%subelem_lcl_connec)) &
              & deallocate(el%subelem_lcl_connec)
              if(allocated(subelem_glb_connec)) &
              & deallocate(subelem_glb_connec)
              allocate(el%subelem(nsub))
              allocate(el%subelem_lcl_connec(nsub))
              allocate(subelem_glb_connec(nsub))
              do j=1, nsub
                ! allocate connec for 6 nodes of coh3d6 sub elems
                allocate(el%subelem_lcl_connec(j)%array(6))
                allocate(subelem_glb_connec(j)%array(6))
                ! initialize these arrays
                el%subelem_lcl_connec(j)%array = 0
                subelem_glb_connec(j)%array    = 0
              end do
              !:::::::::::::::::::::::::!
              !*** calculate lambdas ***
              !:::::::::::::::::::::::::!
              ! find the rel. position of crack pnt on e1 w.r.t its endnode 1
              jcracknode1 = NODES_ON_TOP_EDGES(3,e1)
              jendnode1   = NODES_ON_TOP_EDGES(1,e1)
              jendnode2   = NODES_ON_TOP_EDGES(2,e1)
              call extract(nodes(jendnode1),   x=x1)
              call extract(nodes(jendnode2),   x=x2)
              call extract(nodes(jcracknode1), x=xc)
              el%edge_lambda(e1) = distance(x1,xc,NDIM) / distance(x1,x2,NDIM)
              ! find the rel. position of crack pnt on e2 w.r.t its endnode 1
              jcracknode2 = NODES_ON_TOP_EDGES(3,e2)
              jendnode1   = NODES_ON_TOP_EDGES(1,e2)
              jendnode2   = NODES_ON_TOP_EDGES(2,e2)
              call extract(nodes(jendnode1),   x=x1)
              call extract(nodes(jendnode2),   x=x2)
              call extract(nodes(jcracknode2), x=xc)
              el%edge_lambda(e2) = distance(x1,xc,NDIM) / distance(x1,x2,NDIM)
              !:::::::::::::::::::::::::!
              !*** define sub elm 1 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 3 nodes: e1 nodes 3, and e2 nodes 1 & 3
              ! top 3 nodes: e1 nodes 4, and e2 nodes 1 & 3
              el%subelem_lcl_connec(1)%array(1)=NODES_ON_BOT_EDGES(3,e1)
              el%subelem_lcl_connec(1)%array(2)=NODES_ON_BOT_EDGES(1,e2)
              el%subelem_lcl_connec(1)%array(3)=NODES_ON_BOT_EDGES(3,e2)
              el%subelem_lcl_connec(1)%array(4)=NODES_ON_TOP_EDGES(4,e1)
              el%subelem_lcl_connec(1)%array(5)=NODES_ON_TOP_EDGES(1,e2)
              el%subelem_lcl_connec(1)%array(6)=NODES_ON_TOP_EDGES(3,e2)
              subelem_glb_connec(1)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(1)%array(:))
              ! set sub elem 1
              call set(el%subelem(1), eltype='coh3d6', &
              & connec=subelem_glb_connec(1)%array, istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
              !:::::::::::::::::::::::::!
              !*** define sub elm 2 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 3 nodes: e2 node 3, and e3 nodes 1 & 2
              ! top 3 nodes: e2 node 4, and e3 nodes 1 & 2
              el%subelem_lcl_connec(2)%array(1)=NODES_ON_BOT_EDGES(3,e2)
              el%subelem_lcl_connec(2)%array(2)=NODES_ON_BOT_EDGES(1,e3)
              el%subelem_lcl_connec(2)%array(3)=NODES_ON_BOT_EDGES(2,e3)
              el%subelem_lcl_connec(2)%array(4)=NODES_ON_TOP_EDGES(4,e2)
              el%subelem_lcl_connec(2)%array(5)=NODES_ON_TOP_EDGES(1,e3)
              el%subelem_lcl_connec(2)%array(6)=NODES_ON_TOP_EDGES(2,e3)
              subelem_glb_connec(2)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(2)%array(:))
              ! set sub elem 2
              call set(el%subelem(2),eltype='coh3d6', &
              & connec=subelem_glb_connec(2)%array, istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
              !:::::::::::::::::::::::::!
              !*** define sub elm 3 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 3 nodes: e4 nodes 1 & 2, and e1 node 3
              ! top 3 nodes: e4 nodes 1 & 2, and e1 node 3
              el%subelem_lcl_connec(3)%array(1)=NODES_ON_BOT_EDGES(1,e4)
              el%subelem_lcl_connec(3)%array(2)=NODES_ON_BOT_EDGES(2,e4)
              el%subelem_lcl_connec(3)%array(3)=NODES_ON_BOT_EDGES(3,e1)
              el%subelem_lcl_connec(3)%array(4)=NODES_ON_TOP_EDGES(1,e4)
              el%subelem_lcl_connec(3)%array(5)=NODES_ON_TOP_EDGES(2,e4)
              el%subelem_lcl_connec(3)%array(6)=NODES_ON_TOP_EDGES(3,e1)
              subelem_glb_connec(3)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(3)%array(:))
              ! set sub elem 3
              call set(el%subelem(3),eltype='coh3d6', &
              & connec=subelem_glb_connec(3)%array, istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
              !:::::::::::::::::::::::::!
              !*** define sub elm 4 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 3 nodes: e1 node 3, e2 node 3 and e3 node 2
              ! top 3 nodes: e1 node 3, e2 node 4 and e3 node 2
              el%subelem_lcl_connec(4)%array(1)=NODES_ON_BOT_EDGES(3,e1)
              el%subelem_lcl_connec(4)%array(2)=NODES_ON_BOT_EDGES(3,e2)
              el%subelem_lcl_connec(4)%array(3)=NODES_ON_BOT_EDGES(2,e3)
              el%subelem_lcl_connec(4)%array(4)=NODES_ON_TOP_EDGES(3,e1)
              el%subelem_lcl_connec(4)%array(5)=NODES_ON_TOP_EDGES(4,e2)
              el%subelem_lcl_connec(4)%array(6)=NODES_ON_TOP_EDGES(2,e3)
              subelem_glb_connec(4)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(4)%array(:))
              ! set sub elem 4
              call set(el%subelem(4),eltype='coh3d6', &
              & connec=subelem_glb_connec(4)%array, istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
          ! :::::::::::::::::::::::::::::::::::::::::!
          ! unsupported no. of sub elems, ERROR
          ! :::::::::::::::::::::::::::::::::::::::::!
          case default select_nsub
              istat = STAT_FAILURE
              emsg = 'wrong nsub in update subelem_lcl_connec fCoh3d8'
              call clean_up (subelem_glb_connec, x1, x2, xc)
              return
        end select select_nsub

    case default select_ncrackedges
        istat = STAT_FAILURE
        emsg = 'unexpected n_crackedges value in fCoh3d8 update subelem_lcl_connec!'
        call clean_up (subelem_glb_connec, x1, x2, xc)
        return

    end select select_ncrackedges


    ! deallocate local array

    call clean_up (subelem_glb_connec, x1, x2, xc)



    contains



      pure subroutine clean_up (subelem_glb_connec, x1, x2, xc)
        type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subelem_glb_connec(:)
        real(DP),              allocatable, intent(inout) :: x1(:), x2(:), xc(:)

        if(allocated(subelem_glb_connec)) deallocate(subelem_glb_connec)
        if(allocated(x1)) deallocate(x1)
        if(allocated(x2)) deallocate(x2)
        if(allocated(xc)) deallocate(xc)
      end subroutine clean_up


end subroutine partition_element



pure subroutine update_internal_nodes (el, nodes)
! Purpose:
! - update the components of the internal nodes of this element
use parameter_module, only : ZERO, ONE
use xnode_module, only : xnode, operator(+), operator(*)

  type(fCoh3d8_subelem), intent(in)    :: el
  type(xnode),           intent(inout) :: nodes(NNODE)

  integer     :: Icrackedge, Iinnode, Iendnode1, Iendnode2
  real(DP)    :: lambda
  type(xnode) :: innode, endnode1, endnode2
  integer     :: j

  ! initialize local var.
  Icrackedge = 0
  Iinnode    = 0
  Iendnode1  = 0
  Iendnode2  = 0
  lambda     = ZERO
  j = 0

  do j = 1, NEDGE_TOP
    ! extract the jth cracked edge local index
    Icrackedge  = el%lcl_ID_crack_edges(j)
    ! if it is 0, then no more cracked edges, exit loop
    if (Icrackedge == 0) exit
    ! if it is not 0, then find the internal node on the corresponding
    ! bottom edge, and its two end nodes
    Iinnode     = NODES_ON_BOT_EDGES(3, Icrackedge)
    Iendnode1   = NODES_ON_BOT_EDGES(1, Icrackedge)
    Iendnode2   = NODES_ON_BOT_EDGES(2, Icrackedge)
    ! find the stored lambda of this edge
    lambda      = el%edge_lambda(Icrackedge)
    ! interpolate this internal node with end nodes
    endnode1 = nodes(Iendnode1)
    endnode2 = nodes(Iendnode2)
    innode   = (ONE-lambda) * endnode1 + lambda * endnode2
    ! update this internal node in node list
    nodes(Iinnode) = innode
  end do

end subroutine update_internal_nodes



pure subroutine integrate_assemble_subelem (elem, nodes, material, K_matrix, &
& F_vector, istat, emsg, nofailure)
! Purpose:
! integrate and assemble sub element system arrays
use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, ZERO, &
                      & INTACT, PARTITIONED_FCOHSUB, NDIM
use xnode_module, only : xnode
use cohesive_material_module, only : cohesive_material
use baseCoh_element_module, only : integrate
use global_toolkit_module, only : assembleKF

  ! - passed in variables
  type(fCoh3d8_subelem),    intent(inout) :: elem
  type(xnode),              intent(in)    :: nodes(NNODE)
  type(cohesive_material),  intent(in)    :: material
  real(DP),                 intent(out)   :: K_matrix(NDOF,NDOF), F_vector(NDOF)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure

  ! - local variables
  real(DP), allocatable :: Ki(:,:), Fi(:)
  integer,  allocatable :: dofcnc(:)
  real(DP), allocatable :: Kmat_r(:,:), Fvec_r(:)
  real(DP), allocatable :: Tmatrixfull(:,:)
  logical :: nofail
  integer :: isub, j, k, l, n

  ! initialize intent out and local variables
  K_matrix = ZERO
  F_vector = ZERO
  istat  = STAT_SUCCESS
  emsg   = ''
  nofail = .false.
  isub=0; j=0; k=0; l=0; n=0


  if (present(nofailure)) nofail = nofailure

  ! loop over all sub elements to
  ! integrate their system arrays and assemble together
  if (elem%pstat == INTACT) then
  ! element is not yet partitioned, there should be only 1 sub elem allocated
  ! only integrate subelem 1 without using mnodes

      do isub = 1, size(elem%subelem)

        ! verify that there's only 1 subelem
        if (isub /= 1) then
          istat = STAT_FAILURE
          emsg = 'incompatible no. of subelem with INTACT partition status, &
          & integrate_assemble_subelem, fCoh3d8_subelem_module'
          exit
        end if

        ! call the integrate procedure asoc. with the subelem type
        ! exit the do loop (and if-else subsequently) if an error is encountered
        call integrate(elem%subelem(isub), &
        & nodes(elem%subelem_lcl_connec(isub)%array(:)), material, Ki, Fi, &
        & istat, emsg, nofail)
        if (istat == STAT_FAILURE) exit

        ! if no error, prepare to assemble Ki and Fi to elem's K and F
        if(allocated(dofcnc)) deallocate(dofcnc)
        allocate(dofcnc(size(Fi))); dofcnc = 0

        ! loop over no. of nodes in sub elem i, = NNDRL
        do j = 1, NNDRL
          do l = 1, NDIM
            ! obtain dof indices of the jth node of sub elem i
            dofcnc( (j-1) * NDIM + l ) = &
            & ( elem%subelem_lcl_connec(isub)%array(j) - 1 ) * NDIM + l
          end do
        end do

        ! assemble the sub elem K and F
        call assembleKF(K_matrix, F_vector, Ki, Fi, dofcnc, istat, emsg)
        if (istat == STAT_FAILURE) exit

      end do

  else if (elem%pstat == PARTITIONED_FCOHSUB) then
  ! element is partitioned into subelems

      !:::::::::::::::::::::::::::::::::::::!
      ! integrate and assemble sub elems
      !:::::::::::::::::::::::::::::::::::::!
      do isub = 1, size(elem%subelem)

        ! call the integrate procedure asoc. with the subelem type
        ! exit the do loop (and if-else subsequently) if an error is encountered
        call integrate(elem%subelem(isub), &
        & nodes(elem%subelem_lcl_connec(isub)%array(:)), material, Ki, Fi, &
        & istat, emsg, nofail)
        if (istat == STAT_FAILURE) exit

        ! prepare to assemble sub K and sub F to elem's K and F
        if(allocated(dofcnc)) deallocate(dofcnc)
        allocate(dofcnc(size(Fi))); dofcnc = 0

        ! loop over no. of nodes in sub elem i
        do j = 1, size(elem%subelem_lcl_connec(isub)%array)
          do l = 1, NDIM
            ! dof indices of the jth node of sub elem i
            dofcnc( (j-1) * NDIM + l ) = &
            & ( elem%subelem_lcl_connec(isub)%array(j) - 1 ) * NDIM + l
          end do
        end do

        ! assemble the sub elem K and F
        call assembleKF(K_matrix, F_vector, Ki, Fi, dofcnc, istat, emsg)
        if (istat == STAT_FAILURE) exit

      end do

      !:::::::::::::::::::::::::::::::::::::!
      ! condense the terms of internal nodes
      !:::::::::::::::::::::::::::::::::::::!
      if (istat == STAT_SUCCESS) then

        ! get the Tmatrix for condensing the internal nodes
        call get_Tmatrix (elem, Tmatrixfull)

        ! allocate the condensed K mat and F vec
        n = size(Tmatrixfull(1,:))
        allocate( Kmat_r(n,n), Fvec_r(n) )
        Kmat_r = ZERO
        Fvec_r = ZERO

        ! calculate the condensed K matrix and F vector
        Kmat_r = matmul( matmul( transpose(Tmatrixfull), K_matrix ), Tmatrixfull )
        Fvec_r = matmul( transpose(Tmatrixfull), F_vector )

        ! zero K_matrix and F_vector for re-value
        K_matrix = ZERO
        F_vector = ZERO

        ! copy the condensed K and F into K_matrix and F_vector
        K_matrix(1:n, 1:n) = Kmat_r(:,:)
        F_vector(1:n)      = Fvec_r(:)

      end if

  else
      istat = STAT_FAILURE
      emsg = 'unsupported elem pstat value in integrate_assemble_subelem, &
      & fCoh3d8_subelem_module'

  end if

  ! clean up if above loop exit upon error
  if (istat == STAT_FAILURE) then
    K_matrix = ZERO
    F_vector = ZERO
  end if

  if(allocated(Ki))             deallocate(Ki)
  if(allocated(Fi))             deallocate(Fi)
  if(allocated(dofcnc))         deallocate(dofcnc)
  if(allocated(Tmatrixfull))    deallocate(Tmatrixfull)
  if(allocated(Kmat_r))         deallocate(Kmat_r)
  if(allocated(Fvec_r))         deallocate(Fvec_r)

end subroutine integrate_assemble_subelem



pure subroutine get_Tmatrix (el, Tmatrixfull)
! Purpose:
! - calculate [Tmatrixfull], s.t. {U} = [Tmatrixfull] . {U_condensed}, where
! {U_condensed} is the nodal dof vector without the internal nodal dofs. The
! internal nodal dofs are interpolated from the other nodal dofs.
use parameter_module, only : DP, ZERO, ONE

  type(fCoh3d8_subelem), intent(in)  :: el
  real(DP), allocatable, intent(out) :: Tmatrixfull(:,:)

  ! local var.
  real(DP), allocatable :: Tmatrix(:,:)
  integer  :: Icrackedge, Iinnode, Iendnode1, Iendnode2
  real(DP) :: lambda
  integer  :: j, k, l

  ! initialize intent out and local variables
  Icrackedge = 0
  Iinnode    = 0
  Iendnode1  = 0
  Iendnode2  = 0
  lambda     = ZERO
  j = 0; k = 0; l = 0

  ! allocate Tmatrix relating ALL nodes with only the global nodes
  ! - row dimension = NNODE          : all nodes
  ! - col dimension = NNODE - NNDIN  : all nodes exc. internal nodes
  allocate( Tmatrix(NNODE, NNODE-NNDIN) )
  Tmatrix = ZERO

  ! global nodes' terms remain as is
  do j = 1, NNODE-NNDIN
    Tmatrix(j,j) = ONE
  end do

  ! Populate Tmatrix rows corresponding to internal nodes
  ! find all the cracked edges and the corresponding internal nodes for
  ! condensation
  do j = 1, NEDGE_TOP
    ! extract the jth cracked edge local index
    Icrackedge  = el%lcl_ID_crack_edges(j)
    ! if it is 0, then no more cracked edges, exit loop
    if (Icrackedge == 0) exit
    ! if it is not 0, then find the internal node on the corresponding
    ! bottom edge, and its two end nodes
    Iinnode     = NODES_ON_BOT_EDGES(3, Icrackedge)
    Iendnode1   = NODES_ON_BOT_EDGES(1, Icrackedge)
    Iendnode2   = NODES_ON_BOT_EDGES(2, Icrackedge)
    ! find the stored lambda of this edge
    lambda      = el%edge_lambda(Icrackedge)
    ! put lambda in the Tmatrix to interpolate internal node with end nodes:
    ! U_innode = (1-lambda) * U_endnode1 + lambda * U_endnode2
    Tmatrix(Iinnode,Iendnode1) = ONE - lambda
    Tmatrix(Iinnode,Iendnode2) = lambda
  end do

  ! Populate FULL Tmatrix with Tmatrix
  ! full Tmatrix is basically Tmatrix expanded with NDIM
  ! - row dimension = NNODE * NDIM            : all nodal dofs
  ! - col dimension = (NNODE - NNDIN) * NDIM  : all nodal dofs exc. intl. ones
  allocate(Tmatrixfull(NNODE*NDIM, (NNODE-NNDIN)*NDIM))
  Tmatrixfull = ZERO

  ! copy the terms from Tmatrix to Tmatrixfull;
  ! the term Tmatrix(k,j) is expanded to a diagonal submatrix with rank NDIM
  do j = 1, NNODE-NNDIN
    do k = 1, NNODE
      do l = 1, NDIM
        Tmatrixfull( (k-1)*NDIM+l , (j-1)*NDIM+l ) = Tmatrix(k,j)
      end do
    end do
  end do

end subroutine get_Tmatrix





end module fCoh3d8_subelem_module
