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
use parameter_module,       only : NDIM, DP, ZERO, INT_ALLOC_ARRAY
use baseCoh_element_module, only : baseCoh_element

implicit none

private

! list of private parameters in this module:
! NNDRL     : no. of real     nodes in this element
! NEDGE_SURF: no. of edges in this element on the top/bot surface
! NNDFL     : no. of floating nodes in this element
! NNDIN     : no. of internal nodes in this element
! NNODE     : total no. of nodes in this element
! NDOF      : no. of degree-of-freedom in this element
! NODES_ON_TOP_EDGES(i,j) : this is a matrix to describe the nodal connectivity
!             of the top edges in this element. Each edge is topologically
!             composed of its two end nodes and two floating nodes. This
!             term stores the local nodal index of the i_th node on edge j
! NODES_ON_BOT_EDGES(i,j) : this is a matrix to describe the nodal connectivity
!             of the bottom edges in this element. Each edge is topologically
!             composed of its two end nodes and one internal node (duplicated
!             for the sake of convenience when partitioning the bottom surf).
!             This term stores the local nodal index of the i_th node on edge j
integer, parameter :: NNDRL       = 8,                      &
                   &  NEDGE_SURF  = 4,                      &
                   &  NNDFL       = 2 * NEDGE_SURF,         &
                   &  NNDIN       = 1 * NEDGE_SURF,         &
                   &  NNODE       = NNDRL + NNDFL + NNDIN,  &
                   &  NDOF        = NDIM * NNODE
integer, parameter :: NODES_ON_TOP_EDGES(4,NEDGE_SURF) =    &
& reshape([5,6,9,10,   6,7,11,12,  7,8,13,14,  8,5,15,16], [4,NEDGE_SURF])
integer, parameter :: NODES_ON_BOT_EDGES(4,NEDGE_SURF) =    &
& reshape([1,2,17,17,  2,3,18,18,  3,4,19,19,  4,1,20,20], [4,NEDGE_SURF])


type, public :: fCoh3d8_subelem
  private
  ! list of components of this type:
  ! node_connec         : nodes global connectivity
  ! top_edge_status     : edge status of top surf. edges
  ! edge_lambda         : = (distance btw the crack pnt (if present) on an edge
  !                       and the endnode 1 of the edge) / (length of the edge)
  ! subelem             : list of subelems of type baseCoh
  ! subelem_lcl_connec  : lcl connec of subelems' nodes
  integer  :: node_connec(NNODE)            = 0
  integer  :: top_edge_status(NEDGE_SURF)   = 0
  real(DP) :: edge_lambda(NEDGE_SURF)       = ZERO
  type(baseCoh_element), allocatable :: subelem(:)
  type(INT_ALLOC_ARRAY), allocatable :: subelem_lcl_connec(:)
end type fCoh3d8_subelem

interface empty
    module procedure empty_fCoh3d8_subelem
end interface

interface set
    module procedure set_fCoh3d8_subelem
end interface

interface integrate
    module procedure integrate_fCoh3d8_subelem
end interface

interface extract
    module procedure extract_fCoh3d8_subelem
end interface


public :: empty, set, integrate, extract




contains




pure subroutine empty_fCoh3d8_subelem (elem)

  type(fCoh3d8_subelem), intent(inout) :: elem

  type(fCoh3d8_subelem) :: elem_lcl

  elem = elem_lcl

end subroutine empty_fCoh3d8_subelem



pure subroutine extract_fCoh3d8_subelem (elem, subelem)
! Purpose:
! to extract components of this element, generally used in output
! ** Note: **
! to output this element, subelem array needs to be extracted first, then
! extract the connec, fstat, dm, traction, separation, etc. from each sub elem
use baseCoh_element_module, only : baseCoh_element

  type(fCoh3d8_subelem),                        intent(in)  :: elem
  type(baseCoh_element), allocatable, optional, intent(out) :: subelem(:)

  if(present(subelem)) then
    if(allocated(elem%subelem)) then
      allocate(subelem(size(elem%subelem)))
      subelem=elem%subelem
    end if
  end if

end subroutine extract_fCoh3d8_subelem



pure subroutine set_fCoh3d8_subelem (elem, node_connec, top_edge_status, &
& istat, emsg)
! Purpose:
! to set the element ready for first use
use parameter_module, only : MSGLENGTH, STAT_FAILURE, STAT_SUCCESS,&
                      & INTACT, TRANSITION_EDGE, REFINEMENT_EDGE,  &
                      & CRACK_TIP_EDGE, WEAK_CRACK_EDGE,           &
                      & COH_CRACK_EDGE, STRONG_CRACK_EDGE

  type(fCoh3d8_subelem),    intent(inout) :: elem
  integer,                  intent(in)    :: node_connec(NNODE)
  integer,                  intent(in)    :: top_edge_status(NEDGE_SURF)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! location for emsg
  character(len=MSGLENGTH) :: msgloc
  integer :: n_crackedges

  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ' set, fCoh3d8_subelem module'
  n_crackedges = 0
  
  ! check validity of inputs
  
  ! check node connec
  if ( any(node_connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'node connec indices must be >=1,'//trim(msgloc)
    return
  end if
  
  ! check edge status, see if there's any unexpected edge status value
  if ( any( .not. ( top_edge_status == INTACT          .or.         &
  &                 top_edge_status == TRANSITION_EDGE .or.         &
  &                 top_edge_status == REFINEMENT_EDGE .or.         &
  &                 top_edge_status == CRACK_TIP_EDGE  .or.         &
  &                 top_edge_status == WEAK_CRACK_EDGE .or.         &
  &                 top_edge_status == COH_CRACK_EDGE  .or.         &
  &                 top_edge_status == STRONG_CRACK_EDGE )  )  ) then
    istat = STAT_FAILURE
    emsg  = 'edge status value is NOT recognized,'//trim(msgloc)
    return
  end if
  
  ! check the no. of broken edges; only accepts TWO cracked edges
  n_crackedges = count (top_edge_status >= COH_CRACK_EDGE)
  if (n_crackedges /= 2) then
    istat = STAT_FAILURE
    emsg  = 'no. of cracked edges must be TWO,'//trim(msgloc)
    return
  end if

  ! update to elem
  elem%node_connec     = node_connec
  elem%top_edge_status = top_edge_status

end subroutine set_fCoh3d8_subelem




pure subroutine integrate_fCoh3d8_subelem (elem, nodes, material, K_matrix, &
& F_vector, istat, emsg, nofailure)
! Purpose:
! to integrate this fCoh3d8 subelem and update its internal nodes in nodes array
use parameter_module, only : DP, MSGLENGTH, STAT_FAILURE, STAT_SUCCESS, ZERO
use xnode_module,             only : xnode
use cohesive_material_module, only : cohesive_material

  type(fCoh3d8_subelem),    intent(inout) :: elem
  type(xnode),              intent(inout) :: nodes(NNODE)
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

  ! check nodes and material:
  ! nodes x and u must be allocated, and material must be properly defined
  ! but this should NOT be checked here, they should be ensured correct at
  ! the beginning of analysis.

  ! define local variables

  ! copy intent inout arg. to their local copies
  el = elem
  nodes_lcl = nodes

  ! by default, damage modelling is allowed, unless specified
  if(present(nofailure)) nofail = nofailure


  !**** MAIN CALCULATIONS ****

  ! partition element if first time integration
  if (.not. allocated(elem%subelem)) then
    call partition_element (el, nodes_lcl, istat, emsg)
    if (istat == STAT_FAILURE) then
      call clean_up(K_matrix, F_vector)
      return
    end if
  end if

  ! update internal nodes
  call update_internal_nodes (el, nodes_lcl)

  ! integrate sub elems and assemble system matrices
  call integrate_assemble_subelem(el, nodes_lcl, material, K_matrix, &
  & F_vector, istat, emsg, nofail)
  if (istat == STAT_FAILURE) then
    call clean_up(K_matrix, F_vector)
    return
  end if

  !**** END MAIN CALCULATIONS ****

  ! update to dummy arg. elem before successful return;
  ! for nodes, only need to update the internal nodes stored in the end
  elem  = el
  nodes(NNODE-NNDIN+1 : NNODE) = nodes_lcl(NNODE-NNDIN+1 : NNODE)
  return

  contains


    pure subroutine clean_up (K_matrix, F_vector)
      real(DP), intent(inout) :: K_matrix(:,:), F_vector(:)
      K_matrix = ZERO
      F_vector = ZERO
    end subroutine clean_up

end subroutine integrate_fCoh3d8_subelem



pure subroutine partition_element (el, nodes, istat, emsg)
! Purpose:
! to partition the element into sub elements
! specifically, this subroutine will set the following el's components:
! - edge_lamda
! - subelem_lcl_connec
! - subelem
! according to the edge status of top edges of this element
!** NOTE: ONLY PARTITION WITH TWO CRACKED EDGES IS ALLOWED **
use parameter_module, only : MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                      & INT_ALLOC_ARRAY, DP, ELTYPELENGTH, ZERO,    &
                      & ONE, SMALLNUM, COH_CRACK_EDGE
use xnode_module,           only : xnode, extract
use baseCoh_element_module, only : set
use global_toolkit_module,  only : distance, partition_quad_elem

  ! passed-in variables
  type(fCoh3d8_subelem),    intent(inout) :: el
  type(xnode),              intent(in)    :: nodes(NNODE)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variables
  ! msgloc : location of this procedure, attached at the end of error message
  ! lcl_ID_crack_edges        : local indices of cracked edges
  ! n_crackedges              : no. of cracked edges on top surf.
  ! jcrackedge                : lcl ID of a cracked edge
  ! jcracknode                : fl. node indices on cracked top edges
  ! jendnode1/2               : end node indices of an edge
  ! x1, x2, xc                : coords of crack edge end nodes (x1 & x2) and
  !                             crack tip (xc)
  ! distn1n2, distn1nc        : distance btw 2 endnodes of an edge, and distance
  !                             btw endnode1 and crack point on this edge
  ! lambda                    : = distn1nc / distn1n2
  ! subelem_lcl_connec_top    : local connec of sub elem top    surf. nodes
  ! subelem_lcl_connec_bot    : local connec of sub elem bottom surf. nodes
  ! subelem_glb_connec        : global connec of sub elem nodes
  ! nsub                      : no. of sub elements
  ! subelem_nnode             : no. of nodes in a sub elem
  ! subelem_type              : eltype of sub elem
  !
  character(len=MSGLENGTH) :: msgloc
  integer               :: lcl_ID_crack_edges(NEDGE_SURF)
  ! variables to define edge lambda
  integer               :: n_crackedges, jcrackedge
  integer               :: jcracknode, jendnode1, jendnode2
  real(DP), allocatable :: x1(:), x2(:), xc(:)
  real(DP)              :: distn1n2, distn1nc, lambda
  ! variables to define sub elems
  type(INT_ALLOC_ARRAY), allocatable :: subelem_lcl_connec_top(:)
  type(INT_ALLOC_ARRAY), allocatable :: subelem_lcl_connec_bot(:)
  type(INT_ALLOC_ARRAY), allocatable :: subelem_glb_connec(:)
  integer                            :: nsub, subelem_nnode
  character(len=ELTYPELENGTH)        :: subelem_type
  ! counters
  integer :: i, j, l

  ! initialize intent out and local variables
  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ' partition_element, fCoh3d8_subelem_module'
  lcl_ID_crack_edges = 0
  n_crackedges = 0
  jcrackedge   = 0
  jcracknode   = 0
  jendnode1    = 0
  jendnode2    = 0
  distn1n2     = ZERO
  distn1nc     = ZERO
  lambda       = ZERO
  nsub          = 0
  subelem_nnode = 0
  subelem_type  = ''
  i=0; j=0; l=0

  ! Note: no need to check for input validity or create local copy of intent
  ! inout dummy arg because this is an internal procedure

  ! find the no. of broken edges, and store their local indices in the array
  ! lcl_ID_crack_edges
  do j = 1, NEDGE_SURF
    if (el%top_edge_status(j) >= COH_CRACK_EDGE) then
      n_crackedges = n_crackedges + 1
      lcl_ID_crack_edges (n_crackedges) = j
    end if
  end do

  ! check if n_crackedges is 2
  if (n_crackedges /= 2) then
    istat = STAT_FAILURE
    emsg = 'unexpected no. of cracked edges in'//trim(msgloc)
    return
  end if

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! calculate lambdas (crack point relative loc.) of the cracked edges
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  do j = 1, n_crackedges
    ! find the edge index of the jth cracked edge
    jcrackedge = lcl_ID_crack_edges(j)
    ! find the fl. node (crack point) and the real nodes (endnodes) on this edge
    jcracknode = NODES_ON_TOP_EDGES(3,jcrackedge)
    jendnode1  = NODES_ON_TOP_EDGES(1,jcrackedge)
    jendnode2  = NODES_ON_TOP_EDGES(2,jcrackedge)
    ! extract the coordinates of these nodes
    call extract(nodes(jendnode1),  x=x1)
    call extract(nodes(jendnode2),  x=x2)
    call extract(nodes(jcracknode), x=xc)
    ! calculate distance btw endnode1 and endnode 2
    distn1n2 = distance(x1,x2,NDIM)
    ! calculate distance btw endnode1 and crack point
    distn1nc = distance(x1,xc,NDIM)
    ! distance cannot be ZERO
    if (distn1n2 < SMALLNUM .or. distn1nc < SMALLNUM) then
      istat = STAT_FAILURE
      emsg = 'coords of edge end nodes or fl. nodes are incorrect'//trim(msgloc)
      call clean_up (subelem_glb_connec, subelem_lcl_connec_top, &
      & subelem_lcl_connec_bot, x1, x2, xc)
      return
    end if
    ! calculate lambda
     lambda = distn1nc / distn1n2
    ! lambda must be within (0, 1)
    if ( .not. (SMALLNUM < lambda .and. lambda < ONE-SMALLNUM) ) then
      istat = STAT_FAILURE
      emsg  = 'edge lambda is out of range'//trim(msgloc)
      call clean_up (subelem_glb_connec, subelem_lcl_connec_top, &
      & subelem_lcl_connec_bot, x1, x2, xc)
      return
    end if
    ! update lambda to el components
    el%edge_lambda(jcrackedge) = lambda
  end do

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! call the partition quad elem subroutine from global toolkit module
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  ! pass the TOP surface topology NODES_ON_TOP_EDGES into the subroutine,
  ! this will partition the element top surface into 2D sub elems,
  ! whose nodes (in lcl indices) will be stored in subelem_lcl_connec_top
  call partition_quad_elem (NODES_ON_TOP_EDGES, lcl_ID_crack_edges, &
  & subelem_lcl_connec_top, istat, emsg)
  if (istat == STAT_FAILURE) then
    emsg = emsg//trim(msgloc)
    call clean_up (subelem_glb_connec, subelem_lcl_connec_top, &
    & subelem_lcl_connec_bot, x1, x2, xc)
    return
  end if

  ! pass the BOTTOM surface topology NODES_ON_BOT_EDGES into the subroutine,
  ! this will partition the element bottom surface into 2D sub elems,
  ! whose nodes (in lcl indices) will be stored in subelem_lcl_connec_bot
  call partition_quad_elem (NODES_ON_BOT_EDGES, lcl_ID_crack_edges, &
  & subelem_lcl_connec_bot, istat, emsg)
  if (istat == STAT_FAILURE) then
    emsg = emsg//trim(msgloc)
    call clean_up (subelem_glb_connec, subelem_lcl_connec_top, &
    & subelem_lcl_connec_bot, x1, x2, xc)
    return
  end if

  ! ** NOTE **
  ! the lcl_ID_crack_edges array is applicable to both top and bottom edges; so
  ! the top and bottom surfaces are partitioned in exactly the same way, i.e.,
  ! same number of sub elems, same sub elem types.

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! based on the top/bot surface partitions in subelem_lcl_connec_top/bot,
  ! allocate and populate the el%subelem_lcl_connec array
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  ! extract no. of sub elems
  nsub = size(subelem_lcl_connec_top)

  ! allocate nsub no. of subelem, subelem_lcl_connec and
  ! subelem_glb_connec
  if(allocated(el%subelem))            deallocate(el%subelem)
  if(allocated(el%subelem_lcl_connec)) deallocate(el%subelem_lcl_connec)
  if(allocated(subelem_glb_connec))    deallocate(subelem_glb_connec)
  allocate(el%subelem(nsub))
  allocate(el%subelem_lcl_connec(nsub))
  allocate(subelem_glb_connec(nsub))

  ! allocate & define these arrays
  do j = 1, nsub

    ! determine no. of nodes in sub elem j
    subelem_nnode = 2 * size(subelem_lcl_connec_top(j)%array)

    ! determine sub elem type based on subelem_nnode
    select case (subelem_nnode)
      case (6)
        subelem_type = 'coh3d6'
      case (8)
        subelem_type = 'coh3d8'
      case default
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of nodes for sub elem in'//trim(msgloc)
        call clean_up (subelem_glb_connec, subelem_lcl_connec_top, &
        & subelem_lcl_connec_bot, x1, x2, xc)
        return
    end select

    ! allocate connec for sub elems
    allocate(el%subelem_lcl_connec(j)%array(subelem_nnode))
    allocate(subelem_glb_connec(j)%array(subelem_nnode))

    ! initialize these arrays
    el%subelem_lcl_connec(j)%array = 0
    subelem_glb_connec(j)%array    = 0

    ! populate lcl connec
    ! copy top surf. lcl connec
    el%subelem_lcl_connec(j)%array( subelem_nnode/2 + 1 : subelem_nnode   ) = &
    &  subelem_lcl_connec_top(j)%array(:)
    ! copy bot surf. lcl connec
    el%subelem_lcl_connec(j)%array(                   1 : subelem_nnode/2 ) = &
    &  subelem_lcl_connec_bot(j)%array(:)

    ! populate glb connec through el%node_connec
    subelem_glb_connec(j)%array = el%node_connec(el%subelem_lcl_connec(j)%array)

    ! set this sub element
    call set(el%subelem(j), eltype=trim(subelem_type), &
    & connec=subelem_glb_connec(j)%array, istat=istat, emsg=emsg)
    if (istat == STAT_FAILURE) then
      emsg = emsg//trim(msgloc)
      call clean_up (subelem_glb_connec, subelem_lcl_connec_top, &
      & subelem_lcl_connec_bot, x1, x2, xc)
      return
    end if

  end do


  ! deallocate local alloc arrays before successful return
  call clean_up (subelem_glb_connec, subelem_lcl_connec_top, &
  & subelem_lcl_connec_bot, x1, x2, xc)
  return



  contains


  pure subroutine clean_up (subelem_glb_connec, subelem_lcl_connec_top, &
  & subelem_lcl_connec_bot, x1, x2, xc)
  type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subelem_glb_connec(:)
  type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subelem_lcl_connec_top(:)
  type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subelem_lcl_connec_bot(:)
  real(DP),              allocatable, intent(inout) :: x1(:), x2(:), xc(:)

  if(allocated(subelem_glb_connec)) deallocate(subelem_glb_connec)
  if(allocated(subelem_lcl_connec_top)) deallocate(subelem_lcl_connec_top)
  if(allocated(subelem_lcl_connec_bot)) deallocate(subelem_lcl_connec_bot)
  if(allocated(x1)) deallocate(x1)
  if(allocated(x2)) deallocate(x2)
  if(allocated(xc)) deallocate(xc)
  end subroutine clean_up


end subroutine partition_element



pure subroutine update_internal_nodes (el, nodes)
! Purpose:
! - update the components of the internal nodes of this element
use parameter_module, only : ZERO, ONE, COH_CRACK_EDGE
use xnode_module,     only : xnode, operator(+), operator(*)

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

  do j = 1, NEDGE_SURF
    ! if this top edge is not cracked, cycle to the next
    if (el%top_edge_status(j) < COH_CRACK_EDGE) cycle
    ! extract the local index of the cracked top edge
    Icrackedge  = j
    ! find the internal node on the corresponding
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
                      & NDIM
use xnode_module,             only : xnode
use cohesive_material_module, only : cohesive_material
use baseCoh_element_module,   only : integrate
use global_toolkit_module,    only : assembleKF

  ! - passed in variables
  type(fCoh3d8_subelem),    intent(inout) :: elem
  type(xnode),              intent(in)    :: nodes(NNODE)
  type(cohesive_material),  intent(in)    :: material
  real(DP),                 intent(out)   :: K_matrix(NDOF,NDOF), F_vector(NDOF)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure

  ! - local variables
  character(len=MSGLENGTH) :: msgloc
  real(DP), allocatable :: Ki(:,:), Fi(:)
  integer,  allocatable :: subcnc(:)
  real(DP), allocatable :: Kmat_r(:,:), Fvec_r(:)
  real(DP), allocatable :: Tmatrixfull(:,:)
  logical :: nofail
  integer :: isub, j, k, l, n

  ! initialize intent out and local variables
  K_matrix = ZERO
  F_vector = ZERO
  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ', integrate_assemble_subelem, fCoh3d8_subelem_module'
  nofail = .false.
  isub=0; j=0; k=0; l=0; n=0

  if (present(nofailure)) nofail = nofailure


  !:::::::::::::::::::::::::::::::::::::!
  ! integrate and assemble sub elems
  !:::::::::::::::::::::::::::::::::::::!
  do isub = 1, size(elem%subelem)

    ! extract the local node connec of this sub elem
    if(allocated(subcnc)) deallocate(subcnc)
    allocate(subcnc(size(elem%subelem_lcl_connec(isub)%array)))
    subcnc = elem%subelem_lcl_connec(isub)%array

    ! call the integrate procedure asoc. with the subelem type
    ! exit the do loop (and if-else subsequently) if an error is encountered
    call integrate(elem%subelem(isub), nodes(subcnc), material, Ki, Fi, &
    & istat, emsg, nofail)
    if (istat == STAT_FAILURE) exit

    ! assemble the sub elem K and F
    call assembleKF(K_matrix, F_vector, Ki, Fi, subcnc, NDIM, istat, emsg)
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


  ! clean up if above loop exit upon error
  if (istat == STAT_FAILURE) then
    emsg = emsg//trim(msgloc)
    K_matrix = ZERO
    F_vector = ZERO
  end if

  if(allocated(Ki))             deallocate(Ki)
  if(allocated(Fi))             deallocate(Fi)
  if(allocated(subcnc))         deallocate(subcnc)
  if(allocated(Tmatrixfull))    deallocate(Tmatrixfull)
  if(allocated(Kmat_r))         deallocate(Kmat_r)
  if(allocated(Fvec_r))         deallocate(Fvec_r)

end subroutine integrate_assemble_subelem



pure subroutine get_Tmatrix (el, Tmatrixfull)
! Purpose:
! - calculate [Tmatrixfull], s.t. {U} = [Tmatrixfull] . {U_condensed}, where
! {U_condensed} is the nodal dof vector without the internal nodal dofs. The
! internal nodal dofs are interpolated from the other nodal dofs.
use parameter_module, only : DP, ZERO, ONE, COH_CRACK_EDGE

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
  do j = 1, NEDGE_SURF
    ! if this top edge is not cracked, cycle to the next
    if (el%top_edge_status(j) < COH_CRACK_EDGE) cycle
    ! extract the local index of the cracked top edge
    Icrackedge  = j
    ! find the internal node on the corresponding
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
