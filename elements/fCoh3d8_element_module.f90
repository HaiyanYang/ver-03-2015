module fCoh3d8_subelem_module
!
!  Purpose:
!    define a floating-node coh3d8 sub-element object, with a breakable top
!    surface. 
!
!  topological definition of this element, local nodal indices:
!
!  8__14____(E3)__13___7
!  |\                  |\
!  | \15               | \12
!  |  \(E4)            |  \(E2)
!  |   \16             |   \11
!  |    \              |    \
!  |     \5___9__(E1)__|_10__\6
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
!  bottom surface rl. nodes (counter-clock vise): 1, 2, 3, 4
!  top    surface rl. nodes (counter-clock vise): 5, 6, 7, 8
!  top    surface     edges (counter-clock vise): E1, E2, E3, E4
!
!  nodes of E1 ::    <end nodes: 5, 6>   <fl. nodes:  9, 10>
!  nodes of E2 ::    <end nodes: 6, 7>   <fl. nodes: 11, 12>
!  nodes of E3 ::    <end nodes: 7, 8>   <fl. nodes: 13, 14>
!  nodes of E4 ::    <end nodes: 8, 5>   <fl. nodes: 15, 16>
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    19/04/15  B. Y. Chen            Original code
!
use parameter_module,            only : NST=>NST_COHESIVE, NDIM, DP, MSGLENGTH,&
                                 & STAT_SUCCESS, STAT_FAILURE, ZERO
use global_clock_module,         only : program_clock, GLOBAL_CLOCK, clock_in_sync
use global_material_list_module, only : global_cohesive_list
use global_node_list_module,     only : global_node_list
use global_edge_list_module,     only : global_edge_list
use global_toolkit_module,       only : 
use baseCoh_element_module


implicit none

private

! list of private parameters in this module:
! NNDRL     : no. of real     nodes in this element
! NNDFL     : no. of floating nodes in this element
! NEDGE_TOP : no. of edges in this element on the top surface
! NNODE     : total no. of nodes in this element
! NDOF      : no. of degree-of-freedom in this element
! NODES_ON_TOP_EDGES(i,j) : this is a matrix to describe the nodal connectivity
!             of the top edges in this element. Each edge is topologically
!             composed of its two end nodes and two floating nodes. This 
!             term stores the local nodal index of the i_th node on edge j
integer, parameter :: NNDRL = 8, NEDGE_TOP = 4, NNDFL = 2 * NEDGE_TOP, &
                    & NNODE = NNDRL+NNDFL, NDOF = NDIM * NNODE
integer, parameter :: NODES_ON_TOP_EDGES(4,NEDGE_TOP) = &
         & reshape([5,6,9,10, 6,7,11,12, 7,8,13,14, 8,5,15,16], [4,NEDGE_TOP])
                                          


type, public :: fCoh3d8_subelem
  private
  ! list of components of this type:
  ! pstat               : elem partition status
  ! node_connec         : nodes global connectivity
  ! edge_connec         : edges global connectivity
  ! ID_matlist          : material list index
  ! lcl_ID_crack_edges  : local indices of cracked edges, this array is passed
  !                       from adj. ply elem. it is used to extract no. and ID
  !                       of cracked edges
  ! subelem             : list of subelems of type baseCoh
  ! subelem_lcl_connec  : lcl connec of subelems' nodes
  ! subelem_T           : T matrices of subelems' K and F
  integer :: pstat                         = 0
  integer :: node_connec(NNODE)            = 0
  integer :: edge_connec(NEDGE_TOP)        = 0
  integer :: ID_matlist                    = 0
  integer :: lcl_ID_crack_edges(NEDGE_TOP) = 0
  type(baseCoh_element), allocatable :: subelem(:)
  type(int_alloc_array), allocatable :: subelem_lcl_connec(:)
  type(real_alloc_matrix), allocatable :: subelem_T(:)
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

!~interface precrack
!~    module procedure precrack_fCoh3d8_subelem
!~end interface

interface integrate
    module procedure integrate_fCoh3d8_subelem
end interface

interface extract
    module procedure extract_fCoh3d8_subelem
end interface


public :: empty,set,update,integrate,extract




contains




pure subroutine empty_fCoh3d8_subelem (elem)

  type(fCoh3d8_subelem), intent(inout) :: elem
  
  type(fCoh3d8_subelem) :: elem_lcl
  
  elem = elem_lcl
    
end subroutine empty_fCoh3d8_subelem



pure subroutine set_fCoh3d8_subelem (elem, node_connec, edge_connec, ID_matlist,&
& istat, emsg)
! Purpose:
! to set the element ready for first use

  type(fCoh3d8_subelem),  intent(inout)   :: elem
  integer,                intent(in)      :: node_connec(NNODE)
  integer,                intent(in)      :: edge_connec(NEDGE_TOP)
  integer,                intent(in)      :: ID_matlist
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  
  ! local copy of elem
  type(fCoh3d8_subelem) :: elem_lcl
  ! global_connec of sub element
  integer, allocatable  :: global_connec(:)
  
  istat = STAT_SUCCESS
  emsg  = ''
  
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
  
  if ( ID_matlist < 1 ) then
    istat = STAT_FAILURE
    emsg  = 'ID_matlist must be >=1, set, &
    &fCoh3d8_subelem_module'
    return
  end if
  
  ! update to elem_lcl first
  elem_lcl%node_connec = node_connec
  elem_lcl%edge_connec = edge_connec
  elem_lcl%ID_matlist  = ID_matlist
  
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
  & ID_matlist=elem_lcl%ID_matlist, istat, emsg)
  
  ! if an error is encountered in set, clean up and exit program
  if (istat == STAT_FAILURE) then
    if (allocated(global_connec)) deallocate(global_connec)
    return
  end if

  ! update to dummy arg. elem before successful return
  elem = elem_lcl

end subroutine set_fCoh3d8_subelem



pure subroutine update_fCoh3d8_subelem (elem, lcl_ID_crack_edges, istat, emsg)
! Purpose:
! this subroutine is used solely to update the lcl_ID_crack_edges array of the element

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
& ID_matlist, lcl_ID_crack_edges, newpartition, local_clock, subelem, subelem_lcl_connec)
!** to be cleaned up **
! consider include dm, traction, separation, and remove lcl_ID_crack_edges, 
! newpartition, subelem, subelem_lcl_connec. extract only the useful components

  type(fCoh3d8_subelem),                    intent(in)  :: elem
  integer,                        optional, intent(out) :: pstat
  integer,           allocatable, optional, intent(out) :: node_connec(:)
  integer,           allocatable, optional, intent(out) :: edge_connec(:)
  integer,                        optional, intent(out) :: ID_matlist
  integer,           allocatable, optional, intent(out) :: lcl_ID_crack_edges(:)
  logical,                        optional, intent(out) :: newpartition
  type(program_clock),            optional, intent(out) :: local_clock
  type(sub3d_element),   allocatable, optional, intent(out) :: subelem(:)
  type(int_alloc_array), allocatable, optional, intent(out) :: subelem_lcl_connec(:)

  if(present(pstat)) pstat=elem%pstat
  
  if(present(node_connec)) then 
      allocate(node_connec(NNODE))
      node_connec=elem%node_connec
  end if
  
  if(present(edge_connec)) then 
      allocate(edge_connec(NEDGE_TOP))
      edge_connec=elem%edge_connec
  end if
  
  if(present(ID_matlist)) ID_matlist=elem%ID_matlist

  if(present(lcl_ID_crack_edges)) then 
      allocate(lcl_ID_crack_edges(NEDGE_TOP))
      lcl_ID_crack_edges=elem%lcl_ID_crack_edges
  end if
  
  if(present(newpartition)) newpartition=elem%newpartition
  
  if(present(local_clock)) local_clock=elem%local_clock
  
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



pure subroutine integrate_fCoh3d8_subelem (elem, K_matrix, F_vector, istat,&
& emsg, nofailure)
! Purpose:
! to integrate this fCoh3d8 subelem

  type(fCoh3d8_subelem),    intent(inout) :: elem 
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure

  ! local variables
  type(fCoh3d8_subelem) :: el
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
    esmg  = 'sub element is NOT yet allocated, integrate, fCoh3d8_subelem module'
    return
  end if
  
  ! define local variables
  
  ! copy pass arg. to its local copy
  el = elem
  
  ! by default, damage modelling is allowed, unless specified
  if(present(nofailure)) nofail = nofailure
    

  !**** MAIN CALCULATIONS ****
  !
  ! select what to do based on elem's pstat:
  ! if elem is at INTACT:
  !   - check edge status variables to see if the neighbour top ply has matrix 
  !     crack; if so, update pstat to PARTITIONED_FCOHSUB and partition sub elems.
  !   - integrate and assemble sub elems
  ! if elem is at PARTITIONED_FCOHSUB:
  !   - no need to check edge status variables; no need to update partition.
  !     just integrate and assemble sub elems
  ! other status of pstat variable is NOT allowed
  !
  select case (el%pstat)
  
    case (INTACT)
    ! elem is not yet partitioned
        ! check for edge status, update pstat and sub elem partition
        call edge_status_partition(el, istat, emsg)
        if (istat == STAT_SUCCESS) then
          ! integrate sub elems and assemble system matrices
          call integrate_assemble_subelem(el, K_matrix, F_vector, istat, emsg, nofail)
        end if
    
    case (PARTITIONED_FCOHSUB)
    ! elem is already partitioned
        ! just integrate sub elems and assemble system matrices
        call integrate_assemble_subelem(el, K_matrix, F_vector, istat, emsg, nofail) 
  
    case default
    ! this place should NOT be reached
        istat = STAT_FAILURE
        emsg  = 'unsupported pstat value in fCoh3d8 elem module'
    
  end select
  
  !**** END MAIN CALCULATIONS ****

  ! if above loop is exit with error, clean up and exit program
  if (istat == STAT_FAILURE) then
    K_matrix = ZERO
    F_vector = ZERO
    return
  end if
  
  ! update to dummy arg. elem before successful return
  elem = el

end subroutine integrate_fCoh3d8_subelem



pure subroutine edge_status_partition (el, istat, emsg)
! Purpose:
! update pass arg. el's pstat and its partition w.r.t its edge status variables

  ! passed-in variables
  type(fCoh3d8_subelem),    intent(inout) :: el
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variables
  ! edge_status  : edge status variables of the top edges of this elem, 
  !                extracted from global edge list
  ! n_crackedges : no. of cracked edges in this elem
  ! edge1, edge2 : indices of broken edges
  integer :: edge_status(NEDGE_TOP)
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
  edge_status  = 0
  n_crackedges = 0
  edge1 = 0
  edge2 = 0
  i = 0
  
  ! Note: no need to check for input validity or create local copy of intent  
  ! inout dummy arg because this is an internal procedure
  
  ! no need to update partition if elem is already in final partition
  if(el%pstat == PARTITIONED_FCOHSUB) return

  ! extract edge status variables from glb edge list
  edge_status(:) = global_edge_list(el%edge_connec(:))
  
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
      
  if (el%pstat == PARTITIONED_FCOHSUB) then
  ! update elem partitions into sub elems when it reaches 
  ! PARTITIONED_FCOHSUB status                
    call update_subelem_connec(el, istat, emsg)
  end if
  
  !**** END MAIN CALCULATIONS ****

end subroutine edge_status_partition



pure subroutine update_subelem_connec (el, istat, emsg)
! Purpose:
! to partition the element into sub elements
! specifically, this subroutine will set the following el's components:
! - subelem_lcl_connec
! - subelem_T
! - subelem
! according to the edge status of top edges of this element

  ! passed-in variables
  type(fCoh3d8_subelem),    intent(inout) :: el
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variables
  ! edge_status         : edge status variables of elem top edges
  ! subelem_glb_connec  : global connec of sub elem nodes
  ! ibe1, ibe2          : indices of broken edges
  ! e1, e2, e3, e4      : edge indices, used for partitioning element
  ! nsub                : no. of sub elements
  ! jcracknode1/2       : fl. node indices on 2 cracked top edges
  ! jendnode1/2         : end node indices of an edge
  ! x1, x2, xc          : coords of crack edge end nodes (x1 & x2) and 
  !                       crack tip (xc)
  ! tratio1,2           : ratio of length {xc-x1} over length {x2-x1},
  !                       tratio = |xc-x1|/|x2-x1|
  type(xnode)           :: nodes(NNODE)
  integer               :: edge_status(NEDGE_TOP)
  type(int_alloc_array), allocatable :: subelem_glb_connec(:)
  integer               :: ibe1, ibe2
  integer               :: e1, e2, e3, e4
  integer               :: nsub
  integer               :: jcracknode1, jcracknode2, jendnode1, jendnode2
  real(DP), allocatable :: x1(:), x2(:), xc(:)
  real(DP)              :: tratio1, tratio2
  ! counters
  integer :: i, j, l

  ! initialize intent out and local variables
  istat = STAT_SUCCESS
  emsg  = ''
  edge_status = 0
  ibe1 = 0
  ibe2 = 0
  e1 = 0
  e2 = 0
  e3 = 0
  e4 = 0
  nsub = 0
  jcracknode1 = 0
  jcracknode2 = 0
  jendnode1   = 0
  jendnode2   = 0
  tratio1 = ZERO
  tratio2 = ZERO
  i=0; j=0; l=0
  
  ! extract nodes from global node list
  nodes = global_node_list(el%node_connec)
    
  ! extract edge status variables from glb edge list
  edge_status(:) = global_edge_list(el%edge_connec(:))
  
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
        ibe1 = min( el%lcl_ID_crack_edges(1), el%lcl_ID_crack_edges(2) )
        ibe2 = max( el%lcl_ID_crack_edges(1), el%lcl_ID_crack_edges(2) ) 
        ! ibe1 must be between 1 to 3, and ibe2 between 2 to 4, with ibe2 > ibe1
        if(ibe1<1 .or. ibe1>3 .or. ibe2<2 .or. ibe2>4 .or. ibe2<=ibe1) then
          istat = STAT_FAILURE
          emsg = 'wrong crack edge indices in fCoh3d8 update subelem_lcl_connec &
          & case n_crackedges=2'
          return
        end if
        ! find nsub and e1 - e4
        ! nsub: no. of sub domains
        ! e1 - e4: re-index edges to facilitate partitioning domain
        select_ibe: select case(ibe1)
          case(1) select_ibe
              select case(ibe2)
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
          case(2) select_ibe
              select case(ibe2)
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
          case(3) select_ibe
              if(ibe2==4) then
                  nsub=4
                  e1=3; e2=4; e3=1; e4=2
              else
                  istat = STAT_FAILURE
                  emsg = 'wrong 2nd broken edge in update subelem_lcl_connec fCoh3d8'
                  return                   
              end if    
          case default select_ibe
              istat = STAT_FAILURE
              emsg = 'wrong broken edge in update subelem_lcl_connec fCoh3d8'
              return
        end select select_ibe
        !
        ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
        ! form sub elements based on nsub values
        ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
        select_nsub: select case(nsub)
          ! :::::::::::::::::::::::::::::::::::::::::!
          ! two quad subdomains, two coh3d8 sub elems
          ! :::::::::::::::::::::::::::::::::::::::::!
          case(2) select_nsub
              ! allocate nsub no. of subelem, subelem_lcl_connec and subelem_T
              if(allocated(el%subelem)) &
              & deallocate(el%subelem)
              if(allocated(el%subelem_lcl_connec)) &
              & deallocate(el%subelem_lcl_connec)
              if(allocated(el%subelem_T)) &
              & (el%subelem_T)
              if(allocated(subelem_glb_connec)) &
              & deallocate(subelem_glb_connec)
              allocate(el%subelem(nsub))
              allocate(el%subelem_lcl_connec(nsub)) 
              allocate(el%subelem_T(nsub))
              allocate(subelem_glb_connec(nsub))
              ! allocate&initialize the internal arrays of these arrays
              do j=1, nsub
                ! allocate connec for the 8 numerical nodes of coh3d8 sub elems
                ! allocate Tmatrix for 4 * 4 (bottom surf. 4 material nodes   
                ! interpolated by 4 numerical nodes)
                allocate(el%subelem_lcl_connec(j)%array(8))
                allocate(subelem_glb_connec(j)%array(8))
                allocate(el%subelem_T(j)%matrix(4,4))
                ! initialize these arrays
                el%subelem_lcl_connec(j)%array = 0
                subelem_glb_connec(j)%array    = 0
                el%subelem_T(j)%matrix         = ZERO
              end do
              ! find the rel. position of crack tip on e1 w.r.t its endnode 1
              jcracknode1 = NODES_ON_TOP_EDGES(3,e1)
              jendnode1   = NODES_ON_TOP_EDGES(1,e1)
              jendnode2   = NODES_ON_TOP_EDGES(2,e1)
              call extract(nodes(jendnode1),   x=x1)
              call extract(nodes(jendnode2),   x=x2)
              call extract(nodes(jcracknode1), x=xc)
              tratio1 = distance(x1,xc,NDIM) / distance(x1,x2,NDIM)
              ! find the rel. position of crack tip on e3 w.r.t its endnode 1
              jcracknode2 = NODES_ON_TOP_EDGES(3,e3)
              jendnode1   = NODES_ON_TOP_EDGES(1,e3)
              jendnode2   = NODES_ON_TOP_EDGES(2,e3)
              call extract(nodes(jendnode1),   x=x1)
              call extract(nodes(jendnode2),   x=x2)
              call extract(nodes(jcracknode2), x=xc)
              tratio2 = distance(x1,xc,NDIM) / distance(x1,x2,NDIM)
              !:::::::::::::::::::::::::!
              !*** define sub elm 1 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 4 numerical nodes: all real nodes
              ! top 4 numerical nodes: e1 nodes 1 & 3, and e3 nodes 4 & 2 
              el%subelem_lcl_connec(1)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2
              el%subelem_lcl_connec(1)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
              el%subelem_lcl_connec(1)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
              el%subelem_lcl_connec(1)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
              el%subelem_lcl_connec(1)%array(5)=NODES_ON_TOP_EDGES(1,e1)
              el%subelem_lcl_connec(1)%array(6)=NODES_ON_TOP_EDGES(3,e1)
              el%subelem_lcl_connec(1)%array(7)=NODES_ON_TOP_EDGES(4,e3)
              el%subelem_lcl_connec(1)%array(8)=NODES_ON_TOP_EDGES(2,e3)
              ! define its global connec with global node list
              subelem_glb_connec(1)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(1)%array(:))
              ! update T matrix, 
              ! interpolate bottom surf material nodes with numerical nodes 1-4
              el%subelem_T(1)%matrix = ZERO
              el%subelem_T(1)%matrix(1,1)=ONE
              el%subelem_T(1)%matrix(2,1)=ONE-tratio1
              el%subelem_T(1)%matrix(2,2)=tratio1
              el%subelem_T(1)%matrix(3,3)=ONE-tratio2
              el%subelem_T(1)%matrix(3,4)=tratio2
              el%subelem_T(1)%matrix(4,4)=ONE
              ! set sub element 1
              call set(el%subelem(1), eltype='coh3d8', &
              & connec=subelem_glb_connec(1)%array, ID_matlist=el%ID_matlist, &
              & istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
              !:::::::::::::::::::::::::!
              !*** define sub elm 2 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 4 numerical nodes: all real nodes
              ! top 4 numerical nodes: e1 nodes 4 & 2, and e3 nodes 1 & 3 
              el%subelem_lcl_connec(2)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2
              el%subelem_lcl_connec(2)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
              el%subelem_lcl_connec(2)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
              el%subelem_lcl_connec(2)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
              el%subelem_lcl_connec(2)%array(5)=NODES_ON_TOP_EDGES(4,e1)
              el%subelem_lcl_connec(2)%array(6)=NODES_ON_TOP_EDGES(2,e1)
              el%subelem_lcl_connec(2)%array(7)=NODES_ON_TOP_EDGES(1,e3)
              el%subelem_lcl_connec(2)%array(8)=NODES_ON_TOP_EDGES(3,e3)
              ! define its global connec with global node list
              subelem_glb_connec(2)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(2)%array(:))
              ! update T matrix, 
              ! interpolate bottom surf material nodes with numerical nodes 1-4
              el%subelem_T(2)%matrix=ZERO
              el%subelem_T(2)%matrix(1,1)=ONE-tratio1
              el%subelem_T(2)%matrix(1,2)=tratio1
              el%subelem_T(2)%matrix(2,2)=ONE
              el%subelem_T(2)%matrix(3,3)=ONE
              el%subelem_T(2)%matrix(4,3)=ONE-tratio2
              el%subelem_T(2)%matrix(4,4)=tratio2
              ! set sub element 2
              call set(el%subelem(2), eltype='coh3d8', &
              & connec=subelem_glb_connec(2)%array, ID_matlist=el%ID_matlist, &
              & istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
          ! :::::::::::::::::::::::::::::::::::::::::!
          ! four tri subdomains, four coh3d6 sub elems
          ! :::::::::::::::::::::::::::::::::::::::::!
          case(4) select_nsub
              ! allocate nsub no. of subelem, subelem_lcl_connec and subelem_T
              if(allocated(el%subelem)) &
              & deallocate(el%subelem)
              if(allocated(el%subelem_lcl_connec)) &
              & deallocate(el%subelem_lcl_connec)
              if(allocated(el%subelem_T)) &
              & deallocate(el%subelem_T)
              if(allocated(subelem_glb_connec)) &
              & deallocate(subelem_glb_connec)
              allocate(el%subelem(nsub))
              allocate(el%subelem_lcl_connec(nsub))
              allocate(el%subelem_T(nsub))
              allocate(subelem_glb_connec(nsub))
              do j=1, nsub
                ! allocate 6 numerical nodes for coh3d6 sub elems
                ! allocate Tmatrix for 3 * 4 (bottom surf. 3 material nodes  
                ! interpolated by 4 numerical nodes)
                allocate(el%subelem_lcl_connec(j)%array(6))
                allocate(subelem_glb_connec(j)%array(6))
                allocate(el%subelem_T(j)%matrix(3,4))
                ! initialize these arrays
                el%subelem_lcl_connec(j)%array = 0
                subelem_glb_connec(j)%array    = 0
                el%subelem_T(j)%matrix         = ZERO
              end do
              ! find the rel. position of crack tip on e1 w.r.t its endnode 1
              jcracknode1 = NODES_ON_TOP_EDGES(3,e1)
              jendnode1   = NODES_ON_TOP_EDGES(1,e1)
              jendnode2   = NODES_ON_TOP_EDGES(2,e1)
              call extract(nodes(jendnode1),   x=x1)
              call extract(nodes(jendnode2),   x=x2)
              call extract(nodes(jcracknode1), x=xc)
              tratio1 = distance(x1,xc,NDIM) / distance(x1,x2,NDIM)
              ! find the rel. position of crack tip on e2 w.r.t its endnode 1
              jcracknode2 = NODES_ON_TOP_EDGES(3,e2)
              jendnode1   = NODES_ON_TOP_EDGES(1,e2)
              jendnode2   = NODES_ON_TOP_EDGES(2,e2)
              call extract(nodes(jendnode1),   x=x1)
              call extract(nodes(jendnode2),   x=x2)
              call extract(nodes(jcracknode2), x=xc)
              tratio2 = distance(x1,xc,NDIM) / distance(x1,x2,NDIM)
              !:::::::::::::::::::::::::!
              !*** define sub elm 1 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 4 numerical nodes: all real nodes
              ! top 3 numerical nodes: e1 nodes 4, and e2 nodes 1 & 3
              el%subelem_lcl_connec(1)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2
              el%subelem_lcl_connec(1)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
              el%subelem_lcl_connec(1)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
              el%subelem_lcl_connec(1)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
              el%subelem_lcl_connec(1)%array(5)=NODES_ON_TOP_EDGES(4,e1)
              el%subelem_lcl_connec(1)%array(6)=NODES_ON_TOP_EDGES(1,e2)
              el%subelem_lcl_connec(1)%array(7)=NODES_ON_TOP_EDGES(3,e2)
              subelem_glb_connec(1)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(1)%array(:))
              ! update Tmatrix, 
              ! interpolate bottom surf 3 mat. nodes with numerical nodes 1-4
              el%subelem_T(1)%matrix=ZERO
              el%subelem_T(1)%matrix(1,1)=ONE-tratio1
              el%subelem_T(1)%matrix(1,2)=tratio1
              el%subelem_T(1)%matrix(2,2)=ONE
              el%subelem_T(1)%matrix(3,2)=ONE-tratio2
              el%subelem_T(1)%matrix(3,3)=tratio2
              ! set sub elem 1
              call set(el%subelem(1), eltype='coh3d6', &
              & connec=subelem_glb_connec(1)%array, ID_matlist=el%ID_matlist,&
              & istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
              !:::::::::::::::::::::::::!
              !*** define sub elm 2 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 4 numerical nodes: all real nodes
              ! top 3 numerical nodes: e2 node 4, and e3 nodes 1 & 2
              el%subelem_lcl_connec(2)%array(1)=NODES_ON_TOP_EDGES(1,e2)-NNDRL/2
              el%subelem_lcl_connec(2)%array(2)=NODES_ON_TOP_EDGES(2,e2)-NNDRL/2
              el%subelem_lcl_connec(2)%array(3)=NODES_ON_TOP_EDGES(1,e4)-NNDRL/2
              el%subelem_lcl_connec(2)%array(4)=NODES_ON_TOP_EDGES(2,e4)-NNDRL/2
              el%subelem_lcl_connec(2)%array(5)=NODES_ON_TOP_EDGES(4,e2)
              el%subelem_lcl_connec(2)%array(6)=NODES_ON_TOP_EDGES(1,e3)
              el%subelem_lcl_connec(2)%array(7)=NODES_ON_TOP_EDGES(2,e3) 
              subelem_glb_connec(2)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(2)%array(:)) 
              ! update Tmatrix, 
              ! interpolate bottom surf 3 mat. nodes with numerical nodes 1-4
              el%subelem_T(2)%matrix=ZERO
              el%subelem_T(2)%matrix(1,1)=ONE-tratio2
              el%subelem_T(2)%matrix(1,2)=tratio2
              el%subelem_T(2)%matrix(2,2)=ONE
              el%subelem_T(2)%matrix(3,3)=ONE
              ! set sub elem 2
              call set(el%subelem(2),eltype='coh3d6', &
              & connec=subelem_glb_connec(2)%array, ID_matlist=el%ID_matlist, &
              & istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
              !:::::::::::::::::::::::::!
              !*** define sub elm 3 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 4 numerical nodes: all real nodes
              ! top 3 numerical nodes: e4 nodes 1 & 2, and e1 node 3
              el%subelem_lcl_connec(3)%array(1)=NODES_ON_TOP_EDGES(1,e4)-NNDRL/2
              el%subelem_lcl_connec(3)%array(2)=NODES_ON_TOP_EDGES(2,e4)-NNDRL/2
              el%subelem_lcl_connec(3)%array(3)=NODES_ON_TOP_EDGES(1,e2)-NNDRL/2
              el%subelem_lcl_connec(3)%array(4)=NODES_ON_TOP_EDGES(2,e2)-NNDRL/2
              el%subelem_lcl_connec(3)%array(5)=NODES_ON_TOP_EDGES(1,e4)
              el%subelem_lcl_connec(3)%array(6)=NODES_ON_TOP_EDGES(2,e4)
              el%subelem_lcl_connec(3)%array(7)=NODES_ON_TOP_EDGES(3,e1)
              subelem_glb_connec(3)%array(:) = & 
              & el%node_connec(el%subelem_lcl_connec(3)%array(:))
              ! update Tmatrix, 
              ! interpolate bottom surf 3 mat. nodes with numerical nodes 1-4
              el%subelem_T(3)%matrix=ZERO
              el%subelem_T(3)%matrix(1,1)=ONE
              el%subelem_T(3)%matrix(2,2)=ONE
              el%subelem_T(3)%matrix(3,2)=ONE-tratio1
              el%subelem_T(3)%matrix(3,3)=tratio1
              ! set sub elem 3
              call set(el%subelem(3),eltype='coh3d6', &
              & connec=subelem_glb_connec(3)%array, ID_matlist=el%ID_matlist, &
              & istat=istat, emsg=emsg)
              if (istat == STAT_FAILURE) then
                call clean_up (subelem_glb_connec, x1, x2, xc)
                return
              end if
              !:::::::::::::::::::::::::!
              !*** define sub elm 4 ***
              !:::::::::::::::::::::::::!
              ! define its local connec with parent element nodes
              ! bot 4 numerical nodes: all real nodes
              ! top 3 numerical nodes: e1 node 3, e2 node 4 and e3 node 2
              el%subelem_lcl_connec(4)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2
              el%subelem_lcl_connec(4)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
              el%subelem_lcl_connec(4)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
              el%subelem_lcl_connec(4)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
              el%subelem_lcl_connec(4)%array(5)=NODES_ON_TOP_EDGES(3,e1)
              el%subelem_lcl_connec(4)%array(6)=NODES_ON_TOP_EDGES(4,e2)
              el%subelem_lcl_connec(4)%array(7)=NODES_ON_TOP_EDGES(2,e3)             
              subelem_glb_connec(4)%array(:) = &
              & el%node_connec(el%subelem_lcl_connec(4)%array(:))
              ! update Tmatrix, 
              ! interpolate bottom surf 3 mat. nodes with numerical nodes 1-4
              el%subelem_T(4)%matrix=ZERO
              el%subelem_T(4)%matrix(1,1)=ONE-tratio1
              el%subelem_T(4)%matrix(1,2)=tratio1
              el%subelem_T(4)%matrix(2,2)=ONE-tratio2
              el%subelem_T(4)%matrix(2,3)=tratio2
              el%subelem_T(4)%matrix(3,4)=ONE
              ! set sub elem 4
              call set(el%subelem(4),eltype='coh3d6', &
              & connec=subelem_glb_connec(4)%array, ID_matlist=el%ID_matlist, &
              & istat=istat, emsg=emsg)
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
        type(int_alloc_array), allocatable :: subelem_glb_connec(:)
        real(DP),              allocatable :: x1(:), x2(:), xc(:)
        if(allocated(subelem_glb_connec)) deallocate(subelem_glb_connec)
        if(allocated(x1)) deallocate(x1)
        if(allocated(x2)) deallocate(x2)
        if(allocated(xc)) deallocate(xc)
      end subroutine clean_up


end subroutine update_subelem_connec



pure subroutine integrate_assemble_subelem (elem, K_matrix, F_vector, istat, &
& emsg, nofailure)
! Purpose:
! integrate and assemble sub element system arrays

  ! - passed in variables
  type(fCoh3d8_subelem),    intent(inout) :: elem
  real(DP), 	              intent(inout) :: K_matrix(NDOF,NDOF), F_vector(NDOF)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure
  
  ! - local variables
  real(DP),	allocatable :: Ki(:,:), Fi(:)
  integer,  allocatable :: dofcnc(:)
  real(DP), allocatable :: subKmat(:,:), subFvec(:)
  real(DP), allocatable :: Tmatrix(:,:), Tmatrixfull(:,:)
  type(xnode_alloc_array), allocatable :: subelem_mnodes(:)
  character(len=ELTYPELENGTH) :: subeltype
  logical :: nofail
  integer :: isub, j, k, l
    
  istat  = STAT_SUCCESS
  emsg   = ''
  subeltype = ''
  nofail = .false.
  isub=0; j=0; k=0; l=0
  
  ! empty K and F
  K_matrix = ZERO
  F_vector = ZERO
  
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
        call integrate(elem%subelem(isub), Ki, Fi, istat, emsg, nofail)
        if (istat == STAT_FAILURE) exit
        
        ! if no error, prepare to assemble Ki and Fi to elem's K and F
        if(allocated(dofcnc)) deallocate(dofcnc)
        allocate(dofcnc(size(Fi))); dofcnc = 0 
        
        ! loop over no. of nodes in sub elem i
        do j=1, size(elem%subelem_lcl_connec(isub)%array) 
          do l=1, NDIM
            ! obtain dof indices of the jth node of sub elem i 
            dofcnc((j-1)*NDIM+l)=(elem%subelem_lcl_connec(isub)%array(j)-1)*NDIM+l
          end do
        end do 
        
        call assembleKF(K_matrix, F_vector, Ki, Fi, dofcnc, istat, emsg)
        if (istat == STAT_FAILURE) exit
 
      end do
    
  else if (elem%pstat == PARTITIONED_FCOHSUB) then
  ! element is partitioned into subelems, update mnodes of each subelem and use
  ! them for subelem integration
  
      ! update the material nodes' variables of all sub elems
      call update_mnodes(elem, subelem_mnodes)
      
      ! integrate each sub elem and assemble together
      do isub = 1, size(elem%subelem)
      
        ! call the integrate procedure asoc. with the subelem type, use mnodes
        ! exit the do loop (and if-else subsequently) if an error is encountered
        call integrate(elem%subelem(isub), Ki, Fi, istat, emsg, nofail, &
        & mnodes=subelem_mnodes(isub)%array)
        if (istat == STAT_FAILURE) exit
        
        ! if no error, apply Tmatrix transform to Ki and Fi based on sub elem type 
        call extract(elem%subelem(isub), eltype=subeltype)
        
        select case (trim(subeltype))
          case('coh3d6')
              ! extract stored Tmatrix
              if(allocated(Tmatrix)) deallocate(Tmatrix)
              ! Tmatrix relating bottom surf. 3 mat. nodes with 4 num. nodes
              allocate(Tmatrix(3,4))
              Tmatrix = elem%subelem_T(isub)%matrix
              if(allocated(Tmatrixfull)) deallocate(Tmatrixfull)
              ! full Tmatrix relating 6 mat. nodes with 7 num. nodes
              allocate(Tmatrixfull(6*NDIM,7*NDIM))
              Tmatrixfull=ZERO
              ! Tmatrix corresponding to bottom 3 mat. nodes and 4 num. nodes
              do k=1, 3    
                  do j=1, 4
                      do l=1, NDIM
                          Tmatrixfull((k-1)*NDIM+l,(j-1)*NDIM+l)=Tmatrix(k,j)
                      end do
                  end do   
              end do
              ! Tmatrix corresponding to top 3 mat. nodes and 3 num. nodes
              do k=4, 6    
                  do l=1, NDIM
                      Tmatrixfull((k-1)*NDIM+l,k*NDIM+l)=ONE
                  end do  
              end do
              ! prepare sub elem's K matrix for update
              if(allocated(subKmat)) deallocate(subKmat)
              ! this sub elem's K matrix is for 7 num. nodes
              allocate(subKmat(7*NDIM,7*NDIM)); subKmat=zero
              ! prepare sub elem's F vector for update
              if(allocated(subFvec)) deallocate(subFvec)
              ! this sub elem's F vector is for 7 num. nodes
              allocate(subFvec(7*NDIM)); subFvec=zero
              ! calculate sub elem's K matrix and F vector
              subKmat=matmul(matmul(transpose(Tmatrixfull),Ki),Tmatrixfull)
              subFvec=matmul(transpose(Tmatrixfull),Fi)
            
          case('coh3d8')
              ! extract stored Tmatrix
              if(allocated(Tmatrix)) deallocate(Tmatrix)
              ! Tmatrix relating bottom 4 mat. nodes with 4 num. nodes
              allocate(Tmatrix(4,4))
              Tmatrix = elem%subelem_T(isub)%matrix
              if(allocated(Tmatrixfull)) deallocate(Tmatrixfull)
              ! full Tmatrix relating 8 mat. nodes with 8 num. nodes
              allocate(Tmatrixfull(8*NDIM,8*NDIM))
              Tmatrixfull=zero
              ! Tmatrix corresponding to bottom 4 mat nodes and 4 num nodes
              do k=1, 4
                  do j=1, 4
                      do l=1, NDIM
                          Tmatrixfull((k-1)*NDIM+l,(j-1)*NDIM+l)=Tmatrix(k,j)
                      end do
                  end do   
              end do
              ! Tmatrix corresponding to top 4 mat nodes and 4 num nodes
              do k=5, 8
                  do l=1, NDIM
                      Tmatrixfull((k-1)*NDIM+l,(k-1)*NDIM+l)=ONE
                  end do  
              end do
              ! prepare sub elem's K matrix for update
              if(allocated(subKmat)) deallocate(subKmat)
              ! this sub elem's K matrix is for 8 num. nodes
              allocate(subKmat(8*NDIM,8*NDIM)); subKmat=zero
              ! prepare sub elem's F vector for update
              if(allocated(subFvec)) deallocate(subFvec)
              ! this sub elem's F vector is for 8 num. nodes
              allocate(subFvec(8*NDIM)); subFvec=zero
              ! calculate sub elem's K matrix and F vector
              subKmat=matmul(matmul(transpose(Tmatrixfull),Ki),Tmatrixfull)
              subFvec=matmul(transpose(Tmatrixfull),Fi)

          case default
              istat = STAT_FAILURE
              emsg = 'unsupported sub elem type in integrate_assemble_subelem, &
              & fCoh3d8_subelem_module'
              exit
              
        end select
        
        ! prepare to assemble sub K and sub F to elem's K and F
        if(allocated(dofcnc)) deallocate(dofcnc)
        allocate(dofcnc(size(subFvec))); dofcnc = 0 
        ! loop over no. of nodes in sub elem i
        do j=1, size(elem%subelem_lcl_connec(isub)%array) 
          do l=1, NDIM
            ! dof indices of the jth node of sub elem i 
            dofcnc((j-1)*NDIM+l) = (elem%subelem_lcl_connec(isub)%array(j)-1)*NDIM+l
          end do
        end do 
        call assembleKF(K_matrix, F_vector, subKmat, subFvec, dofcnc, istat, emsg)
        if (istat == STAT_FAILURE) exit
        
      end do
    
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
  if(allocated(subKmat))        deallocate(subKmat)
  if(allocated(subFvec))        deallocate(subFvec)
  if(allocated(Tmatrix))        deallocate(Tmatrix)
  if(allocated(Tmatrixfull))    deallocate(Tmatrixfull)
  if(allocated(subelem_mnodes)) deallocate(subelem_mnodes)
    
end subroutine integrate_assemble_subelem



pure subroutine update_mnodes (elem, subelem_mnodes)

  ! passed-in variables
  type(fCoh3d8_subelem),  intent(in)    :: elem
  integer,                intent(in)    :: i
  type(xnode_alloc_array), allocatable, intent(out) :: subelem_mnodes(:)
  
  
  ! local variables
  type(int_alloc_array), allocatable :: subelem_glb_connec(:)  ! glb cnc of sub elements
  integer :: i, j, l                                  ! counters
  integer :: nsub                             		! no. of sub elements


  real(DP)                :: tratio1, tratio2 		! tratio=|xc-x1|/|x2-x1|
  real(DP), allocatable   :: Tmatrix(:,:)             ! interpolation matrix btw bottom num nodes and mat nodes

!       initialize local variables

    i=0; j=0; l=0
    nsub=0
    tratio1=ZERO; tratio2=ZERO
     

    ! currently only support elem partition PARTITIONED_FCOHSUB
    select case(elem%pstat)
    
      case(PARTITIONED_FCOHSUB)
          nsub=size(elem%subelem)
          
          select case(nsub)
              case(2)
              ! two coh3d8 sub elems

                  allocate(subelem_glb_connec(nsub))
                  
                  ! allocate 8 numerical nodes for sub elems
                  do j=1, nsub
                      allocate(subelem_glb_connec(j)%array(8))
                      subelem_glb_connec(j)%array=0
                      subelem_glb_connec(j)%array(:)=elem%node_connec(elem%subelem_lcl_connec(j)%array(:))
                  end do
                  
                  !*** sub elem 1 mnodes update
                  
                  Tmatrix = elem%subelem_T(1)%matrix
                  
                  ! extract ratio values from Tmatrix
                  tratio1=Tmatrix(2,1)
                  tratio2=Tmatrix(3,3)
                  
                  ! update the mnodes array
                  mnodes(1)=global_node_list(subelem_glb_connec(1)%array(1))
                  mnodes(2)=tratio1*global_node_list(subelem_glb_connec(1)%array(1))+(ONE-tratio1)*global_node_list(subelem_glb_connec(1)%array(2))
                  mnodes(3)=tratio2*global_node_list(subelem_glb_connec(1)%array(3))+(ONE-tratio2)*global_node_list(subelem_glb_connec(1)%array(4))
                  mnodes(4)=global_node_list(subelem_glb_connec(1)%array(4))
                  mnodes(5)=global_node_list(subelem_glb_connec(1)%array(5))
                  mnodes(6)=global_node_list(subelem_glb_connec(1)%array(6))
                  mnodes(7)=global_node_list(subelem_glb_connec(1)%array(7))
                  mnodes(8)=global_node_list(subelem_glb_connec(1)%array(8))                  
                  
                  subelem_mnodes(1)%array=mnodes
                                   
                  
                  !*** sub elm 2 mnodes update
                  
                  Tmatrix = elem%subelem_T(2)%matrix 
                  
                  tratio1=Tmatrix(1,1)
                  tratio2=Tmatrix(4,3)
                  
                  ! update the mnodes array
                  mnodes(1)=tratio1*global_node_list(subelem_glb_connec(2)%array(1))+(ONE-tratio1)*global_node_list(subelem_glb_connec(2)%array(2))
                  mnodes(2)=global_node_list(subelem_glb_connec(2)%array(2))
                  mnodes(3)=global_node_list(subelem_glb_connec(2)%array(3))
                  mnodes(4)=tratio2*global_node_list(subelem_glb_connec(2)%array(3))+(ONE-tratio2)*global_node_list(subelem_glb_connec(2)%array(4))
                  mnodes(5)=global_node_list(subelem_glb_connec(2)%array(5))
                  mnodes(6)=global_node_list(subelem_glb_connec(2)%array(6))
                  mnodes(7)=global_node_list(subelem_glb_connec(2)%array(7))
                  mnodes(8)=global_node_list(subelem_glb_connec(2)%array(8))
                  
                 subelem_mnodes(2)%array=mnodes                
                  
                  
                  
              case(4)
              ! four coh3d6 sub elems


                  if(allocated(subelem_glb_connec)) deallocate(subelem_glb_connec)
                  allocate(subelem_glb_connec(nsub))
                  do j=1, nsub   ! coh3d6 elem cnc
                      allocate(subelem_glb_connec(j)%array(7))
                      subelem_glb_connec(j)%array=0
                      subelem_glb_connec(j)%array(:)=elem%node_connec(elem%subelem_lcl_connec(j)%array(:))
                  end do

                              
                  
                  !*** sub elm 1 connec
                  
                  Tmatrix = elem%subelem_T(1)%matrix
                  
                  tratio1=Tmatrix(1,1)
                  tratio2=Tmatrix(3,2)
                  
                  mnodes(1)=tratio1*global_node_list(subelem_glb_connec(1)%array(1))+(ONE-tratio1)*global_node_list(subelem_glb_connec(1)%array(2))
                  mnodes(2)=global_node_list(subelem_glb_connec(1)%array(2))
                  mnodes(3)=tratio2*global_node_list(subelem_glb_connec(1)%array(2))+(ONE-tratio2)*global_node_list(subelem_glb_connec(1)%array(3))
                  mnodes(4)=global_node_list(subelem_glb_connec(1)%array(5))
                  mnodes(5)=global_node_list(subelem_glb_connec(1)%array(6))
                  mnodes(6)=global_node_list(subelem_glb_connec(1)%array(7))
                  
                  subelem_mnodes(1)%array=mnodes
                  

                  
                  !*** sub elm 2 connec
                  
                  Tmatrix = elem%subelem_T(2)%matrix
                  
                  tratio2=Tmatrix(1,1)
                  
                  mnodes(1)=tratio2*global_node_list(subelem_glb_connec(2)%array(1))+(ONE-tratio2)*global_node_list(subelem_glb_connec(2)%array(2))
                  mnodes(2)=global_node_list(subelem_glb_connec(2)%array(2))
                  mnodes(3)=global_node_list(subelem_glb_connec(2)%array(3))
                  mnodes(4)=global_node_list(subelem_glb_connec(2)%array(5))
                  mnodes(5)=global_node_list(subelem_glb_connec(2)%array(6))
                  mnodes(6)=global_node_list(subelem_glb_connec(2)%array(7))        
                  
                  subelem_mnodes(2)%array=mnodes

                  
                  
                  
                  !*** sub elm 3 connec
                  
                  Tmatrix = elem%subelem_T(3)%matrix
                  
                  tratio1=Tmatrix(3,2)
                  
                  mnodes(1)=global_node_list(subelem_glb_connec(3)%array(1))
                  mnodes(2)=global_node_list(subelem_glb_connec(3)%array(2))
                  mnodes(3)=tratio1*global_node_list(subelem_glb_connec(3)%array(2))+(ONE-tratio1)*global_node_list(subelem_glb_connec(3)%array(3))
                  mnodes(4)=global_node_list(subelem_glb_connec(3)%array(5))
                  mnodes(5)=global_node_list(subelem_glb_connec(3)%array(6))
                  mnodes(6)=global_node_list(subelem_glb_connec(3)%array(7))
    
                  subelem_mnodes(3)%array=mnodes
                  
                  
                  
                  !*** sub elm 4 connec
                  
                  Tmatrix = elem%subelem_T(4)%matrix
                  
                  tratio1=Tmatrix(1,1)
                  tratio2=Tmatrix(2,2)
                  
                  mnodes(1)=tratio1*global_node_list(subelem_glb_connec(4)%array(1))+(ONE-tratio1)*global_node_list(subelem_glb_connec(4)%array(2))
                  mnodes(2)=tratio2*global_node_list(subelem_glb_connec(4)%array(2))+(ONE-tratio2)*global_node_list(subelem_glb_connec(4)%array(3))
                  mnodes(3)=global_node_list(subelem_glb_connec(4)%array(4))
                  mnodes(4)=global_node_list(subelem_glb_connec(4)%array(5))
                  mnodes(5)=global_node_list(subelem_glb_connec(4)%array(6))
                  mnodes(6)=global_node_list(subelem_glb_connec(4)%array(7))          
                  
                  subelem_mnodes(4)%array=mnodes                   
                  
                  
              case default
                  write(msg_file,*)'wrong nsub in update subelem_lcl_connec fCoh3d8'
                  call exit_function
          end select

      case default
         !'partition status not supported in fCoh3d8 update mnodes!'
          
    end select
        
    ! deallocate local array
    
    if(allocated(subelem_glb_connec)) deallocate(subelem_glb_connec)
    if(allocated(mnodes)) deallocate(mnodes)
    if(allocated(Tmatrix)) deallocate(Tmatrix)


end subroutine update_mnodes



pure subroutine update_sdv (elem)

  ! passed-in variables
  type(fCoh3d8_subelem),    intent(inout)   :: elem
  
  ! local variables   
    character(len=eltypelength) ::  subeltype

    ! real and integer failure variables
    real(DP) :: rfvar
    integer :: ifvar
    
    ! coh elem arrays
    type(coh3d6_element),allocatable :: subcoh3d6(:)
    type(coh3d8_element),allocatable :: subcoh3d8(:)
    
    ! intg point arrays
    type(integration_point), allocatable :: igpnt(:) ! intg point array
    
    ! failure variables extracted from ig point sdv array
    type(sdv_array),allocatable :: fsdv(:) 
    
    integer :: i,j,l
    
    
    ! initialize local variables
    subeltype=''
    rfvar=ZERO; ifvar=0
    i=0; j=0; l=0
    

    ! update fCoh3d8 elem sdv for output
    if(.not.allocated(elem%sdv)) then
        allocate(elem%sdv(1))
        allocate(elem%sdv(1)%i(1))  ! pstat
        allocate(elem%sdv(1)%r(1))  ! dm
        elem%sdv(1)%i(1)=0
        elem%sdv(1)%r(1)=0
    end if
                
                
    do i=1,size(elem%subelem)           
    
        ifvar=0
        rfvar=ZERO
        
        ! extract this subelem type
        call extract(elem%subelem(i),eltype=subeltype)
        
        ! extract this subelem intg points based on subelem type
        select case(subeltype)                       
            case('coh3d6')
                call extract(elem%subelem(i),coh3d6=subcoh3d6)
                call extract(subcoh3d6(1),ig_point=igpnt)                        
            case('coh3d8')
                call extract(elem%subelem(i),coh3d8=subcoh3d8)
                call extract(subcoh3d8(1),ig_point=igpnt)                       
            case default
                write(msg_file,*)'subelem type not supported in fCoh3d8 elem!'
                call exit_function
        end select  
        
        ! sum up useful sdv values of all intg points into ifvar and rfvar
        do j=1,size(igpnt)
            call extract(igpnt(j),sdv=fsdv)
            if(allocated(fsdv)) then
                ! update sdv values (equilibrium sdv values, stored in sdv(1))
                if(allocated(fsdv(1)%i)) ifvar=ifvar+fsdv(1)%i(1)
                if(allocated(fsdv(1)%r)) rfvar=rfvar+fsdv(1)%r(1)
                deallocate(fsdv)
            end if
            
        end do 
        ! average fvar values in this sub element
        ifvar=int(ifvar/size(igpnt))
        rfvar=rfvar/size(igpnt)
        
        ! update this fvar to fCoh3d8 sdv
        elem%sdv(1)%i(1)=max(elem%sdv(1)%i(1),ifvar)
        elem%sdv(1)%r(1)=max(elem%sdv(1)%r(1),rfvar)
        
        ! deallocate arrays for re-use
        if(allocated(igpnt))     deallocate(igpnt)
        if(allocated(fsdv))      deallocate(fsdv)
        if(allocated(subcoh3d6)) deallocate(subcoh3d6)
        if(allocated(subcoh3d8)) deallocate(subcoh3d8)

    end do
    
    
    if(allocated(igpnt))     deallocate(igpnt)
    if(allocated(fsdv))      deallocate(fsdv)
    if(allocated(subcoh3d6)) deallocate(subcoh3d6)
    if(allocated(subcoh3d8)) deallocate(subcoh3d8)

end subroutine update_sdv




end module fCoh3d8_subelem_module
