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
  integer :: fstat                  = 0
  integer :: node_connec(NNODE)     = 0
  integer :: edge_connec(NEDGE_TOP) = 0
  integer :: ID_matlist             = 0
  integer :: crack_edges(NEDGE_TOP) = 0
  logical :: newpartition           = .false.
  type(program_clock)                :: local_clock
  type(baseCoh_element), allocatable :: subelem(:)
  type(int_alloc_array), allocatable :: subelem_connec(:)
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
  allocate(elem_lcl%subelem_connec(1))
  allocate(elem_lcl%subelem_connec(1)%array(NNDRL))
  allocate(global_connec(NNDRL))
  
  ! sub elem 1 local connec
  elem_lcl%subelem_connec(1)%array = [(i, i=1,NNDRL)]
  
  ! sub elem 1 global connec
  global_connec(:) = elem_lcl%node_connec( elem_lcl%subelem_connec(1)%array(:) )
  
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



pure subroutine update_fCoh3d8_subelem (elem, crack_edges, istat, emsg)
! Purpose:
! this subroutine is used solely to update the crack_edges array of the element

  type(fCoh3d8_subelem),    intent(inout) :: elem
  integer,                  intent(in)    :: crack_edges(NEDGE_TOP)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  istat = STAT_SUCCESS
  emsg  = ''
  
  ! check validity of inputs
  if ( any(crack_edges < 1) .or. any(crack_edges > NEDGE_TOP) ) then
    istat = STAT_FAILURE
    emsg  = "cracked edges' indices must be within [1, 4], update, &
    &fCoh3d8_subelem_module"
    return
  end if

  elem%crack_edges = crack_edges

end subroutine update_fCoh3d8_subelem



pure subroutine extract_fCoh3d8_subelem (elem, fstat, node_connec, edge_connec,&
& ID_matlist, crack_edges, newpartition, local_clock, subelem, subelem_connec)
!** to be cleaned up **
! consider include dm, traction, separation, and remove crack_edges, 
! newpartition, subelem, subelem_connec. extract only the useful components

  type(fCoh3d8_subelem),                    intent(in)  :: elem
  integer,                        optional, intent(out) :: fstat
  integer,           allocatable, optional, intent(out) :: node_connec(:)
  integer,           allocatable, optional, intent(out) :: edge_connec(:)
  integer,                        optional, intent(out) :: ID_matlist
  integer,           allocatable, optional, intent(out) :: crack_edges(:)
  logical,                        optional, intent(out) :: newpartition
  type(program_clock),            optional, intent(out) :: local_clock
  type(sub3d_element),   allocatable, optional, intent(out) :: subelem(:)
  type(int_alloc_array), allocatable, optional, intent(out) :: subelem_connec(:)

  if(present(fstat)) fstat=elem%fstat
  
  if(present(node_connec)) then 
      allocate(node_connec(NNODE))
      node_connec=elem%node_connec
  end if
  
  if(present(edge_connec)) then 
      allocate(edge_connec(NEDGE_TOP))
      edge_connec=elem%edge_connec
  end if
  
  if(present(ID_matlist)) ID_matlist=elem%ID_matlist

  if(present(crack_edges)) then 
      allocate(crack_edges(NEDGE_TOP))
      crack_edges=elem%crack_edges
  end if
  
  if(present(newpartition)) newpartition=elem%newpartition
  
  if(present(local_clock)) local_clock=elem%local_clock
  
  if(present(subelem)) then
    if(allocated(elem%subelem)) then
      allocate(subelem(size(elem%subelem)))
      subelem=elem%subelem
    end if
  end if
  
  if(present(subelem_connec)) then
    if(allocated(elem%subelem_connec)) then
      allocate(subelem_connec(size(elem%subelem_connec)))
      subelem_connec=elem%subelem_connec
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
  integer :: subelfstat
  logical :: last_converged 
  logical :: nofail
  integer :: i, j, l
  
  ! initialize K & F and local variables
  allocate(K_matrix(NDOF,NDOF), F_vector(NDOF))
  K_matrix = ZERO
  F_vector = ZERO
  istat    = STAT_SUCCESS
  emsg     = ''
  subelfstat = 0
  last_converged = .false.
  nofail         = .false.
  i=0; j=0; l=0
  
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
  
  ! - check if last iteration has converged
  if(.not. clock_in_sync(GLOBAL_CLOCK, el%local_clock)) then
      last_converged = .true.
      el%local_clock = GLOBAL_CLOCK
  end if
    

  !**** MAIN CALCULATIONS ****
  
  ! select what to do based on elem's fstat:
  ! if elem is still INTACT: 
  !   - check to see if there's any delamination onset; if so, update fstat to
  !     FSTAT1_FCOH, indicating the onset of delamination prior to matrix crack
  !     in the neighbour top ply
  !   - check edge status variables to see if the neighbour top ply has matrix 
  !     crack; if so, partition element and update fstat to FSTAT2_FCOH, then
  !     integrate and assemble sub elems
  ! if elem is at FSTAT1_FCOH:
  !   - check edge status variables to see if the neighbour top ply has matrix 
  !     crack; if so, partition element and update fstat to FSTAT2_FCOH, then
  !     integrate and assemble sub elems
  ! if elem is at FSTAT2_FCOH:
  !   - no need to check edge status variables; no need to update partition.
  !     just update nodal values and integrate and assemble sub elems

  select case (el%fstat)
  
    case (INTACT)
        ! check if elem has started to delaminate
        call extract (el%subelem(1), fstat=subelfstat)
        ! update fstat if elem reached delamination onset
        if (subelfstat /= INTACT) then 
            el%fstat = FSTAT1_FCOH
        end if
        ! partition and integrate
        call edge_status_partition(el)
        call integrate_assemble_subelem(el, K_matrix, F_vector, istat, emsg, nofail) 
    
    case (FSTAT1_FCOH)
      ! elem already delaminated
        ! partition and integrate
        call edge_status_partition(el)
        call integrate_assemble_subelem(el, K_matrix, F_vector, istat, emsg, nofail) 
    
    case (FSTAT2_FCOH)
      ! elem already partitioned; need to update mnodes, 
      ! then integrate and assemble
        ! update mnodes
        call update_mnode(el)
        ! integrate and assemble
        call integrate_assemble_subelem(el, K_matrix, F_vector, istat, emsg, nofail) 
  
    case default
      ! this place should NOT be reached
        istat = STAT_FAILURE
        emsg  = 'unsupported elfstat value in fCoh3d8 elem module'
    
  end select
  
  !**** END MAIN CALCULATIONS ****

  ! if above loop is exit on error, clean up and exit program
  if (istat == STAT_FAILURE) then
    K_matrix = ZERO
    F_vector = ZERO
    return
  end if
  
  ! update elem sdv for output 
  if(last_converged) call update_sdv(el)
  
  ! update to dummy arg. elem before successful 
  elem = el

end subroutine integrate_fCoh3d8_subelem



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
  real(DP),	allocatable  :: Ki(:,:), Fi(:)
  integer,  allocatable  :: dofcnc(:)
  logical :: nofail
  integer :: i, j, l
    
  istat  = STAT_SUCCESS
  emsg   = ''
  nofail = .false.
  i=0; j=0; l=0
  
  ! empty K and F
  K_matrix = ZERO
  F_vector = ZERO
  
  if (present(nofailure)) nofail = nofailure

  ! loop over all sub elements to
  ! integrate their system arrays and assemble together
  do i = 1, size(elem%subelem)
  
      call integrate(elem%subelem(i), Ki, Fi, istat, emsg, nofail)
      if (istat == STAT_FAILURE) exit
      
      if(allocated(dofcnc)) deallocate(dofcnc)
      allocate(dofcnc(size(Fi))); dofcnc = 0
      
      ! loop over no. of nodes in sub elem i
      do j=1, size(elem%subelem_connec(i)%array) 
        do l=1, ndim
          ! dof indices of the jth node of sub elem i 
          dofcnc((j-1)*ndim+l) = (elem%subelem_connec(i)%array(j)-1)*ndim+l
        end do
      end do
      
      call assembleKF(K_matrix, F_vector, Ki, Fi, dofcnc, istat, emsg)
      if (istat == STAT_FAILURE) exit
      
      deallocate(Ki)
      deallocate(Fi)
      deallocate(dofcnc)
      
  end do
  
  ! clean up and exit if loop exit upon error
  if (istat == STAT_FAILURE) then
    K_matrix = ZERO
    F_vector = ZERO 
    if(allocated(Ki))     deallocate(Ki)
    if(allocated(Fi))     deallocate(Fi)
    if(allocated(dofcnc)) deallocate(dofcnc)
    return
  end if
  
  if(allocated(Ki))     deallocate(Ki)
  if(allocated(Fi))     deallocate(Fi)
  if(allocated(dofcnc)) deallocate(dofcnc)
    
end subroutine integrate_assemble_subelem



pure subroutine edge_status_partition (el, istat, emsg)

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
  ! *** workings of edge_status, n_crackedges, crack_edges ***
  !
  ! e.g.: element edge 1 and 3 are broken, then:
  !
  ! - n_crackedges=2
  ! - edge_status(1)>0; edge_status(2)=0; edge_status(3)>0; edge_status(4)=0
  ! - crack_edges(1)=1; crack_edges(2)=3; crack_edges(3:)=0
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
    
    ! no need to update partition if elem is already in final partition
    if(el%fstat == FSTAT2_FCOH) return

    ! extract edge status variables from glb edge list
    edge_status(:) = global_edge_list(el%edge_connec(:))
    
    ! find the no. of broken edges
    do i = 1, NEDGE_TOP
      if (el%crack_edges(i)>0) n_crackedges=n_crackedges+1
    end do


    !**** MAIN CALCULATIONS ****
    
    ! update el%fstat w.r.t no. of failed edges and edge status
    ! only update el%fstat to the final partition status, FSTAT2_FCOH, when
    ! the top ply elem has reached its final partition status (MATRIX_CRACK_ELEM
    ! or FIBRE_FAILED_ELEM)
    select case (n_crackedges)
    
      case (0)
      ! adj. ply elem remains INTACT, el's fstat can only be: INTACT, FSTAT1
          if(.not.(el%fstat==INTACT .or. el%fstat==FSTAT1_FCOH)) then
            istat = STAT_FAILURE
            emsg = 'inconsistency btw elfstat and n_crackedges in fCoh3d8 elem!'
            return
          end if
        
      case (1)
      ! adj. ply elem is in transition partition; edge status should be TRANSITION_EDGE
          if(edge_status(el%crack_edges(1)) == TRANSITION_EDGE) then
            ! edge marks the refinement end and the trans elem start in the upper
            ! ply element. for this fCoh elem, no partition is done for this
            ! case
            continue
          else
            ! crack_edges is not correct / edge_status not correct
            istat = STAT_FAILURE
            emsg  = 'wrong edge status for n_crackedges=1 in fCoh3d8'
            return
          end if
        
      case (2)
      ! adj. ply elem could be cracked(final), wake, tip, refinement elem
      ! only update el%fstat if adj. ply elem reaches final partition
          ! store local edge indices to edge 1 and edge2
          edge1 = el%crack_edges(1)
          edge2 = el%crack_edges(2)
          ! sort edge1 and edge2 according to their edge status value
          ! make sure edge_status(edge1) <= edge_status(edge2)
          if (edge_status(edge1) > edge_status(edge2)) then
            edge1 = el%crack_edges(2)
            edge2 = el%crack_edges(1)
          end if
          
          if (edge_status(edge1) <= REFINEMENT_EDGE .and. &
          &   edge_status(edge2) <= REFINEMENT_EDGE) then
          ! ply elem status is refinement elem, no partition for this fCoh elem
            continue
          else if (edge_status(edge1) <= CRACK_TIP_EDGE .and. &
          &        edge_status(edge2) == CRACK_TIP_EDGE) then
          ! ply elem status is crack tip elem, no partition for this fCoh elem
            continue
          else if (edge_status(edge1) < COH_CRACK_EDGE .and. &
          &        edge_status(edge2) >= COH_CRACK_EDGE) then
          ! ply elem status is crack wake elem, no partition for this fCoh elem
            continue
          else if (edge_status(edge1) >= COH_CRACK_EDGE .and. &
          &        edge_status(edge2) >= COH_CRACK_EDGE) then
          ! ply elem status is matrix_crack_elem or fibre_fail_elem, this is the
          ! final partition of the ply elem. 
          ! update fCoh fstat
            el%fstat = FSTAT2_FCOH
          else ! unknown combination
            istat = STAT_FAILURE
            emsg  = 'unknown combination of 2 edge status in fCoh3d8!'
            return
          end if
    
      case default
          istat = STAT_FAILURE
          emsg  = 'unsupported n_crackedges value for edge and el stat update &
          &in fCoh3d8 edge stat partition!'
          return 
    
    end select
        
    if (el%fstat == FSTAT2_FCOH) then
    ! update elem partitions into sub elems when it reaches FSTAT2_FCOH status                
      call update_subelem_connec(elem, edge_status, crack_edges, n_crackedges)
    end if
    
    !**** END MAIN CALCULATIONS ****

end subroutine edge_status_partition


   
pure subroutine update_subelem_connec (elem, edge_status, crack_edges, n_crackedges)

! passed-in variables
type(fCoh3d8_subelem),    intent(inout)   :: elem
integer,                intent(in)      :: edge_status(:), crack_edges(:), n_crackedges



! local variables
type(int_alloc_array), allocatable :: subglbcnc(:)  ! glb cnc of sub elements
integer :: i, j, l                                  ! counters
integer :: ibe, ibe1, ibe2                          ! indices of broken edges
integer :: e1,e2,e3,e4                              ! edge indices, used for partitioning element
integer :: nsub, nbulk                              ! no. of sub elements, and no. of bulk partitions
integer :: jnode, jnode1, jnode2                    ! node index variables

real(DP), allocatable   :: x1(:), x2(:), xc(:)      ! coords of broken edge end nodes (x1 and x2) and crack tip (xc)
real(DP)                :: tratio, tratio1, tratio2 ! tratio=|xc-x1|/|x2-x1|
type(xnode),allocatable :: mnode(:)                 ! material nodes of cohesive elem
real(DP), allocatable   :: Tmatrix(:,:)             ! interpolation matrix btw bottom num nodes and mat nodes


!       initialize local variables

    i=0; j=0; l=0
    e1=0; e2=0; e3=0; e4=0
    ibe=0; ibe1=0; ibe2=0
    nsub=0; nbulk=0
    jnode=0; jnode1=0; jnode2=0
    tratio=ZERO; tratio1=ZERO; tratio2=ZERO
    
    

10      select case (n_crackedges)
    case (0) !- no cracked edge, do nothing
        continue
        
        
    case (1) !- one edge cracked, trans partition
        ! find the index of the broken edge
        ibe=crack_edges(1)

        ! ibe1 must be between 1 to 4
        if(ibe1<1 .or. ibe1>4) then
            write(msg_file,*) 'something wrong in fCoh3d8 update subelem_connec case n_crackedges=1'
            call exit_function
        end if

        ! verify its status variable value
        if(edge_status(ibe)/=TRANSITION_EDGE) then
            write(msg_file,*)'transition partition only accepts edge_status=TRANSITION_EDGE!'
            call exit_function
        end if
        
        ! allocate sub element arrays; in this case, 3 coh3d6 sub3d elements
        nsub=3
        if(allocated(elem%subelem)) deallocate(elem%subelem)
        if(allocated(elem%subelem_connec)) deallocate(elem%subelem_connec)
        if(allocated(subglbcnc)) deallocate(subglbcnc)
        allocate(elem%subelem(nsub))
        allocate(elem%subelem_connec(nsub))
        allocate(subglbcnc(nsub)) 
        ! these coh3d6 elems have 7 num nodes and 6 mat nodes
        do j=1, nsub
            allocate(elem%subelem_connec(j)%array(7))
            allocate(subglbcnc(j)%array(7))
            elem%subelem_connec(j)%array=0
            subglbcnc(j)%array=0    
        end do 
        
        ! allocate 6 material nodes for sub elems
        if(allocated(mnode)) deallocate(mnode)
        allocate(mnode(6))
        
        ! allocate Tmatrix for bottom surface (3 mat nodes interpolated by 4 num nodes)
        if(allocated(Tmatrix)) deallocate(Tmatrix)
        allocate(Tmatrix(3,4))
        
        ! find the neighbouring edges of this broken edge, in counter-clockwise direction
        select case(ibe)
            case (1)
                e1=2;e2=3;e3=4
            case (2)
                e1=3;e2=4;e3=1
            case (3)
                e1=4;e2=1;e3=2
            case (4)
                e1=1;e2=2;e3=3
            case default
                write(msg_file,*)'wrong broken edge in tip partition, subelem_connec'
                call exit_function
        end select
        
        ! find the smaller glb node on the broken edge
        if(elem%node_connec(NODES_ON_TOP_EDGES(3,ibe))<elem%node_connec(NODES_ON_TOP_EDGES(4,ibe))) then
            jnode=NODES_ON_TOP_EDGES(3,ibe)
        else
            jnode=NODES_ON_TOP_EDGES(4,ibe)
        end if
        
        ! find the relative position of crack tip on this edge
        call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(1,ibe))),x=x1)
        call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(2,ibe))),x=x2)
        call extract(global_node_list(elem%node_connec(jnode)),x=xc)
        tratio=distance(x1,xc)/distance(x1,x2)
        
        
        
        !*** sub elm 1 connec; 7 numerical nodes, 6 material nodes     
        
        elem%subelem_connec(1)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2
        elem%subelem_connec(1)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
        elem%subelem_connec(1)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
        elem%subelem_connec(1)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
        elem%subelem_connec(1)%array(5)=NODES_ON_TOP_EDGES(1,e1)
        elem%subelem_connec(1)%array(6)=NODES_ON_TOP_EDGES(2,e1)
        elem%subelem_connec(1)%array(7)=jnode
        
        subglbcnc(1)%array(:)=elem%node_connec(elem%subelem_connec(1)%array(:))
        
        ! update the mnode array
        ! first two mat nodes of bottom surf are the same as num nodes
        mnode(1)=global_node_list(subglbcnc(1)%array(1))
        mnode(2)=global_node_list(subglbcnc(1)%array(2))
        mnode(3)=tratio*global_node_list(subglbcnc(1)%array(4))+(one-tratio)*global_node_list(subglbcnc(1)%array(1))
        mnode(4)=global_node_list(subglbcnc(1)%array(5))
        mnode(5)=global_node_list(subglbcnc(1)%array(6))
        mnode(6)=global_node_list(subglbcnc(1)%array(7))
        
        Tmatrix=ZERO
        Tmatrix(1,1)=one
        Tmatrix(2,2)=one  
        Tmatrix(3,1)=one-tratio
        Tmatrix(3,4)=tratio
        
        call set(elem%subelem(1),eltype='coh3d6',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(1)%array &
        & ,Tmatrix=Tmatrix,mnode=mnode)
        
        
        
        
        !*** sub elm 2 connec; 7 numerical nodes; 6 material nodes
        
        elem%subelem_connec(2)%array(1)=NODES_ON_TOP_EDGES(1,e2)-NNDRL/2
        elem%subelem_connec(2)%array(2)=NODES_ON_TOP_EDGES(2,e2)-NNDRL/2
        elem%subelem_connec(2)%array(3)=NODES_ON_TOP_EDGES(1,ibe)-NNDRL/2
        elem%subelem_connec(2)%array(4)=NODES_ON_TOP_EDGES(2,ibe)-NNDRL/2
        elem%subelem_connec(2)%array(5)=NODES_ON_TOP_EDGES(1,e2)
        elem%subelem_connec(2)%array(6)=NODES_ON_TOP_EDGES(2,e2)
        elem%subelem_connec(2)%array(7)=jnode

        subglbcnc(2)%array(:)=elem%node_connec(elem%subelem_connec(2)%array(:))
        
        ! update the mnode array
        ! first two mat nodes of bottom surf are the same as num nodes
        mnode(1)=global_node_list(subglbcnc(2)%array(1))
        mnode(2)=global_node_list(subglbcnc(2)%array(2))
        mnode(3)=tratio*global_node_list(subglbcnc(2)%array(3))+(one-tratio)*global_node_list(subglbcnc(2)%array(4))
        mnode(4)=global_node_list(subglbcnc(2)%array(5))
        mnode(5)=global_node_list(subglbcnc(2)%array(6))
        mnode(6)=global_node_list(subglbcnc(2)%array(7))
        
        Tmatrix=ZERO
        Tmatrix(1,1)=one
        Tmatrix(2,2)=one  ! first two mat nodes of bottom surf are the same as num nodes
        Tmatrix(3,3)=tratio
        Tmatrix(3,4)=one-tratio
        
        call set(elem%subelem(2),eltype='coh3d6',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(2)%array &
        & ,Tmatrix=Tmatrix,mnode=mnode)
        
        
        
        
        !*** sub elm 3 connec; 7 numerical nodes, 6 material nodes
        
        elem%subelem_connec(3)%array(1)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
        elem%subelem_connec(3)%array(2)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
        elem%subelem_connec(3)%array(3)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2
        elem%subelem_connec(3)%array(4)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2 
        elem%subelem_connec(3)%array(5)=NODES_ON_TOP_EDGES(1,e3)
        elem%subelem_connec(3)%array(6)=NODES_ON_TOP_EDGES(2,e3)
        elem%subelem_connec(3)%array(7)=jnode  
        
        subglbcnc(3)%array(:)=elem%node_connec(elem%subelem_connec(3)%array(:))
        
        ! update the mnode array
        ! first two mat nodes of bottom surf are the same as num nodes
        mnode(1)=global_node_list(subglbcnc(3)%array(1))
        mnode(2)=global_node_list(subglbcnc(3)%array(2))
        mnode(3)=tratio*global_node_list(subglbcnc(3)%array(2))+(one-tratio)*global_node_list(subglbcnc(3)%array(3))
        mnode(4)=global_node_list(subglbcnc(3)%array(5))
        mnode(5)=global_node_list(subglbcnc(3)%array(6))
        mnode(6)=global_node_list(subglbcnc(3)%array(7))
        
        Tmatrix=ZERO
        Tmatrix(1,1)=one
        Tmatrix(2,2)=one  ! first two mat nodes of bottom surf are the same as num nodes
        Tmatrix(3,2)=tratio
        Tmatrix(3,3)=one-tratio
        
        call set(elem%subelem(3),eltype='coh3d6',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(3)%array &
        & ,Tmatrix=Tmatrix,mnode=mnode)
     


    case (2) !- two edges cracked
     
        ibe1=min(crack_edges(1),crack_edges(2))   ! local edge index of 1st broken edge
        ibe2=max(crack_edges(1),crack_edges(2)) 
        
        
        ! ibe1 must be between 1 to 3, and ibe2 between 2 to 4, with ibe2 > ibe1
        if(ibe1<1 .or. ibe1>3 .or. ibe2<2 .or. ibe2>4 .or. ibe2<=ibe1) then
            write(msg_file,*) 'something wrong in fCoh3d8 update subelem_connec case n_crackedges=2'
            call exit_function
        end if

        
        ! determine partition based on the indices of the two broken edges
        !   partition: no. of bulk sub domains
        !   e1 - e4: re-index edges to facilitate partitioning domain
        select case(ibe1)
            case(1)
                select case(ibe2)
                    case(2)
                        nbulk=4
                        e1=1; e2=2; e3=3; e4=4
                    case(3)
                        nbulk=2
                        e1=1; e2=2; e3=3; e4=4
                    case(4)
                        nbulk=4
                        e1=4; e2=1; e3=2; e4=3
                    case default
                        write(msg_file,*)'wrong 2nd broken edge in update subelem_connec fCoh3d8'
                        call exit_function
                end select
            case(2)
                select case(ibe2)
                    case(3)
                        nbulk=4
                        e1=2; e2=3; e3=4; e4=1
                    case(4)
                        nbulk=2
                        e1=2; e2=3; e3=4; e4=1
                    case default
                        write(msg_file,*)'wrong 2nd broken edge in update subelem_connec fCoh3d8'
                        call exit_function
                end select
            case(3)
                if(ibe2==4) then
                    nbulk=4
                    e1=3; e2=4; e3=1; e4=2
                else
                    write(msg_file,*)'wrong 2nd broken edge in update subelem_connec fCoh3d8'
                    call exit_function                    
                end if    
            case default
                write(msg_file,*)'wrong broken edge in update subelem_connec fCoh3d8'
                call exit_function
        end select
        
        
        select case(nbulk)
            case(2)
            ! two quad subdomains

                nsub=2  ! only two coh3d8 sub3d elems

                if(allocated(elem%subelem)) deallocate(elem%subelem)
                if(allocated(elem%subelem_connec)) deallocate(elem%subelem_connec)
                if(allocated(subglbcnc)) deallocate(subglbcnc)
                allocate(elem%subelem(nsub))
                allocate(elem%subelem_connec(nsub))
                allocate(subglbcnc(nsub))
                
                ! allocate 8 numerical nodes for sub elems
                do j=1, nsub
                    allocate(elem%subelem_connec(j)%array(8))
                    allocate(subglbcnc(j)%array(8))
                    elem%subelem_connec(j)%array=0
                    subglbcnc(j)%array=0
                end do
                
                ! allocate 8 material nodes for sub elems
                if(allocated(mnode)) deallocate(mnode)
                allocate(mnode(8))
                
                ! allocate Tmatrix for bottom surface (4 mat nodes interpolated by 4 num nodes)
                if(allocated(Tmatrix)) deallocate(Tmatrix)
                allocate(Tmatrix(4,4))
                
                ! find the smaller glb fl. node on the 1st broken edge
                if(elem%node_connec(NODES_ON_TOP_EDGES(3,e1))<elem%node_connec(NODES_ON_TOP_EDGES(4,e1))) then
                    jnode1=NODES_ON_TOP_EDGES(3,e1)
                else
                    jnode1=NODES_ON_TOP_EDGES(4,e1)
                end if
                
                ! find the relative position of crack tip on this edge
                call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(1,e1))),x=x1)
                call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(2,e1))),x=x2)
                call extract(global_node_list(elem%node_connec(jnode1)),x=xc)
                tratio1=distance(x1,xc)/distance(x1,x2)
                
                ! find the smaller glb fl. node on the 2nd broken edge
                if(elem%node_connec(NODES_ON_TOP_EDGES(3,e3))<elem%node_connec(NODES_ON_TOP_EDGES(4,e3))) then
                    jnode2=NODES_ON_TOP_EDGES(3,e3)
                else
                    jnode2=NODES_ON_TOP_EDGES(4,e3)
                end if
                
                ! find the relative position of crack tip on this edge
                call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(1,e3))),x=x1)
                call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(2,e3))),x=x2)
                call extract(global_node_list(elem%node_connec(jnode2)),x=xc)
                tratio2=distance(x1,xc)/distance(x1,x2)
                
                
                !*** sub elm 1 connec
                elem%subelem_connec(1)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2
                elem%subelem_connec(1)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
                elem%subelem_connec(1)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
                elem%subelem_connec(1)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
                elem%subelem_connec(1)%array(5)=NODES_ON_TOP_EDGES(1,e1)
                elem%subelem_connec(1)%array(6)=NODES_ON_TOP_EDGES(3,e1); if(edge_status(e1)<COH_CRACK_EDGE) elem%subelem_connec(1)%array(6)=jnode1
                elem%subelem_connec(1)%array(7)=NODES_ON_TOP_EDGES(4,e3); if(edge_status(e3)<COH_CRACK_EDGE) elem%subelem_connec(1)%array(7)=jnode2
                elem%subelem_connec(1)%array(8)=NODES_ON_TOP_EDGES(2,e3)
                
                subglbcnc(1)%array(:)=elem%node_connec(elem%subelem_connec(1)%array(:))
                
                ! update the mnode array
                mnode(1)=global_node_list(subglbcnc(1)%array(1))
                mnode(2)=tratio1*global_node_list(subglbcnc(1)%array(1))+(one-tratio1)*global_node_list(subglbcnc(1)%array(2))
                mnode(3)=tratio2*global_node_list(subglbcnc(1)%array(3))+(one-tratio2)*global_node_list(subglbcnc(1)%array(4))
                mnode(4)=global_node_list(subglbcnc(1)%array(4))
                mnode(5)=global_node_list(subglbcnc(1)%array(5))
                mnode(6)=global_node_list(subglbcnc(1)%array(6))
                mnode(7)=global_node_list(subglbcnc(1)%array(7))
                mnode(8)=global_node_list(subglbcnc(1)%array(8))
                
                Tmatrix=ZERO
                Tmatrix(1,1)=one
                Tmatrix(2,1)=tratio1
                Tmatrix(2,2)=one-tratio1
                Tmatrix(3,3)=tratio2
                Tmatrix(3,4)=one-tratio2
                Tmatrix(4,4)=one
                
                call set(elem%subelem(1),eltype='coh3d8',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(1)%array &
                & ,Tmatrix=Tmatrix,mnode=mnode)
                
                
                
                
                
                !*** sub elm 2 connec
                elem%subelem_connec(2)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2
                elem%subelem_connec(2)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
                elem%subelem_connec(2)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
                elem%subelem_connec(2)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
                elem%subelem_connec(2)%array(5)=NODES_ON_TOP_EDGES(4,e1); if(edge_status(e1)<COH_CRACK_EDGE) elem%subelem_connec(2)%array(5)=jnode1
                elem%subelem_connec(2)%array(6)=NODES_ON_TOP_EDGES(2,e1)
                elem%subelem_connec(2)%array(7)=NODES_ON_TOP_EDGES(1,e3)
                elem%subelem_connec(2)%array(8)=NODES_ON_TOP_EDGES(3,e3); if(edge_status(e3)<COH_CRACK_EDGE) elem%subelem_connec(2)%array(8)=jnode2

                subglbcnc(2)%array(:)=elem%node_connec(elem%subelem_connec(2)%array(:))
                
                ! update the mnode array
                mnode(1)=tratio1*global_node_list(subglbcnc(2)%array(1))+(one-tratio1)*global_node_list(subglbcnc(2)%array(2))
                mnode(2)=global_node_list(subglbcnc(2)%array(2))
                mnode(3)=global_node_list(subglbcnc(2)%array(3))
                mnode(4)=tratio2*global_node_list(subglbcnc(2)%array(3))+(one-tratio2)*global_node_list(subglbcnc(2)%array(4))
                mnode(5)=global_node_list(subglbcnc(2)%array(5))
                mnode(6)=global_node_list(subglbcnc(2)%array(6))
                mnode(7)=global_node_list(subglbcnc(2)%array(7))
                mnode(8)=global_node_list(subglbcnc(2)%array(8))
                
                Tmatrix=ZERO
                Tmatrix(1,1)=tratio1
                Tmatrix(1,2)=one-tratio1
                Tmatrix(2,2)=one
                Tmatrix(3,3)=one
                Tmatrix(4,3)=tratio2
                Tmatrix(4,4)=one-tratio2
                
                call set(elem%subelem(2),eltype='coh3d8',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(2)%array &
                & ,Tmatrix=Tmatrix,mnode=mnode)
                

                
                
                
                
            case(4)
            ! four triangular subdomains

                nsub=4

                if(allocated(elem%subelem)) deallocate(elem%subelem)
                if(allocated(elem%subelem_connec)) deallocate(elem%subelem_connec)
                if(allocated(subglbcnc)) deallocate(subglbcnc)
                allocate(elem%subelem(nsub))
                allocate(elem%subelem_connec(nsub))
                allocate(subglbcnc(nsub))
                do j=1, nsub   ! coh3d6 elem cnc
                    allocate(elem%subelem_connec(j)%array(7))
                    allocate(subglbcnc(j)%array(7))
                    elem%subelem_connec(j)%array=0
                    subglbcnc(j)%array=0
                end do

                ! allocate 6 material nodes for sub elems
                if(allocated(mnode)) deallocate(mnode)
                allocate(mnode(6))
                
                ! allocate Tmatrix for bottom surface (3 mat nodes interpolated by 4 num nodes)
                if(allocated(Tmatrix)) deallocate(Tmatrix)
                allocate(Tmatrix(3,4))

                
                ! find the smaller glb fl. node on the 1st broken edge
                if(elem%node_connec(NODES_ON_TOP_EDGES(3,e1))<elem%node_connec(NODES_ON_TOP_EDGES(4,e1))) then
                    jnode1=NODES_ON_TOP_EDGES(3,e1)
                else
                    jnode1=NODES_ON_TOP_EDGES(4,e1)
                end if
                
                ! find the relative position of crack tip on this edge
                call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(1,e1))),x=x1)
                call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(2,e1))),x=x2)
                call extract(global_node_list(elem%node_connec(jnode1)),x=xc)
                tratio1=distance(x1,xc)/distance(x1,x2)
                
                ! find the smaller glb fl. node on the 2nd broken edge
                if(elem%node_connec(NODES_ON_TOP_EDGES(3,e2))<elem%node_connec(NODES_ON_TOP_EDGES(4,e2))) then
                    jnode2=NODES_ON_TOP_EDGES(3,e2)
                else
                    jnode2=NODES_ON_TOP_EDGES(4,e2)
                end if
                
                ! find the relative position of crack tip on this edge
                call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(1,e2))),x=x1)
                call extract(global_node_list(elem%node_connec(NODES_ON_TOP_EDGES(2,e2))),x=x2)
                call extract(global_node_list(elem%node_connec(jnode2)),x=xc)
                tratio2=distance(x1,xc)/distance(x1,x2)
                
                
                
                !*** sub elm 1 connec
                elem%subelem_connec(1)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2 ! upper surf nodes
                elem%subelem_connec(1)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
                elem%subelem_connec(1)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
                elem%subelem_connec(1)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
                elem%subelem_connec(1)%array(5)=NODES_ON_TOP_EDGES(4,e1); if(edge_status(e1)<COH_CRACK_EDGE) elem%subelem_connec(1)%array(5)=jnode1
                elem%subelem_connec(1)%array(6)=NODES_ON_TOP_EDGES(1,e2)
                elem%subelem_connec(1)%array(7)=NODES_ON_TOP_EDGES(3,e2); if(edge_status(e2)<COH_CRACK_EDGE) elem%subelem_connec(1)%array(7)=jnode2
                
                subglbcnc(1)%array(:)=elem%node_connec(elem%subelem_connec(1)%array(:))
                
                mnode(1)=tratio1*global_node_list(subglbcnc(1)%array(1))+(one-tratio1)*global_node_list(subglbcnc(1)%array(2))
                mnode(2)=global_node_list(subglbcnc(1)%array(2))
                mnode(3)=tratio2*global_node_list(subglbcnc(1)%array(2))+(one-tratio2)*global_node_list(subglbcnc(1)%array(3))
                mnode(4)=global_node_list(subglbcnc(1)%array(5))
                mnode(5)=global_node_list(subglbcnc(1)%array(6))
                mnode(6)=global_node_list(subglbcnc(1)%array(7))
                
                Tmatrix=ZERO
                Tmatrix(1,1)=tratio1
                Tmatrix(1,2)=one-tratio1
                Tmatrix(2,2)=one
                Tmatrix(3,2)=tratio2
                Tmatrix(3,3)=one-tratio2
                
                call set(elem%subelem(1),eltype='coh3d6',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(1)%array &
                & ,Tmatrix=Tmatrix,mnode=mnode)

                
                !*** sub elm 2 connec
                elem%subelem_connec(2)%array(1)=NODES_ON_TOP_EDGES(1,e2)-NNDRL/2 ! upper surf nodes
                elem%subelem_connec(2)%array(2)=NODES_ON_TOP_EDGES(2,e2)-NNDRL/2
                elem%subelem_connec(2)%array(3)=NODES_ON_TOP_EDGES(1,e4)-NNDRL/2
                elem%subelem_connec(2)%array(4)=NODES_ON_TOP_EDGES(2,e4)-NNDRL/2
                elem%subelem_connec(2)%array(5)=NODES_ON_TOP_EDGES(4,e2); if(edge_status(e2)<COH_CRACK_EDGE) elem%subelem_connec(2)%array(5)=jnode2
                elem%subelem_connec(2)%array(6)=NODES_ON_TOP_EDGES(1,e3)
                elem%subelem_connec(2)%array(7)=NODES_ON_TOP_EDGES(2,e3) 
                
                subglbcnc(2)%array(:)=elem%node_connec(elem%subelem_connec(2)%array(:))
                
                mnode(1)=tratio2*global_node_list(subglbcnc(2)%array(1))+(one-tratio2)*global_node_list(subglbcnc(2)%array(2))
                mnode(2)=global_node_list(subglbcnc(2)%array(2))
                mnode(3)=global_node_list(subglbcnc(2)%array(3))
                mnode(4)=global_node_list(subglbcnc(2)%array(5))
                mnode(5)=global_node_list(subglbcnc(2)%array(6))
                mnode(6)=global_node_list(subglbcnc(2)%array(7))
                
                Tmatrix=ZERO
                Tmatrix(1,1)=tratio2
                Tmatrix(1,2)=one-tratio2
                Tmatrix(2,2)=one
                Tmatrix(3,3)=one
                
                call set(elem%subelem(2),eltype='coh3d6',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(2)%array &
                & ,Tmatrix=Tmatrix,mnode=mnode)
                
                
                !*** sub elm 3 connec
                elem%subelem_connec(3)%array(1)=NODES_ON_TOP_EDGES(1,e4)-NNDRL/2 ! upper surf nodes
                elem%subelem_connec(3)%array(2)=NODES_ON_TOP_EDGES(2,e4)-NNDRL/2
                elem%subelem_connec(3)%array(3)=NODES_ON_TOP_EDGES(1,e2)-NNDRL/2
                elem%subelem_connec(3)%array(4)=NODES_ON_TOP_EDGES(2,e2)-NNDRL/2
                elem%subelem_connec(3)%array(5)=NODES_ON_TOP_EDGES(1,e4)
                elem%subelem_connec(3)%array(6)=NODES_ON_TOP_EDGES(2,e4)
                elem%subelem_connec(3)%array(7)=NODES_ON_TOP_EDGES(3,e1); if(edge_status(e1)<COH_CRACK_EDGE) elem%subelem_connec(3)%array(7)=jnode1
                
                subglbcnc(3)%array(:)=elem%node_connec(elem%subelem_connec(3)%array(:))
                
                mnode(1)=global_node_list(subglbcnc(3)%array(1))
                mnode(2)=global_node_list(subglbcnc(3)%array(2))
                mnode(3)=tratio1*global_node_list(subglbcnc(3)%array(2))+(one-tratio1)*global_node_list(subglbcnc(3)%array(3))
                mnode(4)=global_node_list(subglbcnc(3)%array(5))
                mnode(5)=global_node_list(subglbcnc(3)%array(6))
                mnode(6)=global_node_list(subglbcnc(3)%array(7))
                
                Tmatrix=ZERO
                Tmatrix(1,1)=one
                Tmatrix(2,2)=one
                Tmatrix(3,2)=tratio1
                Tmatrix(3,3)=one-tratio1
                
                call set(elem%subelem(3),eltype='coh3d6',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(3)%array &
                & ,Tmatrix=Tmatrix,mnode=mnode)
                
                
                !*** sub elm 4 connec
                elem%subelem_connec(4)%array(1)=NODES_ON_TOP_EDGES(1,e1)-NNDRL/2 ! upper surf nodes
                elem%subelem_connec(4)%array(2)=NODES_ON_TOP_EDGES(2,e1)-NNDRL/2
                elem%subelem_connec(4)%array(3)=NODES_ON_TOP_EDGES(1,e3)-NNDRL/2
                elem%subelem_connec(4)%array(4)=NODES_ON_TOP_EDGES(2,e3)-NNDRL/2
                elem%subelem_connec(4)%array(5)=NODES_ON_TOP_EDGES(3,e1); if(edge_status(e1)<COH_CRACK_EDGE) elem%subelem_connec(4)%array(5)=jnode1
                elem%subelem_connec(4)%array(6)=NODES_ON_TOP_EDGES(4,e2); if(edge_status(e2)<COH_CRACK_EDGE) elem%subelem_connec(4)%array(6)=jnode2
                elem%subelem_connec(4)%array(7)=NODES_ON_TOP_EDGES(2,e3)             
                
                subglbcnc(4)%array(:)=elem%node_connec(elem%subelem_connec(4)%array(:))
                
                mnode(1)=tratio1*global_node_list(subglbcnc(4)%array(1))+(one-tratio1)*global_node_list(subglbcnc(4)%array(2))
                mnode(2)=tratio2*global_node_list(subglbcnc(4)%array(2))+(one-tratio2)*global_node_list(subglbcnc(4)%array(3))
                mnode(3)=global_node_list(subglbcnc(4)%array(4))
                mnode(4)=global_node_list(subglbcnc(4)%array(5))
                mnode(5)=global_node_list(subglbcnc(4)%array(6))
                mnode(6)=global_node_list(subglbcnc(4)%array(7))
                
                Tmatrix=ZERO
                Tmatrix(1,1)=tratio1
                Tmatrix(1,2)=one-tratio1
                Tmatrix(2,2)=tratio2
                Tmatrix(2,3)=one-tratio2
                Tmatrix(3,4)=one
                
                call set(elem%subelem(4),eltype='coh3d6',ID_matlist=elem%ID_matlist,glbcnc=subglbcnc(4)%array &
                & ,Tmatrix=Tmatrix,mnode=mnode)
                
                
                
                
            case default
                write(msg_file,*)'wrong nbulk in update subelem_connec fCoh3d8'
                call exit_function
        end select

      
!        case(6)
    ! do not update partition
        

!        case(8)


       
    case default
       write(msg_file,*) 'WARNING: fCoh3d8 update subelem_connec case selection default!'
       
    end select
    
    
    ! deallocate local array
    
    if(allocated(subglbcnc)) deallocate(subglbcnc)
    if(allocated(mnode)) deallocate(mnode)
    if(allocated(Tmatrix)) deallocate(Tmatrix)
    if(allocated(x1)) deallocate(x1)
    if(allocated(x2)) deallocate(x2)
    if(allocated(xc)) deallocate(xc)


end subroutine update_subelem_connec



pure subroutine update_mnode (elem)

  ! passed-in variables
  type(fCoh3d8_subelem),    intent(inout)   :: elem
  
  ! local variables
  type(int_alloc_array), allocatable :: subglbcnc(:)  ! glb cnc of sub elements
  integer :: i, j, l                                  ! counters
  integer :: nsub                             		! no. of sub elements


  real(DP)                :: tratio1, tratio2 		! tratio=|xc-x1|/|x2-x1|
  type(xnode),allocatable :: mnode(:)                 ! material nodes of cohesive elem
  real(DP), allocatable   :: Tmatrix(:,:)             ! interpolation matrix btw bottom num nodes and mat nodes

!       initialize local variables

    i=0; j=0; l=0
    nsub=0
    tratio1=ZERO; tratio2=ZERO

  if(.not.allocated(elem%subelem)) then
    write(msg_file,*)'sub elems not allocated in fCoh3d8 update mnode!'
    call exit_function
  end if

    ! currently only support elem partition FSTAT2_FCOH
    select case(elem%fstat)
    
        case(FSTAT2_FCOH)

            nsub=size(elem%subelem)
            
            select case(nsub)
                    case(2)
                    ! two coh3d8 sub elems

                        allocate(subglbcnc(nsub))
                        
                        ! allocate 8 numerical nodes for sub elems
                        do j=1, nsub
                            allocate(subglbcnc(j)%array(8))
                            subglbcnc(j)%array=0
                            subglbcnc(j)%array(:)=elem%node_connec(elem%subelem_connec(j)%array(:))
                        end do
                        
                        !*** sub elem 1 mnode update
                        
                        call extract(elem%subelem(1),Tmatrix=Tmatrix,mnode=mnode)                                       
                        
                        ! extract ratio values from Tmatrix
                        tratio1=Tmatrix(2,1)
                        tratio2=Tmatrix(3,3)
                        
                        ! update the mnode array
                        mnode(1)=global_node_list(subglbcnc(1)%array(1))
                        mnode(2)=tratio1*global_node_list(subglbcnc(1)%array(1))+(one-tratio1)*global_node_list(subglbcnc(1)%array(2))
                        mnode(3)=tratio2*global_node_list(subglbcnc(1)%array(3))+(one-tratio2)*global_node_list(subglbcnc(1)%array(4))
                        mnode(4)=global_node_list(subglbcnc(1)%array(4))
                        mnode(5)=global_node_list(subglbcnc(1)%array(5))
                        mnode(6)=global_node_list(subglbcnc(1)%array(6))
                        mnode(7)=global_node_list(subglbcnc(1)%array(7))
                        mnode(8)=global_node_list(subglbcnc(1)%array(8))                  
                        
                        call update(elem%subelem(1),mnode=mnode)
                                         
                        
                        !*** sub elm 2 mnode update
                        
                        call extract(elem%subelem(2),Tmatrix=Tmatrix,mnode=mnode) 
                        
                        tratio1=Tmatrix(1,1)
                        tratio2=Tmatrix(4,3)
                        
                        ! update the mnode array
                        mnode(1)=tratio1*global_node_list(subglbcnc(2)%array(1))+(one-tratio1)*global_node_list(subglbcnc(2)%array(2))
                        mnode(2)=global_node_list(subglbcnc(2)%array(2))
                        mnode(3)=global_node_list(subglbcnc(2)%array(3))
                        mnode(4)=tratio2*global_node_list(subglbcnc(2)%array(3))+(one-tratio2)*global_node_list(subglbcnc(2)%array(4))
                        mnode(5)=global_node_list(subglbcnc(2)%array(5))
                        mnode(6)=global_node_list(subglbcnc(2)%array(6))
                        mnode(7)=global_node_list(subglbcnc(2)%array(7))
                        mnode(8)=global_node_list(subglbcnc(2)%array(8))
                        
                        call update(elem%subelem(2),mnode=mnode)                
                        
                        
                        
                    case(4)
                    ! four coh3d6 sub elems


                        if(allocated(subglbcnc)) deallocate(subglbcnc)
                        allocate(subglbcnc(nsub))
                        do j=1, nsub   ! coh3d6 elem cnc
                            allocate(subglbcnc(j)%array(7))
                            subglbcnc(j)%array=0
                            subglbcnc(j)%array(:)=elem%node_connec(elem%subelem_connec(j)%array(:))
                        end do

                                    
                        
                        !*** sub elm 1 connec
                        
                        call extract(elem%subelem(1),Tmatrix=Tmatrix,mnode=mnode)
                        
                        tratio1=Tmatrix(1,1)
                        tratio2=Tmatrix(3,2)
                        
                        mnode(1)=tratio1*global_node_list(subglbcnc(1)%array(1))+(one-tratio1)*global_node_list(subglbcnc(1)%array(2))
                        mnode(2)=global_node_list(subglbcnc(1)%array(2))
                        mnode(3)=tratio2*global_node_list(subglbcnc(1)%array(2))+(one-tratio2)*global_node_list(subglbcnc(1)%array(3))
                        mnode(4)=global_node_list(subglbcnc(1)%array(5))
                        mnode(5)=global_node_list(subglbcnc(1)%array(6))
                        mnode(6)=global_node_list(subglbcnc(1)%array(7))
                        
                        call update(elem%subelem(1),mnode=mnode)
                        
     
                        
                        !*** sub elm 2 connec
                        
                        call extract(elem%subelem(2),Tmatrix=Tmatrix,mnode=mnode)
                        
                        tratio2=Tmatrix(1,1)
                        
                        mnode(1)=tratio2*global_node_list(subglbcnc(2)%array(1))+(one-tratio2)*global_node_list(subglbcnc(2)%array(2))
                        mnode(2)=global_node_list(subglbcnc(2)%array(2))
                        mnode(3)=global_node_list(subglbcnc(2)%array(3))
                        mnode(4)=global_node_list(subglbcnc(2)%array(5))
                        mnode(5)=global_node_list(subglbcnc(2)%array(6))
                        mnode(6)=global_node_list(subglbcnc(2)%array(7))        
                        
                        call update(elem%subelem(2),mnode=mnode)

                        
                        
                        
                        !*** sub elm 3 connec
                        
                        call extract(elem%subelem(3),Tmatrix=Tmatrix,mnode=mnode)
                        
                        tratio1=Tmatrix(3,2)
                        
                        mnode(1)=global_node_list(subglbcnc(3)%array(1))
                        mnode(2)=global_node_list(subglbcnc(3)%array(2))
                        mnode(3)=tratio1*global_node_list(subglbcnc(3)%array(2))+(one-tratio1)*global_node_list(subglbcnc(3)%array(3))
                        mnode(4)=global_node_list(subglbcnc(3)%array(5))
                        mnode(5)=global_node_list(subglbcnc(3)%array(6))
                        mnode(6)=global_node_list(subglbcnc(3)%array(7))
          
                        call update(elem%subelem(3),mnode=mnode)
                        
                        
                        
                        !*** sub elm 4 connec
                        
                        call extract(elem%subelem(4),Tmatrix=Tmatrix,mnode=mnode)
                        
                        tratio1=Tmatrix(1,1)
                        tratio2=Tmatrix(2,2)
                        
                        mnode(1)=tratio1*global_node_list(subglbcnc(4)%array(1))+(one-tratio1)*global_node_list(subglbcnc(4)%array(2))
                        mnode(2)=tratio2*global_node_list(subglbcnc(4)%array(2))+(one-tratio2)*global_node_list(subglbcnc(4)%array(3))
                        mnode(3)=global_node_list(subglbcnc(4)%array(4))
                        mnode(4)=global_node_list(subglbcnc(4)%array(5))
                        mnode(5)=global_node_list(subglbcnc(4)%array(6))
                        mnode(6)=global_node_list(subglbcnc(4)%array(7))          
                        
                        call update(elem%subelem(4),mnode=mnode)                   
                        
                        
                    case default
                        write(msg_file,*)'wrong nbulk in update subelem_connec fCoh3d8'
                        call exit_function
            end select

        case default
            write(msg_file,*)'partition status not supported in fCoh3d8 update mnode!'
            call exit_function
            
    end select
        
    ! deallocate local array
    
    if(allocated(subglbcnc)) deallocate(subglbcnc)
    if(allocated(mnode)) deallocate(mnode)
    if(allocated(Tmatrix)) deallocate(Tmatrix)


end subroutine update_mnode



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
        allocate(elem%sdv(1)%i(1))  ! fstat
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
