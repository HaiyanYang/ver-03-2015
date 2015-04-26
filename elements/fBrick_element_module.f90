module fBrick_element_module
!
!  Purpose:
!    define a floating-node brick element object, with breakable top & bottom
!    surfaces.
!
!  topological definition of this element:
!
!  8__22___(E7)___21___7
!  |\                  |\
!  | \23               | \20
!  |  \(E8)            |  \(E6)     Top    edges (anti-clock wise from front):
!  |   \24             |   \19      E5, E6, E7, E8
!  |    \              |    \
!  |     \5___17__(E5)_|_18__\6
!  |______|____________|      |
! 4\   14 | (E3)  13  3\      |
!   \     |             \     |
!  15\    |            12\    |
! (E4)\   |           (E2)\   |     Bottom edges (anti-clock wise from front):
!    16\  |              11\  |     E1, E2, E3, E4
!       \ |                 \ |
!        \|__________________\|
!         1    9  (E1)  10    2
!
!  :::: bottom surface definition ::::
!  real     nodes on bot edges :  1,  2,  3,  4
!  floating nodes on bot edges :  9, 10, 11, 12, 13, 14, 15, 16
!  nodes of E1 ::    <end nodes: 1, 2>  <fl. nodes:  9, 10> 
!  nodes of E2 ::    <end nodes: 2, 3>  <fl. nodes: 11, 12> 
!  nodes of E3 ::    <end nodes: 3, 4>  <fl. nodes: 13, 14> 
!  nodes of E4 ::    <end nodes: 4, 1>  <fl. nodes: 15, 16> 
!
!  :::: top surface definition ::::
!  real     nodes on top edges :  5,  6,  7,  8
!  floating nodes on top edges : 17, 18, 19, 20, 21, 22, 23, 24
!  nodes of E5 ::    <end nodes: 5, 6>  <fl. nodes: 17, 18> 
!  nodes of E6 ::    <end nodes: 6, 7>  <fl. nodes: 19, 20> 
!  nodes of E7 ::    <end nodes: 7, 8>  <fl. nodes: 21, 22> 
!  nodes of E8 ::    <end nodes: 8, 5>  <fl. nodes: 23, 24> 
!
!
!  topological definition of INTACT element, type brick_element:
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

use parameter_module,       only : NDIM, DP, ZERO, INT_ALLOC_ARRAY
use glb_clock_module,       only : program_clock
use basePly_element_module, only : basePly_element
use baseCoh_element_module, only : baseCoh_element


implicit none

private
  
! NODE AND EDGE PARAMETERS OF THIS ELEMENT
integer, parameter :: NNDRL       = 8,              &
                   &  NEDGE       = 8,              &
                   &  NEDGE_SURF  = 4,              &
                   &  NNDFL       = 2 * NEDGE,      &
                   &  NNODE       = NNDRL + NNDFL,  &
                   &  NDOF        = NDIM * NNODE
  

! NODAL CONNEC OF INTACT  ELEM: REAL NODES ONLY
integer, parameter :: INTACT_ELEM_NODES(8) = [1,2,3,4,5,6,7,8]

integer, parameter :: NODES_ON_EDGES(4,NEDGE)=              &
& reshape([1,2, 9,10,  2,3,11,12,  3,4,13,14,  4,1,15,16,   &
&          5,6,17,18,  6,7,19,20,  7,8,21,22,  8,5,23,24], [4,NEDGE])

integer, parameter :: NODES_ON_TOP_EDGES(4,NEDGE_SURF) =    &
& reshape([5,6,17,18,  6,7,19,20,  7,8,21,22,  8,5,23,24], [4,NEDGE_SURF])

integer, parameter :: NODES_ON_BOT_EDGES(4,NEDGE_SURF) =    &
& reshape([1,2,9,10,   2,3,11,12,  3,4,13,14,  4,1,15,16], [4,NEDGE_SURF])

integer, parameter :: ENDNODES_ON_BOT_EDGES(2,NEDGE_SURF) = &
& reshape([1,2,        2,3,        3,4,        4,1      ], [2,NEDGE_SURF])

type, public :: fBrick_element
  private 

  integer :: curr_status        = 0
  real(DP):: ply_angle          = ZERO
  integer :: node_connec(NNODE) = 0 
  integer :: edge_connec(NEDGE) = 0
  integer :: edge_status_lcl(NEDGE_SURF) = 0
  logical :: newpartition       = .false.
  type(program_clock) :: local_clock
  type(brick_element),   allocatable :: intact_elem
  type(basePly_element), allocatable :: subBulks(:)
  type(baseCoh_element), allocatable :: cohCrack
  type(INT_ALLOC_ARRAY), allocatable :: subBulks_nodes(:)
  type(INT_ALLOC_ARRAY), allocatable :: cohCrack_nodes
  
end type fBrick_element

interface empty
    module procedure empty_fBrick_element
end interface

interface set
    module procedure set_fBrick_element
end interface

interface integrate
    module procedure integrate_fBrick_element
end interface

interface extract
    module procedure extract_fBrick_element
end interface




public :: empty, set, integrate, extract



contains



pure subroutine empty_fBrick_element(elem)

  type(fBrick_element), intent(inout) :: elem
  
  type(fBrick_element) :: elem_lcl
  
  elem = elem_lcl

end subroutine empty_fBrick_element



pure subroutine set_fBrick_element(elem, ply_angle, node_connec, edge_connec, &
& istat, emsg)
! Purpose:
! to set the element ready for first use
use parameter_module,      only : MSGLENGTH, STAT_FAILURE, STAT_SUCCESS
use brick_element_module,  only : set

  type(fBrick_element),     intent(inout) :: elem
  real(DP),                 intent(in)    :: ply_angle
  integer,                  intent(in)    :: node_connec(NNODE)
  integer,                  intent(in)    :: edge_connec(NEDGE)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  
  ! local copy of elem
  type(fBrick_element) :: elem_lcl
  ! global_connec of sub element
  integer, allocatable :: global_connec(:)
  ! location for emsg
  character(len=MSGLENGTH) :: msgloc

  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ' set, fBrick_element_module'
  
  ! check validity of inputs
  
  ! check node connec
  if ( any(node_connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'node connec indices must be >=1,'//trim(msgloc)
    return
  end if
  
  ! check node connec
  if ( any(edge_connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'edge connec indices must be >=1,'//trim(msgloc)
    return
  end if
  
  ! update to elem_lcl first
  elem_lcl%ply_angle   = ply_angle
  elem_lcl%node_connec = node_connec
  elem_lcl%edge_connec = edge_connec
  
  ! allocate INTACT elem
  allocate(elem_lcl%intact_elem)
  allocate(global_connec(NNDRL))

  ! populate the global connec of INTACT element
  global_connec(:) = node_connec(INTACT_ELEM_NODES(:))

  ! set the INTACT element
  call set (elem_lcl%intact_elem, connec=global_connec, ply_angle=ply_angle, &
  & istat=istat, emsg=emsg)

  ! if an error is encountered in set, clean up and exit program
  if (istat == STAT_FAILURE) then
    if (allocated(global_connec)) deallocate(global_connec)
    emsg = emsg//trim(msgloc)
    return
  end if

  ! update to dummy arg. elem before successful return
  elem = elem_lcl

  if (allocated(global_connec)) deallocate(global_connec)

end subroutine set_fBrick_element



pure subroutine extract_fBrick_element(elem, curr_status, edge_status_lcl, &
& intact_elem, subBulks, cohCrack)

  type(fBrick_element), intent(in)  :: elem
  integer,    optional, intent(out) :: curr_status
  integer,    optional, intent(out) :: edge_status_lcl(NEDGE_SURF)
  type(brick_element),   allocatable, optional, intent(out) :: intact_elem
  type(basePly_element), allocatable, optional, intent(out) :: subBulks(:)
  type(baseCoh_element), allocatable, optional, intent(out) :: cohCrack
    
  if(present(curr_status)) curr_status=elem%curr_status
  
  if(present(edge_status_lcl)) edge_status_lcl=elem%edge_status_lcl
  
  if(present(intact_elem)) then
      if(allocated(elem%intact_elem)) then
          allocate(intact_elem)
          intact_elem=elem%intact_elem
      end if
  end if
  
  if(present(subBulks)) then
      if(allocated(elem%subBulks)) then
          allocate(subBulks(size(elem%subBulks)))
          subBulks=elem%subBulks
      end if
  end if
  
  if(present(cohCrack)) then        
      if(allocated(elem%cohCrack)) then
          allocate(cohCrack)
          cohCrack=elem%cohCrack
      end if
  end if

end subroutine extract_fBrick_element



pure subroutine integrate_fBrick_element(elem, nodes, edge_status, lam_mat &
& coh_mat, K_matrix, F_vector, istat, emsg)
use parameter_module, only :
use xnode_module,             only : xnode
use lamina_material_module,   only : lamina_material
use cohesive_material_module, only : cohesive_material
use global_clock_module, only : GLOBAL_CLOCK, clock_in_sync

  type(fBrick_element),     intent(inout) :: elem 
  type(xnode),              intent(inout) :: nodes(NNODE)
  integer,                  intent(inout) :: edge_status(NEDGE)
  type(lamina_material),    intent(in)    :: lam_mat
  type(cohesive_material),  intent(in)    :: coh_mat
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  !:::: local variables ::::
  ! local copy of intent inout variables
  type(fBrick_element) :: el
  type(xnode)          :: nds(NNODE)
  integer              :: egstat(NEDGE)
  ! sub elem K and F
  real(DP), allocatable :: Ki(:,:), Fi(:)
  ! logical control variables
  logical :: last_converged
  logical :: nofailure
  ! error msg location
  character(len=MSGLENGTH) :: msgloc
  
  ! initialize intent out and local variables
  allocate(K_matrix(NDOF,NDOF), F_vector(NDOF))
  K_matrix = ZERO
  F_vector = ZERO
  istat    = STAT_SUCCESS
  emsg     = ''
  last_converged = .false.
  nofailure      = .false.
  msgloc   = ' integrate, fBrick_element module'
  
  ! copy intent inout arg. to its local alias
  el     = elem
  nds    = nodes
  egstat = edge_status
  
  ! - check if last iteration has converged
  if(.not. clock_in_sync(GLOBAL_CLOCK, el%local_clock) then
    last_converged  = .true.
    el%local_clock  = GLOBAL_CLOCK
    el%newpartition = .false.
  end if


  !---------------------------------------------------------------------!
  !       update elem partition using edge status variable
  !---------------------------------------------------------------------!

  
  if (el%curr_status < MATRIX_CRACK_ELEM) then 
  ! if elem has not yet reached final partition, 
  ! check elem edge status variables and update 
  ! elem status and sub elem cnc
  
      ! store current status value
      elstat = el%curr_status  
      
      ! by default, nofailure is false, 
      ! i.e., failure criterion will be assessed
      nofailure = .false.
      
      ! partition elem according to edge status values
      call edge_status_partition (el, edge_status, istat, emsg)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call clean_up()
        return
      end if
   
      ! new partition for this elem
      ! no material degradation/failure criterion partition 
      ! for 1st iteration of new partition
      if (elstat /= el%curr_status) then
        elem%newpartition = .true.
        nofailure         = .true.
      end if 
      
      ! if elem matrix is not yet failed after the edge status partition
      ! integrate and check the failure criterion partition
      if (el%curr_status < MATRIX_CRACK_ELEM) then
        
        call integrate_assemble_subelem (el, nodes, lam_mat, coh_mat, &
        & K_matrix, F_vector, istat, emsg, nofailure)
        if (istat == STAT_FAILURE) then
          emsg = emsg//trim(msgloc)
          call clean_up()
          return
        end if
        
        if(nofailure) then
          ! 1st iteration of new partition, no failure criterion partition 
          ! for stabilization purpose
          continue
        else
          ! 2nd and later iteration of new partition
          !***** check failure criterion *****
          ! failure criterion partitions elem into MATRIX_CRACK_ELEM/FIBRE_FAIL_ELEM partition
          ! if sub elem matrix fails/fibre fails
          call failure_criterion_partition (el)
          if(el%curr_status >= MATRIX_CRACK_ELEM) then
          ! elem partition is updated by failure criterion partition, i.e., 
          ! new partiton is true
            el%newpartition = .true.
            nofailure       = .true.
          end if 
        end if
      end if
         
  end if 
 
  if(el%curr_status == MATRIX_CRACK_ELEM) then
  ! element matrix is already failed
  ! element is already partitioned into 2 bulks and 1 coh, 
  ! integrate and assemble subelems
  ! bulk sub elems may undergo fibre damage and failure 
  ! (only after coh sub elem starts failing)

      ! by default, no material degradation for fibre
      ! until cohesive sub elem starts to fail
      nofailure = .true.
      
      ! during the increment of new partition, 
      ! no fibre material degradation allowed
      if(el%newpartition) then
        continue    ! nofailure remains true
      else
      ! after that increment, fibre failure is considered for 
      ! matrix fail partition only after coh sub elem starts failing
        call extract (el%cohCrack, fstat=subelstat)
        if (subelstat > INTACT) nofailure = .false.
      end if
       
      ! integrate sub elems
      call integrate_assemble_subelem (el, nodes, lam_mat, coh_mat, &
      & K_matrix, F_vector, istat, emsg, nofailure)
      
      ! check if sub bulks have reached fibre failure onset;
      ! if so, update curr status to fibre failure status FIBRE_FAIL_ELEM 
      ! (no need to update partition)
      if (.not. nofailure) then
        do i = 1, size(el%subelem)
          call extract (el%subBulks(i), fstat=subelstat)
          if (subelstat >= FIBRE_ONSET) then
            elem%curr_status = FIBRE_FAIL_ELEM
            goto 10
            exit
          end if 
        end do
      end if
  
  end if
  
  if(elem%curr_status==FIBRE_FAIL_ELEM) then
  ! element fibre is already failed, integrate and assemble subelem
      
      ! during the increment of new partition, no fibre failure allowed
      if (elem%newpartition) then
          nofailure=.true.  ! no fibre failure modelling
      else
      ! after that increment, fibre failure is considered for fibre fail partition
          nofailure=.false.
      end if    
      
      call integrate_assemble_subelem (el, nodes, lam_mat, coh_mat, &
      & K_matrix, F_vector, istat, emsg, nofailure)
      
  end if
 

  
  !---------------------------------------------------------------------!
  !               deallocate local arrays 
  !---------------------------------------------------------------------!
10if(allocated(subglbcnc)) deallocate(subglbcnc)



end subroutine integrate_fBrick_element



pure subroutine edge_status_partition(elem, nodes, edge_status, istat, emsg)
use xnode_module, only : xnode, extract

  ! passed-in variables
  type(fBrick_element),     intent(inout) :: elem
  type(xnode),              intent(inout) :: nodes
  integer,                  intent(inout) :: edge_status
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variable
  character(len=MSGLENGTH)  :: msgloc
  integer                   :: elstatus, eledgestatus(NEDGE_SURF)
  type(REAL_ALLOC_ARRAY)    :: coords(NNODE)
  integer                   :: nfailedge
  integer                   :: ifailedge(NEDGE_SURF)
  
  integer :: jbe1, jbe2, jnode
  integer :: cross_stat
  crosspoint1, crosspoint2, quad_coords, ...

  ! --------------------------------------------------------------------!
  ! *** workings of edgstat, nfailedge, ifedg ***
  !
  ! e.g.: element edge 1 and 3 are broken, then:
  !
  ! - nfailedge=2
  ! - edge_status(1)>0; edge_status(2)=0; edge_status(3)>0; edge_status(4)=0
  ! - ifailedge(1)=1; ifailedge(2)=3; ifailedge(3:)=0
  !
  ! --------------------------------------------------------------------!

  ! initialize intent out and local variables
  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ' edge status partition'
  elstatus     = 0
  eledgestatus = 0
  nfailedge = 0
  ifailedge = 0
  jbe1  = 0
  jbe2  = 0
  jnode = 0
  cross_stat = 0
  i=0; j=0; l=0; k=0
  
  ! check input validity
  ! it is assumed that the top surf. edge status
  ! are the same as the bottom surf. edge status
  if ( any(edge_status(1:NEDGE_SURF) /= edge_status(NEDGE_SURF+1:NEDGE)) ) then
    istat = STAT_FAILURE
    emsg  = 'incompatible top and bot surf edge status values'//trim(msgloc)
    call clean_up ()
    return
  end if
    
  ! extract elem status variable
  elstatus     = elem%curr_status
  ! extract elem edge status variables stored in the elem component
  eledgestatus = elem%edge_status
  ! extract nodal coords from passed in nodes array
  do i=1, NNODE
    call extract (nodes(i), x=coords(i)%array)
  end do
  
  
  !**** MAIN CALCULATIONS ****
  !
  ! - if elem is INTACT, 
  ! if no edge is broken, elem remain intact; 
  ! if edges are broken, find the most critical edge jbe1, ignore the others; 
  ! then from jbe1 crack point, find the other edge crossed by the crack line;
  ! store the edge indices of these two edges (jbe1 and jbe2);
  ! partition elem to trans. elem or ref. elem based on status of jbe1 & jbe2
  !
  ! - if elem is TRANSITION ELEM,
  ! 1st broken edge, jbe1, is already found and stored in eledgestatus, 
  ! find the other edge crossed by the crack line starting from jbe1 crack point
  ! store its edge index jbe2
  ! partition elem to trans. elem or ref. elem based on status of jbe1 & jbe2
  ! 
  ! - if elem is of OTHER STATUS,
  ! then the two broken edges (jbe1 & jbe2) are already stored in eledgestatus;
  ! just update their edge status with the passed in edge status array
  ! then partition elem to other status based on status of jbe1 & jbe2
  !
  jbe1jbe2: select case (elstatus)
  !
  case (INTACT, TRANSITION_ELEM) jbe1jbe2
      ! if elem is INTACT, return of no edge fails; otherwise, find jbe1
      ! as the most critical edge
      if (elstatus == INTACT)
        if ( count(edge_status(1:NEDGE_SURF) > INTACT) == 0 ) then
          ! no broken edge, elem remains intact, do nothing
          return
        else            
          ! find the most damaged edge index jbe1 and ignore the rest
          jbe1 = maxloc(edge_status(1:NEDGE_SURF))
        end if
      ! if elem is TRANSITION, find jbe1 as the stored edge in eledgestatus
      else
        ! check if there's only 1 broken edge stored
        if ( count(eledgestatus > INTACT) /= 1 ) then
          istat = STAT_FAILURE
          emsg  = 'unexpected no. of failed edges for case &
          & elstatus = TRANSITION ELEM'//trim(msgloc)
          call clean_up()
          return
        end if
        ! find the broken edge index
        jbe1 = maxloc(eledgestatus)
        ! check if the edge status is transition edge
        if (eledgestatus(jbe1) /= TRANSITION_EDGE) then
          istat = STAT_FAILURE
          emsg  = 'unexpected edge status of failed edge for case &
          & elstatus = TRANSITION ELEM'//trim(msgloc)
          call clean_up()
          return
        end if
      end if
      ! find the other edge, jbe2, crossed by the crack line which passes
      ! the crack point on edge jbe1 
      ! use the subroutine crack_elem_cracktip2d for this purpose
      ! some inputs of the subroutine need to be prepared
      ! index of a fl. node on this edge
      jnode = NODES_ON_BOT_EDGES(3,jbe1)
      ! coords of this fl. node is the existing crack point coords
      crackpoint1(:) = coords(jnode)%array(1:2)
      ! 2D nodal coords of the bottom quad surface (first 4 nodes of the elem)
      do i = 1, NEDGE_SURF  
        quad_coords(:,i) = coords(i)%array(1:2)
      end do
      ! find the other edge crossed by the crack line and the other crack point
      call crack_elem_cracktip2d (cracktip_point = crackpoint1, &
      & cracktip_edge_index = jbe1, nedge = NEDGE_SURF,         &
      & crack_angle = elem%ply_angle, coords = quad_coords,     &
      & nodes_on_edges = ENDNODES_ON_BOT_EDGES,                 &
      & istat = istat, emsg = emsg,                             &
      & edge_crack_point = crackpoint2, crack_edge_index = jbe2)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call clean_up()
        return
      end if
  !
  case (REFINEMENT_ELEM, CRACK_TIP_ELEM, &
  &     CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM) jbe1jbe2  
      ! once partitioned, only check the status of stored failed edges; 
      ! other failed edges are ignored
      do i = 1, NEDGE_SURF
        ! find the failed edges in stored edge status array
        if (eledgestatus(i) /= INTACT) then
          ! update total no. of failed edges
          nfailedge = nfailedge + 1 
          ! update the indices of failed edges in local array
          ifailedge(nfailedge) = i
        end if
      end do
      ! check if nfailedge == 2
      if (nfailedge /= 2) then
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of failed edges for case elstatus = &
        & REFINEMENT/CRACK TIP/CRACK WAKE/MATRIX CRACK ELEM'//trim(msgloc)
        call clean_up()
        return
      end if
      ! assign jbe1 and jbe2
      jbe1 = ifailedge(1)
      jbe2 = ifailedge(2)
  !
  case default jbe1jbe2
      istat = STAT_FAILURE
      emsg  = 'unexpected elstatus value'//trim(msgloc)
      call clean_up
      return
  end select jbe1jbe2
  
  ! extract the edge status of jbe1 and jbe2 
  ! from passed in edge_status array to local eledgestatus
  eledgestatus(jbe1) = edge_status(jbe1)
  eledgestatus(jbe2) = edge_status(jbe2)
  
  ! update edge status and elem status, based on jbe1 and jbe2 status
  call update_edge_elem_status (eledgestatus, elstatus, jbe1, jbe2)
  
  !**** END MAIN CALCULATIONS ****
 
  !:::::::::::::::::::::::::::::::::!
  ! update the intent inout arguments
  ! when elstatus is updated
  !:::::::::::::::::::::::::::::::::!
 
  if (elstatus /= elem%curr_status) then
 
    ! update edge jbe2 fl. nodes when elem was previously intact
    ! or transition elem (on both surfs)
    if (elem%curr_status == INTACT .or. &
    &   elem%curr_status == TRANSITION_ELEM) then
      ! update two fl. nodes on edge jbe2 of both top and bot surfaces
      jnode = NODES_ON_BOT_EDGES(3,jbe2)
      call update(nodes(jnode), x(1:2)=crackpoint2(:))
      jnode = NODES_ON_BOT_EDGES(4,jbe2)
      call update(nodes(jnode), x(1:2)=crackpoint2(:))
      jnode = NODES_ON_TOP_EDGES(3,jbe2)
      call update(nodes(jnode), x(1:2)=crackpoint2(:))
      jnode = NODES_ON_TOP_EDGES(4,jbe2)
      call update(nodes(jnode), x(1:2)=crackpoint2(:))
    end if
    
    ! update to passed in edge_status array (both surfs)
    edge_status(jbe1) = eledgestatus(jbe1)
    edge_status(jbe2) = eledgestatus(jbe2)
    edge_status(jbe1+NEDGE_SURF) = eledgestatus(jbe1)
    edge_status(jbe2+NEDGE_SURF) = eledgestatus(jbe2)
    
    ! update to elem components
    elem%curr_status = elstatus
    elem%edge_status = eledgestatus
  
  end if
  
  ! exit program
  call clean_up()
  return


  contains


  pure subroutine update_edge_elem_status (eledgestatus, elstatus, ibe1, ibe2)
  ! update the edge status of the two broken edges and elem status
    integer, intent(inout) :: elstatus, eledgestatus(NEDGE_SURF)
    integer, intent(in)    :: ibe1, ibe2
    
    integer :: jbe1, jbe2
    jbe1 = 0
    jbe2 = 0
    
    ! store the more critical edge of the two in jbe1
    if (eledgestatus(ibe2) > eledgestatus(ibe1) then
      jbe1 = ibe2
      jbe2 = ibe1
    else
      jbe1 = ibe1
      jbe2 = ibe2
    end if
  
    select case (eledgestatus(jbe1))
    ! if 1st broken edge is transition edge, the 2nd can only be transition
    ! edge or intact edge
    case (TRANSITION_EDGE)
        if (eledgestatus(jbe2) == TRANSITION_EDGE) then
          eledgestatus(jbe1) = REFINEMENT_EDGE
          eledgestatus(jbe2) = REFINEMENT_EDGE
          elstatus = REFINEMENT_ELEM
        else if (eledgestatus(jbe2) == INTACT) then
          elstatus = TRANSITION_ELEM
        end if
    ! if 1st broken edge is refinement edge, the 2nd can only be refinement,
    ! transition or intact edge
    case (REFINEMENT_EDGE)
        if (eledgestatus(jbe2) == REFINEMENT_EDGE .or. &
        &   eledgestatus(jbe2) == TRANSITION_EDGE) then
          elstatus = REFINEMENT_ELEM
        else if (eledgestatus(jbe2) == INTACT) then
          eledgestatus(jbe2) = TRANSITION_EDGE
          elstatus = REFINEMENT_ELEM
        end if
    ! if 1st broken edge is crack tip edge, the 2nd can only be crack tip, 
    ! refinement, transition or intact edge
    case (CRACK_TIP_EDGE)
        if (eledgestatus(jbe2) == CRACK_TIP_EDGE) then
          eledgestatus(jbe1) = COH_CRACK_EDGE
          eledgestatus(jbe2) = COH_CRACK_EDGE
          elstatus = MATRIX_CRACK_ELEM
        else if (eledgestatus(jbe2) == REFINEMENT_EDGE) then
          elstatus = CRACK_TIP_ELEM
        else if (eledgestatus(jbe2) == TRANSITION_EDGE .or. &
        &        eledgestatus(jbe2) == INTACT) then
          eledgestatus(jbe2) = REFINEMENT_EDGE
          elstatus = CRACK_TIP_ELEM
        end if
    ! if 1st broken edge is coh crack edge, the 2nd can only be coh crack,
    ! crack tip, refinement, transition or intact edge
    case (COH_CRACK_EDGE)
        if (eledgestatus(jbe2) == COH_CRACK_EDGE) then
          elstatus = MATRIX_CRACK_ELEM
        else if (eledgestatus(jbe2) == CRACK_TIP_EDGE) then
          elstatus = CRACK_WAKE_ELEM
        else if (eledgestatus(jbe2) == REFINEMENT_EDGE .or. &
        &        eledgestatus(jbe2) == TRANSITION_EDGE .or. &
        &        eledgestatus(jbe2) == INTACT) then
          eledgestatus(jbe2) = CRACK_TIP_EDGE
          elstatus = CRACK_WAKE_ELEM
        end if
    end select

  end subroutine update_edge_elem_status


  pure subroutine clean_up ()
  
  end subroutine clean_up

end subroutine edge_status_partition



pure subroutine failure_criterion_partition(elem)
! this subroutine update elem status & partition to MATRIX_CRACK_ELEM if any sub elem is nolonger INTACT

! passed in variables
type(fBrick_element) :: elem

! extracted variables, from glb libraries
integer :: edgstat(NEDGE)               ! status variable array of element edges
type(real_alloc_array) :: coord(NNODE)  ! nodal coord arrays to store the coords of elem nodes extracted from glb node lib

real(DP)    :: xelm(NDIM,NNODE)         ! local copy of element nodal coords

real(DP)    :: theta



! local variable

integer :: nfailedge        ! no. of failed edges in the element
integer :: ifedg(NEDGE)     ! index of failed edges in the element
integer :: elstat           ! local copy of elem curr status
integer :: jbe1,jbe2, jnode ! indices of 2 broken edges, and a variable to hold a glb node index
integer :: iscross          ! indicator of line intersection; >0 if two lines intersect
integer :: i, j, l, k       ! counters
integer :: subelstat        ! sub elem status
logical :: failed           ! true if elem is failed (one of the sub elems reached failure onset)
logical :: fbfail           ! true if elem has fibre failure (INTACT elem reached fibre failure onset)

  
real(DP) :: xp1, yp1, xp2, yp2      ! (x,y) of point 1 and point 2 on the crack line
real(DP) :: x1, y1, z1, x2, y2, z2  ! (x,y) of node 1 and node 2 of an element edge
real(DP) :: xct, yct, zct           ! (x,y) of a crack tip on an edge, i.e., intersection of crack line & edge
real(DP) :: xo, yo                  ! (x,y) of element centroid
real(DP) :: detlc,a1,b1,a2,b2,xmid,ymid
real(DP) :: xmid1,ymid1,zmid1,xmid2,ymid2,zmid2
real(DP) :: xct1,yct1,zct1,xct2,yct2,zct2


! --------------------------------------------------------------------!
!       *** workings of edgstat, nfailedge, ifedg ***
!
!       e.g.: element edge 1 and 3 are broken, then:
!
!           - nfailedge=2
!           - edgstat(1)>0; edgstat(2)=0; edgstat(3)>0; edgstat(4)=0
!           - ifedg(1)=1; ifedg(2)=3; ifedg(3:)=0
!
! --------------------------------------------------------------------!

    ! subroutine only partitions elstat < MATRIX_CRACK_ELEM elems
    if(elem%curr_status>=MATRIX_CRACK_ELEM) return ! elem already failed, no need to proceed

    ! initialize local variables
    
    edgstat=0; xelm=zero; theta=zero
      
    nfailedge=0; ifedg=0; elstat=0
    jbe1=0; jbe2=0; jnode=0
    iscross=0
    i=0; j=0; l=0; k=0
    subelstat=0
    
    failed=.false.
       
    xp1=zero; yp1=zero
    xp2=zero; yp2=zero

    x1=zero; y1=zero; z1=zero
    x2=zero; y2=zero; z2=zero
    
    xct=zero; yct=zero; zct=zero
    
    xo=zero; yo=zero; detlc=zero

    a1=zero; b1=zero; a2=zero; b2=zero; xmid=zero; ymid=zero
    xmid1=zero; ymid1=zero; zmid1=zero; xmid2=zero; ymid2=zero; zmid2=zero
    xct1=zero; yct1=zero; zct1=zero; xct2=zero; yct2=zero; zct2=zero
    
    


!-----------------------------------------------------------------------!
!               EXTRACTION INTERFACE
!           extract variables from global libraries
!-----------------------------------------------------------------------!
    ! extract edge status variables from glb edge library
    edgstat(:)=lib_edge(elem%edge_connec(:))
    
    ! extract nodal coords from glb node library
    do i=1, NNODE
        call extract(lib_node(elem%node_connec(i)),x=coord(i)%array)
        if(allocated(coord(i)%array)) xelm(:,i)=coord(i)%array(:)
    end do
    

    ! extract material orientation (fibre angle)
    theta=elem%ply_angle
    
    elstat=elem%curr_status

    ifedg=elem%ifailedge


!-----------------------------------------------------------------------!
!           check failure criterion on all sub elements
!   elem is judged failed a.l.a. 1 sub elem has reached failure onset
!   elem is judged fibre failed if elstat=INTACT and subelstat=fibre_onset
!-----------------------------------------------------------------------!

    do i=1, size(elem%subelem)
        call extract(elem%subelem(i),curr_status=subelstat)
        if(subelstat>INTACT) failed=.true.
        if(failed) exit
    end do

    if (elstat==INTACT .and. subelstat>=fibre_onset) then
        fbfail=.true.
    end if
!-----------------------------------------------------------------------!

   
        
    if(failed) then
    ! element failed, find newly broken edge and update edge status vars
        
    
      ! find xp1, yp1
        if(elstat .eq. INTACT) then ! element originally INTACT
            ! find centroid
            xo=quarter*(xelm(1,1)+xelm(1,2)+xelm(1,3)+xelm(1,4))
            yo=quarter*(xelm(2,1)+xelm(2,2)+xelm(2,3)+xelm(2,4))
            xp1=xo
            yp1=yo
            ! find xp2, yp2
            xp2=xp1+cos(theta/halfcirc*pi)
            yp2=yp1+sin(theta/halfcirc*pi)
            do i=1,NEDGE/2 
                ! tip corods of edge i
                x1=xelm(1,topo(1,i))
                y1=xelm(2,topo(1,i))
                z1=xelm(3,topo(1,i))
                x2=xelm(1,topo(2,i))
                y2=xelm(2,topo(2,i))
                z2=xelm(3,topo(2,i))
                iscross=0
                xct=zero
                yct=zero
                call klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct)
                if (iscross.gt.0) then
                    edgstat(i)=cohcrack  ! edge i cracked
                    nfailedge=nfailedge+1
                    ifedg(nfailedge)=i   ! store failed edge indices
                    zct=half*(z1+z2)     ! z1 should be == to z2 (shell bottom plane)
                    xelm(:,topo(3,i))=[xct,yct,zct] ! store c tip coords on edge i fl. nd 1
                    xelm(:,topo(4,i))=[xct,yct,zct] ! store c tip coords on edge i fl. nd 2
                 endif
                 if(nfailedge==2) exit ! found 2 broken edges already
            end do


            ! badly shaped element may have large angles; e.g.: 3 or all edges almost parallel to crack, then numerical error may
            ! prevent the algorithm from finding any broken edge or only one broken edge

            if(nfailedge==0) then
            ! use trial lines: connecting midpoints of two edges and find the most parrallel-to-crack one

                ! the following algorithm will find two broken edges
                nfailedge=2
                ! line equation constants of the crack
                a2=sin(theta/halfcirc*pi)
                b2=-cos(theta/halfcirc*pi)
                ! connecting midpoints of two edges and find the most parrallel one
                do i=1,NEDGE/2-1
                    ! tip coords of edge i
                    x1=xelm(1,topo(1,i))
                    y1=xelm(2,topo(1,i))
                    z1=xelm(3,topo(1,i))
                    x2=xelm(1,topo(2,i))
                    y2=xelm(2,topo(2,i))
                    z2=xelm(3,topo(2,i))
                    ! mid point of edge i
                    xmid1=half*(x1+x2)
                    ymid1=half*(y1+y2)
                    zmid1=half*(z1+z2)
                    ! loop over midpoints of other edges
                    do j=i+1,NEDGE/2
                        ! tip coords of edge j
                        x1=xelm(1,topo(1,j))
                        y1=xelm(2,topo(1,j))
                        z1=xelm(3,topo(1,j))
                        x2=xelm(1,topo(2,j))
                        y2=xelm(2,topo(2,j))
                        z2=xelm(3,topo(2,j))
                        ! mid point of edge j
                        xmid2=half*(x1+x2)
                        ymid2=half*(y1+y2)
                        zmid2=half*(z1+z2)
                        ! line equation constants of midpoint1-midpoint2
                        a1=ymid1-ymid2
                        b1=xmid2-xmid1
                        ! initialize detlc and intersection info
                        if(i==1 .and. j==2) then
                            detlc=a1*b2-a2*b1                               
                            ifedg(1)=i   ! store failed edge indices
                            ifedg(2)=j
                            xct1=xmid1
                            yct1=ymid1
                            zct1=zmid1
                            xct2=xmid2
                            yct2=ymid2
                            zct2=zmid2
                        end if
                        ! find the most parallel trial line and update stored info
                        if(abs(a1*b2-a2*b1)<abs(detlc)) then
                            ifedg(1)=i   ! store failed edge indices
                            ifedg(2)=j
                            xct1=xmid1
                            yct1=ymid1
                            zct1=zmid1
                            xct2=xmid2
                            yct2=ymid2
                            zct2=zmid2
                        endif
                    end do
                end do
                edgstat(ifedg(1:2))=cohcrack
                xelm(:,topo(3,ifedg(1)))=[xct1,yct1,zct1] ! store c tip coords on edge i fl. nd 1
                xelm(:,topo(4,ifedg(1)))=[xct1,yct1,zct1] ! store c tip coords on edge i fl. nd 2 
                xelm(:,topo(3,ifedg(2)))=[xct2,yct2,zct2] ! store c tip coords on edge i fl. nd 1
                xelm(:,topo(4,ifedg(2)))=[xct2,yct2,zct2] ! store c tip coords on edge i fl. nd 2 

            else if(nfailedge==1) then                
            ! use trial lines: connecting the existing crack tip (xp1,yp1) to the midpoints of the other 3 edges

                ! the following algorithm will find the second broken edge
                nfailedge=2
                ! existing crack tip coords
                xp1=xct
                yp1=yct
                ! crack line equation constants
                a2=sin(theta/halfcirc*pi)
                b2=-cos(theta/halfcirc*pi)
                ! find the midpoint which forms the most parrallel-to-crack line with (xp1,yp1)
                do i=1,NEDGE/2 
                    if (i==ifedg(1)) cycle
                    ! tip coords of edge i
                    x1=xelm(1,topo(1,i))
                    y1=xelm(2,topo(1,i))
                    z1=xelm(3,topo(1,i))
                    x2=xelm(1,topo(2,i))
                    y2=xelm(2,topo(2,i))
                    z2=xelm(3,topo(2,i))
                    ! mid point of edge i
                    xmid=half*(x1+x2)
                    ymid=half*(y1+y2)
                    ! line equation constants of midpoint-(xp1,yp1)
                    a1=ymid-yp1
                    b1=xp1-xmid
                    ! initialize detlc and intersection info
                    if(i==1 .or. (ifedg(1)==1 .and. i==2)) then
                        detlc=a1*b2-a2*b1
                        ifedg(2)=i   ! store failed edge indices
                        xct=xmid
                        yct=ymid
                        zct=half*(z1+z2)     ! z1 should be == to z2 (shell bottom plane)
                    end if
                    ! find the most parallel trial line and update stored info
                    if(abs(a1*b2-a2*b1)<abs(detlc)) then
                        ifedg(2)=i   ! store failed edge indices
                        xct=xmid
                        yct=ymid
                        zct=half*(z1+z2)     ! z1 should be == to z2 (shell bottom plane)
                    endif
                end do
                edgstat(ifedg(2))=cohcrack
                xelm(:,topo(3,ifedg(2)))=[xct,yct,zct] ! store c tip coords on edge i fl. nd 1
                xelm(:,topo(4,ifedg(2)))=[xct,yct,zct] ! store c tip coords on edge i fl. nd 2 
            end if                

            !***** project edge status variables to top edges *****
            
            nfailedge=4
            ifedg(3)=ifedg(1)+NEDGE/2
            ifedg(4)=ifedg(2)+NEDGE/2
            edgstat(ifedg(3))=edgstat(ifedg(1))
            edgstat(ifedg(4))=edgstat(ifedg(2))
            xelm(1:2,topo(3,ifedg(3)))=xelm(1:2,topo(3,ifedg(1)))
            xelm(1:2,topo(4,ifedg(3)))=xelm(1:2,topo(4,ifedg(1)))
            xelm(1:2,topo(3,ifedg(4)))=xelm(1:2,topo(3,ifedg(2)))
            xelm(1:2,topo(4,ifedg(4)))=xelm(1:2,topo(4,ifedg(2)))
            
            xelm(3,topo(3,ifedg(3)))=half*(xelm(3,topo(1,ifedg(3)))+xelm(3,topo(2,ifedg(3))))
            xelm(3,topo(4,ifedg(3)))=half*(xelm(3,topo(1,ifedg(3)))+xelm(3,topo(2,ifedg(3))))
            xelm(3,topo(3,ifedg(4)))=half*(xelm(3,topo(1,ifedg(4)))+xelm(3,topo(2,ifedg(4))))
            xelm(3,topo(4,ifedg(4)))=half*(xelm(3,topo(1,ifedg(4)))+xelm(3,topo(2,ifedg(4))))
          
        else ! element already partitioned

            ! find (xp1,yp1), (xp2,yp2)
            if (elstat .eq. eltrans) then ! transition elm, already has an edge partitioned

                nfailedge=2
                edgstat(ifedg(1:2))=cohcrack ! cohesive crack
                
                j=ifedg(1)

                ! find xp1, yp1
                xp1=xelm(1,topo(3,j))           
                yp1=xelm(2,topo(3,j))
             
                ! find xp2, yp2
                xp2=xp1+cos(theta/halfcirc*pi)
                yp2=yp1+sin(theta/halfcirc*pi)
                 
                 ! next, find the edge crossed by the crack line
                do i=1,NEDGE/2
                    if (i.eq.j) cycle ! the already cracked edge,go to next edge
                    ! tip corods of edge i
                    x1=xelm(1,topo(1,i))
                    y1=xelm(2,topo(1,i))
                    z1=xelm(3,topo(1,i))
                    x2=xelm(1,topo(2,i))
                    y2=xelm(2,topo(2,i))
                    z2=xelm(3,topo(2,i))
                    iscross=0
                    xct=zero
                    yct=zero
                    call klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct)
                    if (iscross .gt. 0) then
                        edgstat(i)=cohcrack  ! edge i cracked
                        nfailedge=nfailedge+1
                        ifedg(nfailedge)=i   ! 2nd broken edge index is i
                        zct=half*(z1+z2)     ! z1 should be == to z2 (shell bottom plane)
                        xelm(:,topo(3,i))=[xct,yct,zct] ! store c tip coords on edge i fl. nd 1
                        xelm(:,topo(4,i))=[xct,yct,zct] ! store c tip coords on edge i fl. nd 2
                    end if
                    if(nfailedge .eq. 3) exit ! found 2 broken edges already
                end do

                ! badly shaped element may have large angles; e.g.: 3 or all edges almost parallel to crack, then numerical error may
                ! prevent the algorithm from finding any broken edge or only one broken edge

                if(nfailedge == 2) then
                ! use trial lines: connecting the existing crack tip (xp1,yp1) to the midpoints of the other 3 edges
                    ! crack line equation constants
                    a2=sin(theta/halfcirc*pi)
                    b2=-cos(theta/halfcirc*pi)
                    ! find the midpoint which forms the most parrallel-to-crack line with (xp1,yp1)
                    do i=1,NEDGE/2 
                        if (i==ifedg(1)) cycle
                        ! tip coords of edge i
                        x1=xelm(1,topo(1,i))
                        y1=xelm(2,topo(1,i))
                        z1=xelm(3,topo(1,i))
                        x2=xelm(1,topo(2,i))
                        y2=xelm(2,topo(2,i))
                        z2=xelm(3,topo(2,i))
                        ! mid point of edge i
                        xmid=half*(x1+x2)
                        ymid=half*(y1+y2)
                        ! line equation constants of midpoint-(xp1,yp1)
                        a1=ymid-yp1
                        b1=xp1-xmid
                        ! initialize detlc and intersection info
                        if(i==1 .or. (ifedg(1)==1 .and. i==2)) then
                            detlc=a1*b2-a2*b1
                            ifedg(3)=i   ! store failed edge indices
                            xct=xmid
                            yct=ymid
                            zct=half*(z1+z2)     ! z1 should be == to z2 (shell bottom plane)
                        end if
                        ! find the most parallel trial line and update stored info
                        if(abs(a1*b2-a2*b1)<abs(detlc)) then
                            ifedg(3)=i   ! store failed edge indices
                            xct=xmid
                            yct=ymid
                            zct=half*(z1+z2)     ! z1 should be == to z2 (shell bottom plane)
                        end if
                    end do
                    edgstat(ifedg(3))=cohcrack
                    xelm(:,topo(3,ifedg(3)))=[xct,yct,zct] ! store c tip coords on edge i fl. nd 1
                    xelm(:,topo(4,ifedg(3)))=[xct,yct,zct] ! store c tip coords on edge i fl. nd 2 
                end if

                
                !***** project edge status variables to top edges *****
                
                nfailedge=4
                ifedg(4)=ifedg(3)+NEDGE/2
                edgstat(ifedg(4))=cohcrack
                xelm(1:2,topo(3,ifedg(4)))=xelm(1:2,topo(3,ifedg(3)))
                xelm(1:2,topo(4,ifedg(4)))=xelm(1:2,topo(4,ifedg(3)))
                xelm(3,topo(3,ifedg(4)))=half*(xelm(3,topo(1,ifedg(4)))+xelm(3,topo(2,ifedg(4))))
                xelm(3,topo(4,ifedg(4)))=half*(xelm(3,topo(1,ifedg(4)))+xelm(3,topo(2,ifedg(4))))
                
            else if (elstat.gt.eltrans .and. elstat.lt.MATRIX_CRACK_ELEM) then 
            ! element already has two edges partitioned; just update edge status var to coh crack status

                nfailedge=4
                ! update edge status var
                edgstat(ifedg(1:4))=cohcrack ! cohesive crack 
       
            else
                write(msg_file,*) 'unsupported elstat value for failure!'
                call exit_function          
            end if
          
        end if
      
       
      
        ! update the elstat to failed status value
        if(fbfail) then
        ! if INTACT elem reaches fibre failure onset, new elem partition is named FIBRE_FAIL_ELEM
            elstat=FIBRE_FAIL_ELEM
        else
            elstat=MATRIX_CRACK_ELEM
        end if
       
!       update element curr_status and sub-element cnc matrices
        
        elem%curr_status=elstat 
        elem%ifailedge=ifedg 

        call update_subcnc(elem,edgstat,ifedg,nfailedge)


!-----------------------------------------------------------------------!
!                   UPDATE INTERFACE
!               update global libraries
!-----------------------------------------------------------------------!

        ! update glb edge array, only the broken edges' status variables
        lib_edge(elem%edge_connec(ifedg(1:nfailedge)))=edgstat(ifedg(1:nfailedge))
        
        ! update glb node array, only the broken edges' fl. node coord
        do i=1, nfailedge
            j=ifedg(i)          ! local index of broken edge
            
            l=topo(3,j)          ! elem lcl index of broken edge fl. node 1
            coord(l)%array(:)=xelm(:,l)
            k=elem%node_connec(l)   ! global index of broken edge fl. node 1
            call update(lib_node(k),x=coord(l)%array)
            
            l=topo(4,j)          ! elem lcl index of broken edge fl. node 2
            coord(l)%array(:)=xelm(:,l)
            k=elem%node_connec(l)   ! global index of broken edge fl. node 2
            call update(lib_node(k),x=coord(l)%array)
        end do



    end if

    ! deallocate local arrays
    
    
  return  
  end subroutine failure_criterion_partition


   
pure subroutine update_subcnc(elem,edgstat,ifedg,nfailedge)

! passed-in variables
type(fBrick_element),    intent(inout)   :: elem
integer,                intent(in)      :: edgstat(:), ifedg(:), nfailedge



! local variables
type(int_alloc_array), allocatable :: subglbcnc(:)  ! glb cnc of sub elements
integer :: i, j, l                                  ! counters
integer :: ibe, ibe1, ibe2                          ! indices of broken edges
integer :: e1,e2,e3,e4                              ! edge indices, used for partitioning element
integer :: nsub, nbulk                              ! no. of sub elements, and no. of bulk partitions
integer :: jnode, jnode1, jnode2, jnode3, jnode4    ! node index variables
logical :: iscoh



!       initialize local variables

    i=0; j=0; l=0
    e1=0; e2=0; e3=0; e4=0
    ibe=0; ibe1=0; ibe2=0
    nsub=0; nbulk=0
    jnode=0; jnode1=0; jnode2=0; jnode3=0; jnode4=0
    iscoh=.false.
    


10      select case (nfailedge)
    case (0) !- no cracked edge, do nothing
        continue
        
        
    case (2) !- one edge cracked, trans partition
        ! find the index of the broken edge
        ibe=ifedg(1)

        ! verify its status variable value
        if(edgstat(ibe)/=TRANSITION_EDGE) then
            write(msg_file,*)'transition partition only accepts edgstat=TRANSITION_EDGE!'
            call exit_function
        end if
        
        ! allocate sub element arrays; in this case, 3 wedge elements
        nsub=3
        if(allocated(elem%subelem)) deallocate(elem%subelem)
        if(allocated(elem%subcnc)) deallocate(elem%subcnc)
        if(allocated(subglbcnc)) deallocate(subglbcnc)
        allocate(elem%subelem(nsub))
        allocate(elem%subcnc(nsub))
        allocate(subglbcnc(nsub))
        do j=1, nsub
            allocate(elem%subcnc(j)%array(6))
            allocate(subglbcnc(j)%array(6))
            elem%subcnc(j)%array=0
            subglbcnc(j)%array=0
        end do             
        
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
                write(msg_file,*)'wrong broken edge in tip partition, subcnc'
                call exit_function
        end select
        
        ! find the smaller glb node on the broken edge
        if(elem%node_connec(topo(3,ibe))<elem%node_connec(topo(4,ibe))) then
            jnode=topo(3,ibe)
        else
            jnode=topo(4,ibe)
        end if
        
        ! sub elm 1 connec
        elem%subcnc(1)%array(1)=topo(1,e1)
        elem%subcnc(1)%array(2)=topo(2,e1)
        elem%subcnc(1)%array(3)=jnode
        
        elem%subcnc(1)%array(4:5)=elem%subcnc(1)%array(1:2)+NNDRL/2 ! upper surf nodes
        elem%subcnc(1)%array(6)=elem%subcnc(1)%array(3)+NNDFL/2
        
        subglbcnc(1)%array(:)=elem%node_connec(elem%subcnc(1)%array(:))
        
        ! sub elm 2 connec
        elem%subcnc(2)%array(1)=topo(1,e2)
        elem%subcnc(2)%array(2)=topo(2,e2)
        elem%subcnc(2)%array(3)=jnode 

        elem%subcnc(2)%array(4:5)=elem%subcnc(2)%array(1:2)+NNDRL/2 ! upper surf nodes
        elem%subcnc(2)%array(6)=elem%subcnc(2)%array(3)+NNDFL/2 

        subglbcnc(2)%array(:)=elem%node_connec(elem%subcnc(2)%array(:))
        
        ! sub elm 3 connec
        elem%subcnc(3)%array(1)=topo(1,e3)
        elem%subcnc(3)%array(2)=topo(2,e3)
        elem%subcnc(3)%array(3)=jnode  
        
        elem%subcnc(3)%array(4:5)=elem%subcnc(3)%array(1:2)+NNDRL/2 ! upper surf nodes
        elem%subcnc(3)%array(6)=elem%subcnc(3)%array(3)+NNDFL/2

        subglbcnc(3)%array(:)=elem%node_connec(elem%subcnc(3)%array(:))

        ! create sub elements
        call set(elem%subelem(1),eltype='wedge', matkey=elem%bulkmat, ply_angle=elem%ply_angle, glbcnc=subglbcnc(1)%array)
        call set(elem%subelem(2),eltype='wedge', matkey=elem%bulkmat, ply_angle=elem%ply_angle, glbcnc=subglbcnc(2)%array)
        call set(elem%subelem(3),eltype='wedge', matkey=elem%bulkmat, ply_angle=elem%ply_angle, glbcnc=subglbcnc(3)%array)
     


    case (4) !- two edges cracked
     
        ibe1=min(ifedg(1),ifedg(2),ifedg(3),ifedg(4))   ! local edge index of 1st broken edge
        do i=1, 4
            if(ifedg(i)<=NEDGE/2) ibe2=max(ibe2,ifedg(i))! local edge index of 2nd broken edge
        end do  
        
        

        ! ibe1 must be between 1 to 3, and ibe2 between 2 to 4, with ibe2 > ibe1
        if(ibe1>=NEDGE/2 .or. ibe1==0 .or. ibe2<=1 .or. ibe2<=ibe1) then
            write(msg_file,*) 'something wrong in fBrick update subcnc case nfailedge=4'
            call exit_function
        end if

        if(edgstat(ibe1)==cohcrack .or. edgstat(ibe2)==cohcrack) iscoh=.true.

        
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
                        write(msg_file,*)'wrong 2nd broken edge in update subcnc fBrick'
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
                        write(msg_file,*)'wrong 2nd broken edge in update subcnc fBrick'
                        call exit_function
                end select
            case(3)
                if(ibe2==4) then
                    nbulk=4
                    e1=3; e2=4; e3=1; e4=2
                else
                    write(msg_file,*)'wrong 2nd broken edge in update subcnc fBrick'
                    call exit_function                    
                end if    
            case default
                write(msg_file,*)'wrong broken edge in update subcnc fBrick'
                call exit_function
        end select
        
        
        select case(nbulk)
            case(2)
            ! two brick bulk subdomains
                if(iscoh) then
                    nsub=3  ! two brick and one coh2d
                else
                    nsub=2  ! only two brick
                end if
                if(allocated(elem%subelem)) deallocate(elem%subelem)
                if(allocated(elem%subcnc)) deallocate(elem%subcnc)
                if(allocated(subglbcnc)) deallocate(subglbcnc)
                allocate(elem%subelem(nsub))
                allocate(elem%subcnc(nsub))
                allocate(subglbcnc(nsub))
                do j=1, nsub
                    allocate(elem%subcnc(j)%array(8))
                    allocate(subglbcnc(j)%array(8))
                    elem%subcnc(j)%array=0
                    subglbcnc(j)%array=0
                end do
                
                ! find the smaller glb fl. node on the broken edge
                if(elem%node_connec(topo(3,e1))<elem%node_connec(topo(4,e1))) then
                    jnode1=topo(3,e1)
                else
                    jnode1=topo(4,e1)
                end if
                
                ! find the smaller glb fl. node on the broken edge
                if(elem%node_connec(topo(3,e3))<elem%node_connec(topo(4,e3))) then
                    jnode3=topo(3,e3)
                else
                    jnode3=topo(4,e3)
                end if
                
                ! sub elm 1 connec
                elem%subcnc(1)%array(1)=topo(1,e1)
                elem%subcnc(1)%array(2)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(1)%array(2)=jnode1
                elem%subcnc(1)%array(3)=topo(4,e3); if(edgstat(e3)<cohcrack) elem%subcnc(1)%array(3)=jnode3
                elem%subcnc(1)%array(4)=topo(2,e3)
                
                elem%subcnc(1)%array(5)=elem%subcnc(1)%array(1)+NNDRL/2
                elem%subcnc(1)%array(6)=elem%subcnc(1)%array(2)+NNDFL/2
                elem%subcnc(1)%array(7)=elem%subcnc(1)%array(3)+NNDFL/2
                elem%subcnc(1)%array(8)=elem%subcnc(1)%array(4)+NNDRL/2
                
                subglbcnc(1)%array(:)=elem%node_connec(elem%subcnc(1)%array(:))
                
                ! sub elm 2 connec
                elem%subcnc(2)%array(1)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(2)%array(1)=jnode1
                elem%subcnc(2)%array(2)=topo(2,e1)
                elem%subcnc(2)%array(3)=topo(1,e3)
                elem%subcnc(2)%array(4)=topo(3,e3); if(edgstat(e3)<cohcrack) elem%subcnc(2)%array(4)=jnode3
                
                elem%subcnc(2)%array(5)=elem%subcnc(2)%array(1)+NNDFL/2
                elem%subcnc(2)%array(6)=elem%subcnc(2)%array(2)+NNDRL/2
                elem%subcnc(2)%array(7)=elem%subcnc(2)%array(3)+NNDRL/2
                elem%subcnc(2)%array(8)=elem%subcnc(2)%array(4)+NNDFL/2

                subglbcnc(2)%array(:)=elem%node_connec(elem%subcnc(2)%array(:))
                
                ! create sub bulk elements
                
                call set(elem%subelem(1),eltype='brick',matkey=elem%bulkmat,ply_angle=elem%ply_angle,&
                &glbcnc=subglbcnc(1)%array)
                call set(elem%subelem(2),eltype='brick',matkey=elem%bulkmat,ply_angle=elem%ply_angle,&
                &glbcnc=subglbcnc(2)%array)
                
                
                if(iscoh) then
                
                ! sub elm 3 connec
                elem%subcnc(3)%array(2)=topo(4,e3); if(edgstat(e3)<cohcrack) elem%subcnc(3)%array(2)=jnode3
                elem%subcnc(3)%array(1)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(1)=jnode1
                
                elem%subcnc(3)%array(3)=elem%subcnc(3)%array(2)+NNDFL/2
                elem%subcnc(3)%array(4)=elem%subcnc(3)%array(1)+NNDFL/2
                
                elem%subcnc(3)%array(6)=topo(3,e3); if(edgstat(e3)<cohcrack) elem%subcnc(3)%array(6)=jnode3
                elem%subcnc(3)%array(5)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(5)=jnode1
                
                elem%subcnc(3)%array(7)=elem%subcnc(3)%array(6)+NNDFL/2
                elem%subcnc(3)%array(8)=elem%subcnc(3)%array(5)+NNDFL/2
                
                subglbcnc(3)%array(:)=elem%node_connec(elem%subcnc(3)%array(:))
                
                
                call set(elem%subelem(3),eltype='coh3d8',matkey=elem%cohmat,glbcnc=subglbcnc(3)%array)
                
                
                end if
                
            case(4)
            ! four triangular bulk subdomains
                if(iscoh) then
                    nsub=5
                else
                    nsub=4
                end if
                if(allocated(elem%subelem)) deallocate(elem%subelem)
                if(allocated(elem%subcnc)) deallocate(elem%subcnc)
                if(allocated(subglbcnc)) deallocate(subglbcnc)
                allocate(elem%subelem(nsub))
                allocate(elem%subcnc(nsub))
                allocate(subglbcnc(nsub))
                do j=1, 4   ! bulk elem cnc
                    allocate(elem%subcnc(j)%array(6))
                    allocate(subglbcnc(j)%array(6))
                    elem%subcnc(j)%array=0
                    subglbcnc(j)%array=0
                end do
                if(iscoh) then
                    allocate(elem%subcnc(5)%array(8))
                    allocate(subglbcnc(5)%array(8))
                    elem%subcnc(5)%array=0
                    subglbcnc(5)%array=0
                end if
                
                ! find the smaller glb fl. node on the broken edge
                if(elem%node_connec(topo(3,e1))<elem%node_connec(topo(4,e1))) then
                    jnode1=topo(3,e1)
                else
                    jnode1=topo(4,e1)
                end if
                
                ! find the smaller glb fl. node on the broken edge
                if(elem%node_connec(topo(3,e2))<elem%node_connec(topo(4,e2))) then
                    jnode2=topo(3,e2)
                else
                    jnode2=topo(4,e2)
                end if
                
                ! sub elm 1 connec
                elem%subcnc(1)%array(1)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(1)%array(1)=jnode1
                elem%subcnc(1)%array(2)=topo(1,e2)
                elem%subcnc(1)%array(3)=topo(3,e2); if(edgstat(e2)<cohcrack) elem%subcnc(1)%array(3)=jnode2

                elem%subcnc(1)%array(4)=elem%subcnc(1)%array(1)+NNDFL/2 ! upper surf nodes
                elem%subcnc(1)%array(5)=elem%subcnc(1)%array(2)+NNDRL/2
                elem%subcnc(1)%array(6)=elem%subcnc(1)%array(3)+NNDFL/2
                
                subglbcnc(1)%array(:)=elem%node_connec(elem%subcnc(1)%array(:))
                
                ! sub elm 2 connec
                elem%subcnc(2)%array(1)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(2)%array(1)=jnode2
                elem%subcnc(2)%array(2)=topo(1,e3)
                elem%subcnc(2)%array(3)=topo(2,e3)
                
                elem%subcnc(2)%array(4)=elem%subcnc(2)%array(1)+NNDFL/2 ! upper surf nodes
                elem%subcnc(2)%array(5)=elem%subcnc(2)%array(2)+NNDRL/2
                elem%subcnc(2)%array(6)=elem%subcnc(2)%array(3)+NNDRL/2
                
                subglbcnc(2)%array(:)=elem%node_connec(elem%subcnc(2)%array(:))
                
                ! sub elm 3 connec
                elem%subcnc(3)%array(1)=topo(1,e4)
                elem%subcnc(3)%array(2)=topo(2,e4)
                elem%subcnc(3)%array(3)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(3)=jnode1
                
                elem%subcnc(3)%array(4)=elem%subcnc(3)%array(1)+NNDRL/2 ! upper surf nodes
                elem%subcnc(3)%array(5)=elem%subcnc(3)%array(2)+NNDRL/2
                elem%subcnc(3)%array(6)=elem%subcnc(3)%array(3)+NNDFL/2
                
                subglbcnc(3)%array(:)=elem%node_connec(elem%subcnc(3)%array(:))
                
                ! sub elm 4 connec
                elem%subcnc(4)%array(1)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(4)%array(1)=jnode1
                elem%subcnc(4)%array(2)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(4)%array(2)=jnode2
                elem%subcnc(4)%array(3)=topo(2,e3)

                elem%subcnc(4)%array(4)=elem%subcnc(4)%array(1)+NNDFL/2 ! upper surf nodes
                elem%subcnc(4)%array(5)=elem%subcnc(4)%array(2)+NNDFL/2
                elem%subcnc(4)%array(6)=elem%subcnc(4)%array(3)+NNDRL/2
                
                subglbcnc(4)%array(:)=elem%node_connec(elem%subcnc(4)%array(:))
                
                ! create sub bulk elements
                call set(elem%subelem(1),eltype='wedge',matkey=elem%bulkmat,&
                &ply_angle=elem%ply_angle,glbcnc=subglbcnc(1)%array)
                call set(elem%subelem(2),eltype='wedge',matkey=elem%bulkmat,&
                &ply_angle=elem%ply_angle,glbcnc=subglbcnc(2)%array)
                call set(elem%subelem(3),eltype='wedge',matkey=elem%bulkmat,&
                &ply_angle=elem%ply_angle,glbcnc=subglbcnc(3)%array)
                call set(elem%subelem(4),eltype='wedge',matkey=elem%bulkmat,&
                &ply_angle=elem%ply_angle,glbcnc=subglbcnc(4)%array)
                
                if(iscoh) then
                
                ! sub elm 5 connec
                elem%subcnc(5)%array(5)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(5)%array(5)=jnode1
                elem%subcnc(5)%array(6)=topo(3,e2); if(edgstat(e2)<cohcrack) elem%subcnc(5)%array(6)=jnode2
                
                elem%subcnc(5)%array(7)=elem%subcnc(5)%array(6)+NNDFL/2
                elem%subcnc(5)%array(8)=elem%subcnc(5)%array(5)+NNDFL/2
                
                elem%subcnc(5)%array(1)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(5)%array(1)=jnode1
                elem%subcnc(5)%array(2)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(5)%array(2)=jnode2   
                
                elem%subcnc(5)%array(3)=elem%subcnc(5)%array(2)+NNDFL/2
                elem%subcnc(5)%array(4)=elem%subcnc(5)%array(1)+NNDFL/2
                
                subglbcnc(5)%array(:)=elem%node_connec(elem%subcnc(5)%array(:))
                
                call set(elem%subelem(5),eltype='coh3d8', matkey=elem%cohmat, glbcnc=subglbcnc(5)%array)
                
                end if
                
            case default
                write(msg_file,*)'wrong nbulk in update subcnc fBrick'
                call exit_function
        end select

      
!        case(6)
    ! do not update partition
        

!        case(8)


       
    case default
       write(msg_file,*) 'WARNING: fBrick update subcnc case selection default!'
       
    end select
    
    
    ! deallocate local array
    
    if(allocated(subglbcnc)) deallocate(subglbcnc)


end subroutine update_subcnc



pure subroutine integrate_assemble_subelems (elem, nodes, lam_mat, coh_mat, &
& K_matrix, F_vector, istat, emsg, nofailure)
!---------------------------------------------------------------------!
!       integrate and assemble sub element system arrays
!---------------------------------------------------------------------!     
  ! - passed in variables   
  type(fBrick_element), intent(inout)	    :: elem
  real(DP), 	intent(inout)			:: K_matrix(:,:), F_vector(:)
    logical, intent(in)                     :: nofailure
  ! - local variables
  real(DP),	allocatable           	:: Ki(:,:), Fi(:)   ! sub_elem K matrix and F vector
  integer, 		allocatable 			:: dofcnc(:)
    integer :: i,j,l
    character(len=eltypelength) ::  subeltype
    
    i=0;j=0;l=0
    
    
    
    ! empty K and F for reuse
    K_matrix=zero; F_vector=zero
    
    ! integrate sub elements and assemble into global matrix
    do i=1, size(elem%subelem)
    
        call extract(elem%subelem(i),eltype=subeltype)
        
        if(subeltype=='coh3d6' .or. subeltype=='coh3d8') then
            call integrate(elem%subelem(i),Ki,Fi)
        else
            call integrate(elem%subelem(i),Ki,Fi,nofailure)  
        end if
        
        if(allocated(dofcnc)) deallocate(dofcnc)
        allocate(dofcnc(size(Fi))); dofcnc=0
        
        do j=1, size(elem%subcnc(i)%array) ! no. of nodes in sub elem i
            do l=1, NDIM
                ! dof indices of the jth node of sub elem i 
                dofcnc((j-1)*NDIM+l)=(elem%subcnc(i)%array(j)-1)*NDIM+l
            end do
        end do
        
        call assembleKF(K_matrix,F_vector,Ki,Fi,dofcnc)
        
        deallocate(Ki)
        deallocate(Fi)
        deallocate(dofcnc)
        
    end do 
    
    if(allocated(Ki)) deallocate(Ki)
    if(allocated(Fi)) deallocate(Fi)
    if(allocated(dofcnc)) deallocate(dofcnc)
    
end subroutine integrate_assemble




end module fBrick_element_module
