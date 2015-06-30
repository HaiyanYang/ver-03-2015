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
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    29/06/15  B. Y. Chen            Original code
!
use parameter_module,       only : NDIM, DP, ZERO, INT_ALLOC_ARRAY
use global_clock_module,    only : program_clock
use brick_element_module,   only : brick_element
use basePly_element_module, only : basePly_element
use cohCrack_element_module,only : cohCrack_element


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

  integer :: curr_status                  = 0
  real(DP):: ply_angle                    = ZERO
  integer :: node_connec(NNODE)           = 0
  integer :: edge_status_lcl(NEDGE_SURF)  = 0
  logical :: newpartition                 = .false.
  type(program_clock)                :: local_clock
  type(brick_element),   allocatable :: intact_elem
  type(basePly_element), allocatable :: subBulks(:)
  type(cohCrack_element),allocatable :: cohCrack
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



pure subroutine empty_fBrick_element (elem)

  type(fBrick_element), intent(inout) :: elem

  type(fBrick_element) :: elem_lcl

  elem = elem_lcl

end subroutine empty_fBrick_element



pure subroutine set_fBrick_element (elem, ply_angle, node_connec, istat, emsg)
! Purpose:
! to set the element ready for first use
use parameter_module,      only : MSGLENGTH, STAT_FAILURE, STAT_SUCCESS
use brick_element_module,  only : set

  type(fBrick_element),     intent(inout) :: elem
  real(DP),                 intent(in)    :: ply_angle
  integer,                  intent(in)    :: node_connec(NNODE)
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

  ! update to elem_lcl first
  elem_lcl%ply_angle   = ply_angle
  elem_lcl%node_connec = node_connec

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



pure subroutine extract_fBrick_element (elem, curr_status, edge_status_lcl, &
& intact_elem, subBulks, cohCrack)

  type(fBrick_element), intent(in)  :: elem
  integer,    optional, intent(out) :: curr_status
  integer,    optional, intent(out) :: edge_status_lcl(NEDGE_SURF)
  type(brick_element),   allocatable, optional, intent(out) :: intact_elem
  type(basePly_element), allocatable, optional, intent(out) :: subBulks(:)
  type(cohCrack_element),allocatable, optional, intent(out) :: cohCrack

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



pure subroutine integrate_fBrick_element (elem, nodes, edge_status, lam_mat, &
& coh_mat, K_matrix, F_vector, istat, emsg)
use parameter_module,         only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,&
                              & ZERO,    INTACT, TRANSITION_ELEM,              &
                              & REFINEMENT_ELEM,  CRACK_TIP_ELEM,              &
                              & CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM
use fnode_module,             only : fnode
use lamina_material_module,   only : lamina_material
use cohesive_material_module, only : cohesive_material
use global_clock_module,      only : GLOBAL_CLOCK, clock_in_sync

  type(fBrick_element),     intent(inout) :: elem
  type(fnode),              intent(inout) :: nodes(NNODE)
  integer,                  intent(inout) :: edge_status(NEDGE)
  type(lamina_material),    intent(in)    :: lam_mat
  type(cohesive_material),  intent(in)    :: coh_mat
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  !:::: local variables ::::
  ! local copy of intent inout variables
  type(fBrick_element)      :: el
  type(fnode)               :: nds(NNODE)
  integer                   :: egstatus(NEDGE)
  integer                   :: elstatus
  ! logical control variables
  logical                   :: nofailure
  ! error msg location
  character(len=MSGLENGTH)  :: msgloc

  ! initialize intent out and local variables
  allocate(K_matrix(NDOF,NDOF), F_vector(NDOF))
  K_matrix        = ZERO
  F_vector        = ZERO
  istat           = STAT_SUCCESS
  emsg            = ''
  egstatus        = 0
  elstatus        = 0
  nofailure       = .false.
  msgloc          = ' integrate, fBrick_element module'

  ! copy intent inout arg. to its local alias
  el       = elem
  nds      = nodes
  egstatus = edge_status

  ! check if last iteration has converged; if so, the elem's current partition
  ! has lead to the convergence of the last increment, and hence it is NOT a
  ! new partition but a converged partition, el%newpartition = .false.
  if (.not. clock_in_sync(GLOBAL_CLOCK, el%local_clock)) then
    el%local_clock  = GLOBAL_CLOCK
    el%newpartition = .false.
  end if


  !---------------------------------------------------------------------!
  !********** MAIN CALCULATIONS **********
  !---------------------------------------------------------------------!

  elstatuscase: select case (el%curr_status)
  !
  ! if elem has not yet reached final partition, then in each ITERATION,
  ! check for EDGE_STATUS_PARTITION of the elem first, and update newpartition:
  ! -> if el partition is UPDATED by edge status, el newpartition is set TRUE
  ! -> if el partition is UNCHANGED by edge status, el newpartition is UNCHANGED
  !
  ! after the edge status partition, calculate based on el partition status:
  ! -> if newpartition is true (current partition is not converged), then:
  !     - integrate subelems with NOfailure
  !     - do NOT check for FAILURE_CRITERION_PARTITION
  ! -> if newpartition is false (current partition has converged), then:
  !     - integrate subelems with failure
  !     - check for FAILURE_CRITERION_PARTITION, and:
  !       * if its partition is UPDATED by the failure criterion,
  !         newpartition becomes TRUE and re-integrate subelems with NOfailure
  !
  case (INTACT, TRANSITION_ELEM, REFINEMENT_ELEM, &
  &              CRACK_TIP_ELEM, CRACK_WAKE_ELEM) elstatuscase

      !***** check edge status variables *****
      ! save existing el status for comparison purpose later
      elstatus = el%curr_status
      ! update elem status according to edge status values
      ! and update the edge status values
      ! and partition the elem
      call edge_status_partition (el, nds, egstatus, istat, emsg)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call clean_up (K_matrix, F_vector)
        return
      end if
      ! if elstat is changed, then this elem has a new partition
      ! from edge status update; el%newpartition variable = TRUE
      if (elstatus /= el%curr_status) el%newpartition = .true.

      !***** select what to do based on el partition status *****
      newPartition: select case (el%newpartition)
      ! if the current elem partition is a new partition (not converged),
      ! no failure in subelems and no failure_criterion_partition is done
      ! just integrate and assemble subelems
      case (.true.) newPartition
          ! suppress subelem failure for increment of new elem partition
          nofailure = .true.
          ! integrate the subelems and assemle their K and F
          call integrate_assemble_subelems (el, nds, lam_mat, coh_mat, &
          & K_matrix, F_vector, istat, emsg, nofailure)
          if (istat == STAT_FAILURE) then
            emsg = emsg//trim(msgloc)
            call clean_up (K_matrix, F_vector)
            return
          end if
      ! if the current elem partition is a converged partition,
      ! then failure in subelems is allowed, and failure_criterion_partition
      ! needs to be checked
      case (.false.) newPartition
          ! allow failure in subelems
          nofailure = .false.
          ! integrate the subelems before failure criterion partition
          call integrate_assemble_subelems (el, nds, lam_mat, coh_mat, &
          & K_matrix, F_vector, istat, emsg, nofailure)
          if (istat == STAT_FAILURE) then
            emsg = emsg//trim(msgloc)
            call clean_up (K_matrix, F_vector)
            return
          end if
          !***** check failure criterion *****
          ! failure criterion partitions elem of any status directly into
          ! MATRIX_CRACK_ELEM partition if the failure criterion judges
          ! any subelem reaches MATRIX/FIBRE failure onset
          call failure_criterion_partition (el, nds, egstatus, istat, emsg)
          if (istat == STAT_FAILURE) then
            emsg = emsg//trim(msgloc)
            call clean_up (K_matrix, F_vector)
            return
          end if
          ! check to see if the elem partition is updated by failure criterion
          ! NOTE: any update goes straight to MATRIX_CRACK_ELEM status
          ! if updated, newpartition becomes TRUE and
          ! subelems need to be re-integrated with failure suppressed
          if (el%curr_status == MATRIX_CRACK_ELEM) then
            ! this is a new element partition not yet converged
            el%newpartition = .true.
            ! suppress failure for subelem during new partition increment
            nofailure = .true.
            ! re-integrate the subelems and re-assemle their K and F
            call integrate_assemble_subelems (el, nds, lam_mat, coh_mat, &
            & K_matrix, F_vector, istat, emsg, nofailure)
            if (istat == STAT_FAILURE) then
              emsg = emsg//trim(msgloc)
              call clean_up (K_matrix, F_vector)
              return
            end if
          end if
      end select newPartition

  ! element matrix is already failed
  ! element is already partitioned into subBulks and one cohCrack subelems,
  ! just integrate and assemble the subelems
  case (MATRIX_CRACK_ELEM) elstatuscase

      ! if the current partition is a new partition (not yet converged),
      ! suppress failure in subelems
      nofailure = el%newpartition

      ! integrate and assemble sub elems
      call integrate_assemble_subelems (el, nds, lam_mat, coh_mat, &
      & K_matrix, F_vector, istat, emsg, nofailure)


  case default elstatuscase

      istat = STAT_FAILURE
      emsg  = 'unexpected elstatus value in'//trim(msgloc)
      call clean_up (K_matrix, F_vector)
      return

  end select elstatuscase


  !---------------------------------------------------------------------!
  !********** END MAIN CALCULATIONS **********
  !---------------------------------------------------------------------!

  ! update to intent inout args before successful return
  elem        = el
  nodes       = nds
  edge_status = egstatus

  return


  contains


  pure subroutine clean_up (K_matrix, F_vector)
    real(DP), intent(inout) :: K_matrix(:,:), F_vector(:)
    K_matrix = ZERO
    F_vector = ZERO
  end subroutine clean_up

end subroutine integrate_fBrick_element



pure subroutine edge_status_partition (elem, nodes, edge_status, istat, emsg)
use parameter_module,      only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,  &
                          & REAL_ALLOC_ARRAY, ZERO, INTACT,                   &
                          & TRANSITION_EDGE, REFINEMENT_EDGE,                 &
                          & CRACK_TIP_EDGE,  COH_CRACK_EDGE,                  &
                          & TRANSITION_ELEM, REFINEMENT_ELEM,                 &
                          & CRACK_TIP_ELEM,  CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM
use fnode_module,          only : fnode, extract, update
use global_toolkit_module, only : crack_elem_cracktip2d

  ! passed-in variables
  type(fBrick_element),     intent(inout) :: elem
  type(fnode),              intent(inout) :: nodes(NNODE)
  integer,                  intent(inout) :: edge_status(NEDGE)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variable
  character(len=MSGLENGTH)  :: msgloc
  integer                   :: elstatus, eledgestatus_lcl(NEDGE_SURF)
  type(REAL_ALLOC_ARRAY)    :: coords(NNODE)
  integer                   :: nfailedge
  integer                   :: ifailedge(NEDGE_SURF)
  real(DP)                  :: crackpoint1(2), crackpoint2(2)
  real(DP)                  :: botsurf_coords(2,NEDGE_SURF)
  integer                   :: jbe1, jbe2, jnode
  integer                   :: loc(1)
  integer                   :: i

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
  istat             = STAT_SUCCESS
  emsg              = ''
  msgloc            = ' edge status partition'
  elstatus          = 0
  eledgestatus_lcl  = 0
  nfailedge         = 0
  ifailedge         = 0
  crackpoint1       = ZERO
  crackpoint2       = ZERO
  botsurf_coords    = ZERO
  jbe1              = 0
  jbe2              = 0
  jnode             = 0
  loc               = 0
  i = 0

  ! check input validity
  ! it is assumed that the top surf. edge status
  ! are the same as the bottom surf. edge status
  if ( any(edge_status(1:NEDGE_SURF) /= edge_status(NEDGE_SURF+1:NEDGE)) ) then
    istat = STAT_FAILURE
    emsg  = 'incompatible top and bot surf edge status values'//trim(msgloc)
    return
  end if

  ! extract elem status variable
  elstatus          = elem%curr_status
  ! extract elem edge status variables stored in the elem component
  eledgestatus_lcl  = elem%edge_status_lcl
  ! extract nodal coords from passed in nodes array
  do i = 1, NNODE
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
  ! 1st broken edge, jbe1, is already found and stored in eledgestatus_lcl,
  ! find the other edge crossed by the crack line starting from jbe1 crack point
  ! store its edge index jbe2
  ! partition elem to trans. elem or ref. elem based on status of jbe1 & jbe2
  !
  ! - if elem is of OTHER STATUS,
  ! then the two broken edges (jbe1 & jbe2) are already stored in eledgestatus_lcl;
  ! just update their edge status with the passed in edge status array
  ! then partition elem to other status based on status of jbe1 & jbe2
  !
  jbe1jbe2: select case (elstatus)

  case (INTACT, TRANSITION_ELEM) jbe1jbe2
      ! if elem is INTACT, return of no edge fails; otherwise, find jbe1
      ! as the most critical edge
      if (elstatus == INTACT) then
        if ( count(edge_status(1:NEDGE_SURF) > INTACT) == 0 ) then
          ! no broken edge, elem remains intact, do nothing
          return
        else
          ! find the most damaged edge index jbe1 and ignore the rest
          !jbe1 = maxloc(edge_status(1:NEDGE_SURF)) will not work, RHS is array
          loc  = maxloc(edge_status(1:NEDGE_SURF))
          jbe1 = loc(1)
        end if
      ! if elem is TRANSITION, find jbe1 as the stored edge in eledgestatus_lcl
      else
        ! check if there's only 1 broken edge stored
        if ( count(eledgestatus_lcl > INTACT) /= 1 ) then
          istat = STAT_FAILURE
          emsg  = 'unexpected no. of failed edges for case &
          & elstatus = TRANSITION ELEM'//trim(msgloc)
          return
        end if
        ! find the broken edge index
        loc  = maxloc(eledgestatus_lcl)
        jbe1 = loc(1)
        ! check if the edge status is transition edge
        if (eledgestatus_lcl(jbe1) /= TRANSITION_EDGE) then
          istat = STAT_FAILURE
          emsg  = 'unexpected edge status of failed edge for case &
          & elstatus = TRANSITION ELEM'//trim(msgloc)
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
        botsurf_coords(:,i) = coords(i)%array(1:2)
      end do
      ! find the other edge crossed by the crack line and the other crack point
      call crack_elem_cracktip2d (cracktip_point = crackpoint1,               &
      &                      cracktip_edge_index = jbe1,                      &
      &                                    nedge = NEDGE_SURF,                &
      &                              crack_angle = elem%ply_angle,            &
      &                                   coords = botsurf_coords,            &
      &                           nodes_on_edges = ENDNODES_ON_BOT_EDGES,     &
      &                                    istat = istat,                     &
      &                                     emsg = emsg,                      &
      &                         edge_crack_point = crackpoint2,               &
      &                         crack_edge_index = jbe2)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        return
      end if
      ! in the future, consider adding here the algorithm to find the top surf
      ! jbe2 crackpoint according to the fracture plane

  case (REFINEMENT_ELEM, CRACK_TIP_ELEM, &
  &     CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM) jbe1jbe2
      ! once partitioned, only check the status of stored failed edges;
      ! other failed edges are ignored
      do i = 1, NEDGE_SURF
        ! find the failed edges in stored edge status array
        if (eledgestatus_lcl(i) /= INTACT) then
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
        return
      end if
      ! assign jbe1 and jbe2
      jbe1 = ifailedge(1)
      jbe2 = ifailedge(2)

  case default jbe1jbe2
      istat = STAT_FAILURE
      emsg  = 'unexpected elstatus value'//trim(msgloc)
      return

  end select jbe1jbe2

  ! extract the edge status of jbe1 and jbe2
  ! from passed in edge_status array and update to eledgestatus_lcl
  eledgestatus_lcl(jbe1) = edge_status(jbe1)
  eledgestatus_lcl(jbe2) = edge_status(jbe2)

  ! update edge status and elem status, based on jbe1 and jbe2 status
  call update_edge_elem_status (eledgestatus_lcl, elstatus, jbe1, jbe2)

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
      ! update the coords of two fl. nodes on edge jbe2
      ! of both top and bot surfaces, assuming a perpendicular matrix crack
      jnode = NODES_ON_BOT_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_BOT_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_TOP_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_TOP_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10
    end if
10  nd_update_err : if (istat == STAT_FAILURE) then
      emsg = emsg//trim(msgloc)
      return
    end if nd_update_err

    ! update to passed in edge_status array (both surfs)
    edge_status(jbe1) = eledgestatus_lcl(jbe1)
    edge_status(jbe2) = eledgestatus_lcl(jbe2)
    edge_status(jbe1+NEDGE_SURF) = eledgestatus_lcl(jbe1)
    edge_status(jbe2+NEDGE_SURF) = eledgestatus_lcl(jbe2)

    ! update to elem components
    elem%curr_status     = elstatus
    elem%edge_status_lcl = eledgestatus_lcl

    ! update elem partition
    call partition_elem (elem, istat, emsg)
    if (istat == STAT_FAILURE) then
      emsg = emsg//trim(msgloc)
      return
    end if

  end if

  ! exit program
  return


  contains


  pure subroutine update_edge_elem_status (eledgestatus_lcl, elstatus, &
  & ibe1, ibe2)
  ! update the edge status of the two broken edges and elem status
    integer, intent(inout) :: eledgestatus_lcl(NEDGE_SURF), elstatus
    integer, intent(in)    :: ibe1, ibe2

    integer :: jbe1, jbe2
    jbe1 = 0
    jbe2 = 0

    ! store the more critical edge of the two in jbe1
    if (eledgestatus_lcl(ibe2) > eledgestatus_lcl(ibe1)) then
      jbe1 = ibe2
      jbe2 = ibe1
    else
      jbe1 = ibe1
      jbe2 = ibe2
    end if

    select case (eledgestatus_lcl(jbe1))
    ! if 1st broken edge is transition edge, the 2nd can only be transition
    ! edge or intact edge
    case (TRANSITION_EDGE)
        if (eledgestatus_lcl(jbe2) == TRANSITION_EDGE) then
          eledgestatus_lcl(jbe1) = REFINEMENT_EDGE
          eledgestatus_lcl(jbe2) = REFINEMENT_EDGE
          elstatus = REFINEMENT_ELEM
        else if (eledgestatus_lcl(jbe2) == INTACT) then
          elstatus = TRANSITION_ELEM
        end if
    ! if 1st broken edge is refinement edge, the 2nd can only be refinement,
    ! transition or intact edge
    case (REFINEMENT_EDGE)
        if (eledgestatus_lcl(jbe2) == REFINEMENT_EDGE .or. &
        &   eledgestatus_lcl(jbe2) == TRANSITION_EDGE) then
          elstatus = REFINEMENT_ELEM
        else if (eledgestatus_lcl(jbe2) == INTACT) then
          eledgestatus_lcl(jbe2) = TRANSITION_EDGE
          elstatus = REFINEMENT_ELEM
        end if
    ! if 1st broken edge is crack tip edge, the 2nd can only be crack tip,
    ! refinement, transition or intact edge
    case (CRACK_TIP_EDGE)
        if (eledgestatus_lcl(jbe2) == CRACK_TIP_EDGE) then
          eledgestatus_lcl(jbe1) = COH_CRACK_EDGE
          eledgestatus_lcl(jbe2) = COH_CRACK_EDGE
          elstatus = MATRIX_CRACK_ELEM
        else if (eledgestatus_lcl(jbe2) == REFINEMENT_EDGE) then
          elstatus = CRACK_TIP_ELEM
        else if (eledgestatus_lcl(jbe2) == TRANSITION_EDGE .or. &
        &        eledgestatus_lcl(jbe2) == INTACT) then
          eledgestatus_lcl(jbe2) = REFINEMENT_EDGE
          elstatus = CRACK_TIP_ELEM
        end if
    ! if 1st broken edge is coh crack edge, the 2nd can only be coh crack,
    ! crack tip, refinement, transition or intact edge
    case (COH_CRACK_EDGE)
        if (eledgestatus_lcl(jbe2) == COH_CRACK_EDGE) then
          elstatus = MATRIX_CRACK_ELEM
        else if (eledgestatus_lcl(jbe2) == CRACK_TIP_EDGE) then
          elstatus = CRACK_WAKE_ELEM
        else if (eledgestatus_lcl(jbe2) == REFINEMENT_EDGE .or. &
        &        eledgestatus_lcl(jbe2) == TRANSITION_EDGE .or. &
        &        eledgestatus_lcl(jbe2) == INTACT) then
          eledgestatus_lcl(jbe2) = CRACK_TIP_EDGE
          elstatus = CRACK_WAKE_ELEM
        end if
    end select

  end subroutine update_edge_elem_status


end subroutine edge_status_partition



pure subroutine failure_criterion_partition (elem, nodes, edge_status, &
& istat, emsg)
! Purpose:
! this subroutine updates elem status & partition to MATRIX_CRACK_ELEM
! if any sub elem is nolonger INTACT
use parameter_module,       only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                          & REAL_ALLOC_ARRAY, ZERO, INTACT,                   &
                          & TRANSITION_EDGE, COH_CRACK_EDGE,                  &
                          & TRANSITION_ELEM, REFINEMENT_ELEM,                 &
                          & CRACK_TIP_ELEM,  CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM
use fnode_module,           only : fnode, extract, update
use brick_element_module,   only : extract
use basePly_element_module, only : extract
use cohCrack_element_module,only : extract
use global_toolkit_module,  only : crack_elem_centroid2d, crack_elem_cracktip2d

  ! passed-in variables
  type(fBrick_element),     intent(inout) :: elem
  type(fnode),              intent(inout) :: nodes(NNODE)
  integer,                  intent(inout) :: edge_status(NEDGE)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variable
  character(len=MSGLENGTH)  :: msgloc
  integer                   :: subelstatus
  type(REAL_ALLOC_ARRAY)    :: coords(NNODE)
  integer                   :: nfailedge
  integer                   :: ifailedge(NEDGE_SURF)
  real(DP)                  :: crackpoints(2,2), crackpoint1(2), crackpoint2(2)
  real(DP)                  :: botsurf_coords(2,NEDGE_SURF)
  integer                   :: crackedges(2), jbe1, jbe2, jnode
  integer                   :: loc(1)
  logical                   :: failed
  integer                   :: i

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
  istat             = STAT_SUCCESS
  emsg              = ''
  msgloc            = ' failure criterion partition'
  subelstatus       = 0
  nfailedge         = 0
  ifailedge         = 0
  crackpoints       = ZERO
  crackpoint1       = ZERO
  crackpoint2       = ZERO
  botsurf_coords    = ZERO
  crackedges        = 0
  jbe1              = 0
  jbe2              = 0
  jnode             = 0
  loc               = 0
  failed            = .false.
  i = 0

  ! check input validity
  ! it is assumed that the top surf. edge status
  ! are the same as the bottom surf. edge status
  if ( any(edge_status(1:NEDGE_SURF) /= edge_status(NEDGE_SURF+1:NEDGE)) ) then
    istat = STAT_FAILURE
    emsg  = 'incompatible top and bot surf edge status values'//trim(msgloc)
    return
  end if

  ! extract nodal coords from passed in nodes array
  do i = 1, NNODE
    call extract (nodes(i), x=coords(i)%array)
  end do


  !**** MAIN CALCULATIONS ****
  !
  ! - if elem is INTACT, check intact_elem fstat:
  ! if NO matrix/fibre failure onset, elem remain intact;
  ! if matrix/fibre failure onset, crack elem from centroid:
  !   * find the edge indices of the two cracked edges (jbe1 and jbe2)
  !
  ! - if elem is TRANSITION ELEM, check subBulks fstat:
  ! if NO matrix/fibre failure onset, elem partition remain unchanged;
  ! if matrix/fibre failure onset in any subBulk, then partition elem from
  ! existing crack point:
  !   * 1st crack edge, jbe1, is already found and stored in eledgestatus_lcl,
  !   * find the 2nd crack edge crossed by the crack line starting from crack
  !     point on edge jbe1, and store its edge index jbe2
  !
  ! - if elem is REFINEMENT/CRACK TIP ELEM, check subBulks fstat:
  ! if NO matrix/fibre failure onset, elem partition remain unchanged;
  ! if matrix/fibre failure onset in any subBulk, then partition elem from
  ! the two crack edges (jbe1 & jbe2) stored in eledgestatus_lcl;
  !
  ! - if elem is CRACK WAKE ELEM, check subBulks and cohCrack fstat:
  ! if NO matrix/fibre/coh failure onset, elem partition remain unchanged;
  ! if matrix/fibre/coh failure onset in any subBulk/cohCrack, then partition
  ! elem from the two crack edges (jbe1 & jbe2) stored in eledgestatus_lcl;

  ! set failed to false before the update in main loop
  failed = .false.

  select case (elem%curr_status)

  ! elem is INTACT, check intact_elem fstat:
  case (INTACT)
      ! check expected no. of cracked edges
      if ( count(elem%edge_status_lcl > INTACT) /= 0 ) then
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of cracked edges for case &
        &INTACT in'//trim(msgloc)
        return
      end if
      ! proceed if no error
      call extract(elem%intact_elem, fstat=subelstatus)
      ! if matrix/fibre failure onset, crack elem from centroid and
      ! find the edge indices of the two cracked edges (jbe1 and jbe2)
      ! the crack_elem_centroid2d subroutine is used for this purpose
      ! some inputs need to be prepared
      if (subelstatus /= INTACT) then
        ! update failed to true
        failed = .true.
        ! 2D nodal coords of the bottom quad surface (first 4 nodes of the elem)
        do i = 1, NEDGE_SURF
          botsurf_coords(:,i) = coords(i)%array(1:2)
        end do
        ! use the crack_elem_centroid2d subroutine to find 2 cross points and 2
        ! crack edges of the elem, with crack passing elem centroid
        call crack_elem_centroid2d (nedge = NEDGE_SURF,             &
        &                     crack_angle = elem%ply_angle,         &
        &                          coords = botsurf_coords,         &
        &                  nodes_on_edges = ENDNODES_ON_BOT_EDGES,  &
        &                           istat = istat,                  &
        &                            emsg = emsg,                   &
        &               edge_crack_points = crackpoints,            &
        &              crack_edge_indices = crackedges)
        if (istat == STAT_FAILURE) then
          emsg = emsg//trim(msgloc)
          return
        end if
        jbe1 = crackedges(1)
        jbe2 = crackedges(2)
        ! in the future, considering extracting also the matrix crack angle
        ! w.r.t the shell plane, and partition both the top and the bot surfs
      end if

  ! elem is TRANSITION ELEM, check subBulks fstat
  case (TRANSITION_ELEM)
      ! check expected no. of cracked edges
      if ( count(elem%edge_status_lcl /= INTACT) /= 1 ) then
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of cracked edges for case &
        &TRANSITION ELEM in'//trim(msgloc)
        return
      end if
      ! index of the 1st failed edge
      !jbe1  = maxloc(elem%edge_status_lcl) will not work, RHS returns array
      loc = maxloc(elem%edge_status_lcl)
      jbe1 = loc(1)
      ! check expected status of cracked edge
      if (elem%edge_status_lcl(jbe1) /= TRANSITION_EDGE) then
        istat = STAT_FAILURE
        emsg  = 'unexpected cracked edge status for case &
        &TRANSITION ELEM in'//trim(msgloc)
        return
      end if
      ! proceed if no error encounterd
      do i = 1, size(elem%subBulks)
        call extract(elem%subBulks(i), fstat=subelstatus)
        if (subelstatus /= INTACT) then
          ! update failed to true
          failed = .true.
          exit
        end if
      end do
      ! if matrix/fibre failure onset in any subBulk, then partition elem from
      ! existing crack point; find the other edge, jbe2, crossed by the crack
      ! line which passes the crack point on edge jbe1.
      ! use the subroutine crack_elem_cracktip2d for this purpose
      if (failed) then
        ! some inputs of the subroutine need to be prepared
        ! index of a fl. node on this edge
        jnode = NODES_ON_BOT_EDGES(3,jbe1)
        ! coords of this fl. node is the existing crack point coords
        crackpoint1(:) = coords(jnode)%array(1:2)
        ! 2D nodal coords of the bottom quad surface (first 4 nodes of the elem)
        do i = 1, NEDGE_SURF
          botsurf_coords(:,i) = coords(i)%array(1:2)
        end do
        ! find the 2nd edge crossed by the crack line and the 2nd crack point
        call crack_elem_cracktip2d (cracktip_point = crackpoint1,              &
        &                      cracktip_edge_index = jbe1,                     &
        &                                    nedge = NEDGE_SURF,               &
        &                              crack_angle = elem%ply_angle,           &
        &                                   coords = botsurf_coords,           &
        &                           nodes_on_edges = ENDNODES_ON_BOT_EDGES,    &
        &                                    istat = istat,                    &
        &                                     emsg = emsg,                     &
        &                         edge_crack_point = crackpoint2,              &
        &                         crack_edge_index = jbe2)
        if (istat == STAT_FAILURE) then
          emsg = emsg//trim(msgloc)
          return
        end if
        ! in the future, consider adding here the algorithm to find the top surf
        ! jbe2 crackpoint according to the fracture plane
      end if

  ! elem is REFINEMENT/CRACK TIP/CRACK WAKE ELEM,
  ! check subBulks fstat and cohCrack fstat (CRACK WAKE ELEM ONLY)
  case (REFINEMENT_ELEM, CRACK_TIP_ELEM, CRACK_WAKE_ELEM)
      ! extract failed edges stored in elem
      do i = 1, NEDGE_SURF
        ! find the failed edges in stored edge status array
        if (elem%edge_status_lcl(i) /= INTACT) then
          ! update total no. of failed edges
          nfailedge = nfailedge + 1
          ! update the indices of failed edges in local array
          ifailedge(nfailedge) = i
        end if
      end do
      ! check if nfailedge == 2
      if (nfailedge /= 2) then
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of failed edges for case &
        & REFINEMENT/CRACK TIP/CRACK WAKE ELEM'//trim(msgloc)
        return
      end if
      ! assign jbe1 and jbe2
      jbe1 = ifailedge(1)
      jbe2 = ifailedge(2)

      ! check if subelems are failed
      ! if elem is CRACK WAKE ELEM, check if coh crack reaches failure onset
      if (elem%curr_status == CRACK_WAKE_ELEM) then
        call extract(elem%cohCrack, fstat=subelstatus)
        if (subelstatus /= INTACT) then
          ! update failed to true
          failed = .true.
        end if
      end if
      ! if no coh crack or no failure onset in coh crack, continue to check
      ! if subBulks have reached failure onset
      if (.not. failed) then
        do i = 1, size(elem%subBulks)
          call extract(elem%subBulks(i), fstat=subelstatus)
          if (subelstatus /= INTACT) then
            ! update failed to true
            failed = .true.
            exit
          end if
        end do
      end if

  case default
      istat = STAT_FAILURE
      emsg  = 'unexpected elem curr status in'//trim(msgloc)
      return

  end select

  ! if matrix/fibre/coh failure onset in any intact_elem/subBulk/cohCrack, then:
  ! - update the passed-in nodes array, coords of fl. nodes on jbe1 & jbe2
  ! - update the passed-in edge status array, status of jbe1 & jbe2
  ! - update the edge status of the elem
  ! - update the curr status of the elem
  ! - update the partition   of the elem
  if (failed) then

    !:::: update passed-in nodes array ::::!

    ! update edges jbe1 and jbe2 fl. nodes when elem was previously intact
    ! (on both surfs), coords update to crackpoints
    if (elem%curr_status == INTACT) then
      ! update the coords of two fl. nodes on edge jbe1
      ! of both top and bot surfaces, now assuming a perpendicular matrix crack
      jnode = NODES_ON_BOT_EDGES(3,jbe1)
      coords(jnode)%array(1:2)=crackpoints(:,1)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10
      
      jnode = NODES_ON_BOT_EDGES(4,jbe1)
      coords(jnode)%array(1:2)=crackpoints(:,1)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_TOP_EDGES(3,jbe1)
      coords(jnode)%array(1:2)=crackpoints(:,1)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_TOP_EDGES(4,jbe1)
      coords(jnode)%array(1:2)=crackpoints(:,1)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      ! update the coords of two fl. nodes on edge jbe2
      ! of both top and bot surfaces, now assuming a perpendicular matrix crack
      jnode = NODES_ON_BOT_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoints(:,2)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_BOT_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoints(:,2)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_TOP_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoints(:,2)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_TOP_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoints(:,2)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

    ! update edge jbe2 fl. nodes when elem was previously transition elem
    ! (on both surfs), coords update to crackpoint2
    else if (elem%curr_status == TRANSITION_ELEM) then
      ! update the coords of two fl. nodes on edge jbe2
      ! of both top and bot surfaces, now assuming a perpendicular matrix crack
      jnode = NODES_ON_BOT_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_BOT_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_TOP_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

      jnode = NODES_ON_TOP_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      call update(nodes(jnode), istat, emsg, x=coords(jnode)%array)
      if (istat == STAT_FAILURE) goto 10

    end if
10  nd_update_err : if (istat == STAT_FAILURE) then
      emsg = emsg//trim(msgloc)
      return
    end if nd_update_err

    !:::: update passed-in edge_status array (both surfs) ::::!
    edge_status(jbe1)             = COH_CRACK_EDGE
    edge_status(jbe2)             = COH_CRACK_EDGE
    edge_status(jbe1+NEDGE_SURF)  = COH_CRACK_EDGE
    edge_status(jbe2+NEDGE_SURF)  = COH_CRACK_EDGE

    !:::: update elem components ::::!
    elem%edge_status_lcl(jbe1)    = COH_CRACK_EDGE
    elem%edge_status_lcl(jbe2)    = COH_CRACK_EDGE
    elem%curr_status              = MATRIX_CRACK_ELEM

    !:::: update elem partition ::::!
    call partition_elem (elem, istat, emsg)
    if (istat == STAT_FAILURE) then
      emsg = emsg//trim(msgloc)
      return
    end if

  end if

end subroutine failure_criterion_partition



pure subroutine partition_elem (el, istat, emsg)
use parameter_module, only : MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,         &
                      & INT_ALLOC_ARRAY, ELTYPELENGTH, INTACT,              &
                      & TRANSITION_ELEM, REFINEMENT_ELEM,   CRACK_TIP_ELEM, &
                      & CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM, CRACK_TIP_EDGE, &
                      & COH_CRACK_EDGE,  TRANSITION_EDGE
use global_toolkit_module,  only : partition_quad_elem
use basePly_element_module, only : set
use cohCrack_element_module,only : set

  ! passed-in variables
  type(fBrick_element),     intent(inout) :: el
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variables
  character(len=MSGLENGTH)           :: msgloc
  integer                            :: nfailedge, ifailedge(NEDGE_SURF)
  integer                            :: nds_botsurf(4,NEDGE_SURF)
  integer                            :: nds_topsurf(4,NEDGE_SURF)
  integer                            :: flnode1, flnode2
  type(INT_ALLOC_ARRAY), allocatable :: subBulks_nds_top(:)
  type(INT_ALLOC_ARRAY), allocatable :: subBulks_nds_bot(:)
  type(INT_ALLOC_ARRAY), allocatable :: subBulks_glb_connec(:)
  integer                            :: nsub, subnnode
  character(len=ELTYPELENGTH)        :: subeltype
  integer                            :: cohcrack_nds_top(4)
  integer                            :: cohcrack_nds_bot(4)
  integer                            :: cohcrack_glb_connec(8)
  integer :: i, j

  ! initialize intent out and local variables
  istat               = STAT_SUCCESS
  emsg                = ''
  msgloc              = ' partition_elem'
  nfailedge           = 0
  ifailedge           = 0
  nds_botsurf         = 0
  nds_topsurf         = 0
  flnode1             = 0
  flnode2             = 0
  nsub                = 0
  subnnode            = 0
  subeltype           = ''
  cohcrack_nds_top    = 0
  cohcrack_nds_bot    = 0
  cohcrack_glb_connec = 0
  i = 0; j = 0

  ! deallocate intact elem first
  if(allocated(el%intact_elem)) deallocate(el%intact_elem)

  ! find the failed edges in stored edge status array
  do i = 1, NEDGE_SURF
    if (el%edge_status_lcl(i) /= INTACT) then
      ! update total no. of failed edges
      nfailedge = nfailedge + 1
      ! update the indices of failed edges in local array
      ifailedge(nfailedge) = i
    end if
  end do

  ! copy the top & bot edges topology to local arrays for changes later
  nds_botsurf = NODES_ON_BOT_EDGES
  nds_topsurf = NODES_ON_TOP_EDGES

  !********** define nds_bot/topsurf based on elem's curr status **********
  select case (el%curr_status)

  case (TRANSITION_ELEM)
  ! partition into subBulks only, for the purpose of mesh refinement
  ! only 1 fl. node is used on the crack edge to form the partition
      ! check edge status
      if (nfailedge /= 1) then
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of cracked edges, trans elem,'//trim(msgloc)
        return
      end if
      if (el%edge_status_lcl(ifailedge(1)) /= TRANSITION_EDGE) then
        istat = STAT_FAILURE
        emsg  = 'unexpected status of cracked edge, trans elem,'//trim(msgloc)
        return
      end if
      ! make the 2nd fl. node same as the 1st on the crack edge, so that when
      ! partitioned, a refinement node, rather than a crack point, is formed
      ! the 1st fl. node is the one with a smaller global node index
      flnode1 = NODES_ON_BOT_EDGES(3,ifailedge(1))
      flnode2 = NODES_ON_BOT_EDGES(4,ifailedge(1))
      if (el%node_connec(flnode2) < el%node_connec(flnode1)) then
        nds_botsurf(3:4,ifailedge(1)) = flnode2
      else
        nds_botsurf(3:4,ifailedge(1)) = flnode1
      end if
      flnode1 = NODES_ON_TOP_EDGES(3,ifailedge(1))
      flnode2 = NODES_ON_TOP_EDGES(4,ifailedge(1))
      if (el%node_connec(flnode2) < el%node_connec(flnode1)) then
        nds_topsurf(3:4,ifailedge(1)) = flnode2
      else
        nds_topsurf(3:4,ifailedge(1)) = flnode1
      end if

  case (REFINEMENT_ELEM, CRACK_TIP_ELEM)
  ! partition into subBulks only, for the purpose of mesh refinement
  ! only 1 fl. node is used on each crack edge to form the partition
      ! check edge status
      if (nfailedge /= 2) then
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of cracked edges, ref/tip elem,'//trim(msgloc)
        return
      end if
      if (maxval(el%edge_status_lcl(ifailedge(1:2))) > CRACK_TIP_EDGE .or. &
      &   minval(el%edge_status_lcl(ifailedge(1:2))) < TRANSITION_EDGE) then
        istat = STAT_FAILURE
        emsg  = 'unexpected status of cracked edges, ref/tip elem,'//trim(msgloc)
        return
      end if
      ! make the 2nd fl. node same as the 1st on crack edges, so that when
      ! partitioned, a refinement line, rather than an actual crack, is formed
      ! the 1st fl. node is the one with a smaller global node index
      do j = 1, 2
          flnode1 = NODES_ON_BOT_EDGES(3,ifailedge(j))
          flnode2 = NODES_ON_BOT_EDGES(4,ifailedge(j))
          if (el%node_connec(flnode2) < el%node_connec(flnode1)) then
            nds_botsurf(3:4,ifailedge(j)) = flnode2
          else
            nds_botsurf(3:4,ifailedge(j)) = flnode1
          end if
          flnode1 = NODES_ON_TOP_EDGES(3,ifailedge(j))
          flnode2 = NODES_ON_TOP_EDGES(4,ifailedge(j))
          if (el%node_connec(flnode2) < el%node_connec(flnode1)) then
            nds_topsurf(3:4,ifailedge(j)) = flnode2
          else
            nds_topsurf(3:4,ifailedge(j)) = flnode1
          end if
      end do

  case (CRACK_WAKE_ELEM)
  ! partition into TWO subBulks and ONE cohCrack, with one crack edge being
  ! the CRACK_TIP_EDGE and the other the COH_CRACK_EDGE
  ! only 1 fl. node is used on the CRACK_TIP_EDGE to form the partition
      ! check edge status
      if (nfailedge /= 2) then
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of cracked edges, wake elem,'//trim(msgloc)
        return
      end if
      if (maxval(el%edge_status_lcl(ifailedge(1:2))) /= COH_CRACK_EDGE .or. &
      &   minval(el%edge_status_lcl(ifailedge(1:2))) /= CRACK_TIP_EDGE) then
        istat = STAT_FAILURE
        emsg  = 'unexpected status of cracked edges, wake elem,'//trim(msgloc)
        return
      end if
      ! make the 2nd fl. node same as the 1st on crack tip edge, so that
      ! when partitioned, it is a refinement node, not an edge crack point
      ! the 1st fl. node is the one with a smaller global node index
      do j = 1, 2
        if (el%edge_status_lcl(ifailedge(j)) == CRACK_TIP_EDGE) then
          flnode1 = NODES_ON_BOT_EDGES(3,ifailedge(j))
          flnode2 = NODES_ON_BOT_EDGES(4,ifailedge(j))
          if (el%node_connec(flnode2) < el%node_connec(flnode1)) then
            nds_botsurf(3:4,ifailedge(j)) = flnode2
          else
            nds_botsurf(3:4,ifailedge(j)) = flnode1
          end if
          flnode1 = NODES_ON_TOP_EDGES(3,ifailedge(j))
          flnode2 = NODES_ON_TOP_EDGES(4,ifailedge(j))
          if (el%node_connec(flnode2) < el%node_connec(flnode1)) then
            nds_topsurf(3:4,ifailedge(j)) = flnode2
          else
            nds_topsurf(3:4,ifailedge(j)) = flnode1
          end if
          exit
        end if
      end do

  case (MATRIX_CRACK_ELEM)
  ! partition into TWO subBulks and ONE cohCrack, no modification of
  ! nds_topsurf and nds_botsurf arrays
      ! check edge status
      if (nfailedge /= 2) then
        istat = STAT_FAILURE
        emsg  = 'unexpected no. of cracked edges, crack elem,'//trim(msgloc)
        return
      end if
      if (any( el%edge_status_lcl(ifailedge(1:2)) /= COH_CRACK_EDGE )) then
        istat = STAT_FAILURE
        emsg  = 'unexpected status of cracked edges, crack elem,'//trim(msgloc)
        return
      end if

  case default
      istat = STAT_FAILURE
      emsg  = 'unexpected elem curr status in'//trim(msgloc)
      return

  end select

  !********** partition the top and bottom surfaces first **********
  ! call the partition quad elem subroutine from global toolkit module
  ! pass the TOP surface topology nds_topsurf into the subroutine,
  ! this will partition the element top surface into 2D subBulk elems,
  ! whose nodes (in lcl indices) will be stored in subBulks_nds_top;
  ! nodes of the 2D coh crack will be stored in cohcrack_nds_top.
  call partition_quad_elem (nds_topsurf, ifailedge, subBulks_nds_top, &
  & istat, emsg, cohcrack_nds_top)
  if (istat == STAT_FAILURE) then
    emsg = emsg//trim(msgloc)
    call clean_up (subBulks_glb_connec, subBulks_nds_top, subBulks_nds_bot)
    return
  end if
  ! pass the BOTTOM surface topology nds_botsurf into the subroutine,
  ! this will partition the element bottom surface into 2D subBulk elems,
  ! whose nodes (in lcl indices) will be stored in subBulks_nds_bot;
  ! nodes of the 2D coh crack will be stored in cohcrack_nds_bot.
  call partition_quad_elem (nds_botsurf, ifailedge, subBulks_nds_bot, &
  & istat, emsg, cohcrack_nds_bot)
  if (istat == STAT_FAILURE) then
    emsg = emsg//trim(msgloc)
    call clean_up (subBulks_glb_connec, subBulks_nds_top, subBulks_nds_bot)
    return
  end if
  ! ** NOTE **
  ! the ifailedge array is applicable to both top and bottom edges; so
  ! the top and bottom surfaces are partitioned in exactly the same way, i.e.,
  ! same number of sub elems, same sub elem types.


  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! based on the top/bot surface partitions in subBulks_nds_top/bot,
  ! allocate and populate the el%subBulks_nodes array and set el%subBulks
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  ! extract no. of sub elems
  nsub = size(subBulks_nds_top)

  ! allocate nsub no. of subBulks, subBulks_nodes and
  ! subBulks_glb_connec
  if(allocated(el%subBulks))         deallocate(el%subBulks)
  if(allocated(el%subBulks_nodes))   deallocate(el%subBulks_nodes)
  if(allocated(subBulks_glb_connec)) deallocate(subBulks_glb_connec)
  allocate(el%subBulks(nsub))
  allocate(el%subBulks_nodes(nsub))
  allocate(subBulks_glb_connec(nsub))

  ! allocate & define subBulks
  do_subbulks: do j = 1, nsub
      ! determine no. of nodes in sub elem j
      subnnode = 2 * size(subBulks_nds_top(j)%array)
      ! determine sub elem type based on subnnode
      select case (subnnode)
        case (6)
          subeltype = 'wedge'
        case (8)
          subeltype = 'brick'
        case default
          istat = STAT_FAILURE
          emsg  = 'unexpected no. of nodes for sub elem in'//trim(msgloc)
          call clean_up (subBulks_glb_connec, subBulks_nds_top, subBulks_nds_bot)
          return
      end select
      ! allocate connec for sub elems
      allocate(el%subBulks_nodes(j)%array(subnnode))
      allocate(subBulks_glb_connec(j)%array(subnnode))
      ! initialize these arrays
      el%subBulks_nodes(j)%array   = 0
      subBulks_glb_connec(j)%array = 0
      ! populate lcl connec
      ! copy top surf. lcl connec
      el%subBulks_nodes(j)%array( subnnode/2 + 1 : subnnode   ) = &
      &  subBulks_nds_top(j)%array(:)
      ! copy bot surf. lcl connec
      el%subBulks_nodes(j)%array(              1 : subnnode/2 ) = &
      &  subBulks_nds_bot(j)%array(:)
      ! populate glb connec through el%node_connec
      subBulks_glb_connec(j)%array = el%node_connec(el%subBulks_nodes(j)%array)
      ! set this sub element
      call set(el%subBulks(j), eltype=trim(subeltype),               &
      & connec=subBulks_glb_connec(j)%array, ply_angle=el%ply_angle, &
      & istat=istat, emsg=emsg)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call clean_up (subBulks_glb_connec, subBulks_nds_top, subBulks_nds_bot)
        return
      end if
  end do do_subbulks

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! based on the el's current partition status,
  ! allocate and populate the el%cohCrack_nodes array and set el%cohCrack
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  if (el%curr_status == CRACK_WAKE_ELEM .or. &
  &   el%curr_status == MATRIX_CRACK_ELEM) then
  ! cohCrack is present for this partition
      ! allocate arrays for the coh crack
      if (.not. allocated(el%cohCrack)) allocate(el%cohCrack)
      if (.not. allocated(el%cohCrack_nodes)) then
        allocate(el%cohCrack_nodes)
        ! allocate 8 nodes for cohCrack coh crack subelem
        allocate(el%cohCrack_nodes%array(8))
        el%cohCrack_nodes%array = 0
      end if
      ! assign values to the cohcrack nodes from top&bot suf 2d cohcrack nodes
!      el%cohCrack_nodes%array(1:2) = cohcrack_nds_top(1:2)
!      el%cohCrack_nodes%array(3:4) = cohcrack_nds_bot(2:1)
!      el%cohCrack_nodes%array(5:6) = cohcrack_nds_top(4:3)
!      el%cohCrack_nodes%array(7:8) = cohcrack_nds_bot(3:4)

      el%cohCrack_nodes%array(1) = cohcrack_nds_top(1)
      el%cohCrack_nodes%array(2) = cohcrack_nds_top(2)
      el%cohCrack_nodes%array(3) = cohcrack_nds_bot(2)
      el%cohCrack_nodes%array(4) = cohcrack_nds_bot(1)
      el%cohCrack_nodes%array(5) = cohcrack_nds_top(4)
      el%cohCrack_nodes%array(6) = cohcrack_nds_top(3)
      el%cohCrack_nodes%array(7) = cohcrack_nds_bot(3)
      el%cohCrack_nodes%array(8) = cohcrack_nds_bot(4)
      ! populate glb connec through el%node_connec
      cohcrack_glb_connec = el%node_connec(el%cohCrack_nodes%array)
      ! set this sub element
      call set(el%cohCrack, connec=cohCrack_glb_connec, istat=istat, emsg=emsg)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call clean_up (subBulks_glb_connec, subBulks_nds_top, subBulks_nds_bot)
        return
      end if
  end if


  ! deallocate local alloc arrays before successful return
  call clean_up (subBulks_glb_connec, subBulks_nds_top, subBulks_nds_bot)

  return

  contains

  pure subroutine clean_up (subBulks_glb_connec, subBulks_nds_top, &
  & subBulks_nds_bot)
    type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subBulks_glb_connec(:)
    type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subBulks_nds_top(:)
    type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subBulks_nds_bot(:)

    if (allocated(subBulks_glb_connec)) deallocate(subBulks_glb_connec)
    if (allocated(subBulks_nds_top))    deallocate(subBulks_nds_top)
    if (allocated(subBulks_nds_bot))    deallocate(subBulks_nds_bot)

  end subroutine clean_up

end subroutine partition_elem



pure subroutine integrate_assemble_subelems (elem, nodes, lam_mat, coh_mat, &
& K_matrix, F_vector, istat, emsg, nofailure)
! Purpose :
! integrate and assemble sub element system arrays
use parameter_module, only : MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,         &
                      & DP, NDIM, ZERO, ONE, PENALTY_STIFFNESS, INTACT,     &
                      & TRANSITION_ELEM, REFINEMENT_ELEM,   CRACK_TIP_ELEM, &
                      & CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM, CRACK_TIP_EDGE
use fnode_module,             only : fnode
use lamina_material_module,   only : lamina_material
use cohesive_material_module, only : cohesive_material
use brick_element_module,     only : integrate
use basePly_element_module,   only : integrate
use cohCrack_element_module,    only : integrate
use global_toolkit_module,    only : assembleKF

  ! - passed in variables
  type(fBrick_element),     intent(inout) :: elem
  type(fnode),              intent(in)    :: nodes(NNODE)
  type(lamina_material),    intent(in)    :: lam_mat
  type(cohesive_material),  intent(in)    :: coh_mat
  real(DP),                 intent(out)   :: K_matrix(NDOF,NDOF), F_vector(NDOF)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure

  ! - local variables
  character(len=MSGLENGTH)    :: msgloc
  logical                     :: nofail
  integer                     :: ifailedge(NEDGE_SURF), nfailedge
  integer                     :: node1, node2
  integer                     :: i, j

  ! initialize intent out and local variables
  K_matrix  = ZERO
  F_vector  = ZERO
  istat     = STAT_SUCCESS
  emsg      = ''
  msgloc    = ' integrate_assemble_subelems'
  nofail    = .false.
  ifailedge = 0
  nfailedge = 0
  node1     = 0
  node2     = 0
  i = 0
  j = 0

  if (present(nofailure)) nofail = nofailure

  ! find the failed edges in stored edge status array
  do i = 1, NEDGE_SURF
    if (elem%edge_status_lcl(i) /= INTACT) then
      ! update total no. of failed edges
      nfailedge = nfailedge + 1
      ! update the indices of failed edges in local array
      ifailedge(nfailedge) = i
    end if
  end do

  !********** select integration procedures based on elem status **********
  select case (elem%curr_status)

  case (INTACT)
      ! check if intact elem is allocated
      if (.not. allocated(elem%intact_elem)) then
        istat = STAT_FAILURE
        emsg  = 'intact elem NOT allocated for case INTACT in'//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if
      ! integrate and assemble the intact elem
      call integrate_assemble_intact_elem (elem, nodes, lam_mat, &
      & K_matrix, F_vector, istat, emsg, nofail)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if

  case (TRANSITION_ELEM, REFINEMENT_ELEM, CRACK_TIP_ELEM)
  ! integrate the subBulks and constrain fl. nodes of crack edges
      ! check if subBulks are allocated
      if (.not. allocated(elem%subBulks) .or. &
      &   .not. allocated(elem%subBulks_nodes)) then
        istat = STAT_FAILURE
        emsg  = 'subBulks NOT allocated for trans/ref/tip elem in'//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if
      ! integrate and assmeble subBulks
      call integrate_assemble_subBulks (elem, nodes, lam_mat, &
      & K_matrix, F_vector, istat, emsg, nofail)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if
      ! constrain fl. nodes on cracked edges
      if (elem%curr_status == TRANSITION_ELEM) then
        ! impose constrain btw 2 fl. nodes on the 1 crack edge of both surfs
        node1 = NODES_ON_BOT_EDGES(3, ifailedge(1))
        node2 = NODES_ON_BOT_EDGES(4, ifailedge(1))
        call tie_two_nodes (K_matrix, node1, node2)
        node1 = NODES_ON_TOP_EDGES(3, ifailedge(1))
        node2 = NODES_ON_TOP_EDGES(4, ifailedge(1))
        call tie_two_nodes (K_matrix, node1, node2)
      else
        ! impose constrain btw 2 fl. nodes on the 2 crack edges of both surfs
        node1 = NODES_ON_BOT_EDGES(3, ifailedge(1))
        node2 = NODES_ON_BOT_EDGES(4, ifailedge(1))
        call tie_two_nodes (K_matrix, node1, node2)
        node1 = NODES_ON_BOT_EDGES(3, ifailedge(2))
        node2 = NODES_ON_BOT_EDGES(4, ifailedge(2))
        call tie_two_nodes (K_matrix, node1, node2)
        node1 = NODES_ON_TOP_EDGES(3, ifailedge(1))
        node2 = NODES_ON_TOP_EDGES(4, ifailedge(1))
        call tie_two_nodes (K_matrix, node1, node2)
        node1 = NODES_ON_TOP_EDGES(3, ifailedge(2))
        node2 = NODES_ON_TOP_EDGES(4, ifailedge(2))
        call tie_two_nodes (K_matrix, node1, node2)
      end if

  case (CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM)
  ! integrate both subBulks and cohCrack and
  ! constrain fl. nodes of crack_tip_edge for wake elem
      ! check if subBulks are allocated
      if (.not. allocated(elem%subBulks) .or. &
      &   .not. allocated(elem%subBulks_nodes)) then
        istat = STAT_FAILURE
        emsg  = 'subBulks NOT allocated for wake/crack elem in'//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if
      ! check if cohCrack is allocated
      if (.not. allocated(elem%cohCrack) .or. &
      &   .not. allocated(elem%cohCrack_nodes)) then
        istat = STAT_FAILURE
        emsg  = 'cohCrack NOT allocated for wake/crack elem in'//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if
      ! integrate and assmeble subBulks
      call integrate_assemble_subBulks (elem, nodes, lam_mat, &
      & K_matrix, F_vector, istat, emsg, nofail)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if
      ! integrate and assmeble cohCrack
      call integrate_assemble_cohCrack (elem, nodes, coh_mat, &
      & K_matrix, F_vector, istat, emsg, nofail)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if
      ! constrain fl. nodes on crack tip edge in case of wake elem
      if (elem%curr_status == CRACK_WAKE_ELEM) then
        do j=1, 2
          if ( elem%edge_status_lcl(ifailedge(j)) == CRACK_TIP_EDGE ) then
            ! impose constrain btw 2 fl. nodes on the crack edge of both surfs
            node1 = NODES_ON_BOT_EDGES(3, ifailedge(j))
            node2 = NODES_ON_BOT_EDGES(4, ifailedge(j))
            call tie_two_nodes (K_matrix, node1, node2)
            node1 = NODES_ON_TOP_EDGES(3, ifailedge(j))
            node2 = NODES_ON_TOP_EDGES(4, ifailedge(j))
            call tie_two_nodes (K_matrix, node1, node2)
            exit
          end if
        end do
      end if

  case default
      istat = STAT_FAILURE
      emsg  = 'unexpect elem curr status in'//trim(msgloc)
      call zeroKF (K_matrix, F_vector)
      return
  end select

  return


  contains


  pure subroutine zeroKF (K_matrix, F_vector)
    real(DP), intent(inout) :: K_matrix(:,:), F_vector(:)
    K_matrix = ZERO
    F_vector = ZERO
  end subroutine zeroKF


  pure subroutine tie_two_nodes (K_matrix, node1, node2)
    real(DP), intent(inout) :: K_matrix(:,:)
    integer,  intent(in)    :: node1, node2

    real(DP) :: Ktie(2,2)
    integer :: i, nrow, ncol

    Ktie = ZERO
    i    = 0
    nrow = 0
    ncol = 0

    ! form tie constrain stiffness matrix
    ! [Ktie] * {u_node1, u_node2} = {ZERO, ZERO}
    Ktie(1,1) =  ONE
    Ktie(1,2) = -ONE
    Ktie(2,1) = -ONE
    Ktie(2,2) =  ONE

    ! use the penalty stiffness defined in parameter module
    Ktie = PENALTY_STIFFNESS * Ktie

    ! add the Ktie terms to the corresponding K_matrix terms

    do i = 1, NDIM
      nrow = (node1-1)*NDIM+i
      ncol = (node1-1)*NDIM+i
      K_matrix(nrow,ncol) = K_matrix(nrow,ncol) + Ktie(1,1)
    end do

    do i = 1, NDIM
      nrow = (node1-1)*NDIM+i
      ncol = (node2-1)*NDIM+i
      K_matrix(nrow,ncol) = K_matrix(nrow,ncol) + Ktie(1,2)
    end do

    do i = 1, NDIM
      nrow = (node2-1)*NDIM+i
      ncol = (node1-1)*NDIM+i
      K_matrix(nrow,ncol) = K_matrix(nrow,ncol) + Ktie(2,1)
    end do

    do i = 1, NDIM
      nrow = (node2-1)*NDIM+i
      ncol = (node2-1)*NDIM+i
      K_matrix(nrow,ncol) = K_matrix(nrow,ncol) + Ktie(2,2)
    end do

  end subroutine tie_two_nodes


  pure subroutine clean_up (Ki, Fi)
    real(DP), allocatable, intent(inout) :: Ki(:,:), Fi(:)
    if(allocated(Ki)) deallocate(Ki)
    if(allocated(Fi)) deallocate(Fi)
  end subroutine clean_up


  pure subroutine integrate_assemble_intact_elem (elem, nodes, lam_mat, &
  & K_matrix, F_vector, istat, emsg, nofail)
    type(fBrick_element),     intent(inout) :: elem
    type(fnode),              intent(in)    :: nodes(NNODE)
    type(lamina_material),    intent(in)    :: lam_mat
    real(DP),                 intent(out)   :: K_matrix(NDOF,NDOF)
    real(DP),                 intent(out)   :: F_vector(NDOF)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    logical,                  intent(in)    :: nofail

    real(DP), allocatable :: Ki(:,:), Fi(:)

    ! only integrate the intact elem
    call integrate(elem%intact_elem, nodes(INTACT_ELEM_NODES), lam_mat, &
    & Ki, Fi, istat, emsg, nofail)
    if (istat == STAT_FAILURE) then
      call clean_up (Ki, Fi)
      return
    end if

    ! assemble the intact elem K and F
    call assembleKF(K_matrix, F_vector, Ki, Fi, INTACT_ELEM_NODES, NDIM, &
    & istat, emsg)
    if (istat == STAT_FAILURE) then
      call clean_up (Ki, Fi)
      return
    end if

    call clean_up (Ki, Fi)
  end subroutine integrate_assemble_intact_elem


  pure subroutine integrate_assemble_subBulks (elem, nodes, lam_mat, &
  & K_matrix, F_vector, istat, emsg, nofail)
    type(fBrick_element),     intent(inout) :: elem
    type(fnode),              intent(in)    :: nodes(NNODE)
    type(lamina_material),    intent(in)    :: lam_mat
    real(DP),                 intent(out)   :: K_matrix(NDOF,NDOF)
    real(DP),                 intent(out)   :: F_vector(NDOF)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    logical,                  intent(in)    :: nofail

    real(DP), allocatable :: Ki(:,:), Fi(:)
    integer :: isub

    do isub = 1, size(elem%subBulks)

        ! call the integrate procedure asoc. with the subelem type
        call integrate (elem%subBulks(isub),      &
        & nodes(elem%subBulks_nodes(isub)%array), &
        & lam_mat, Ki, Fi, istat, emsg, nofail)
        if (istat == STAT_FAILURE) then
          call clean_up (Ki, Fi)
          return
        end if

        ! assemble the sub elem K and F
        call assembleKF(K_matrix, F_vector, Ki, Fi, &
        & elem%subBulks_nodes(isub)%array, NDIM, istat, emsg)
        if (istat == STAT_FAILURE) then
          call clean_up (Ki, Fi)
          return
        end if

    end do

    call clean_up (Ki, Fi)
  end subroutine integrate_assemble_subBulks


  pure subroutine integrate_assemble_cohCrack (elem, nodes, coh_mat, &
  & K_matrix, F_vector, istat, emsg, nofail)
    type(fBrick_element),     intent(inout) :: elem
    type(fnode),              intent(in)    :: nodes(NNODE)
    type(cohesive_material),  intent(in)    :: coh_mat
    real(DP),                 intent(out)   :: K_matrix(NDOF,NDOF)
    real(DP),                 intent(out)   :: F_vector(NDOF)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    logical,                  intent(in)    :: nofail

    real(DP), allocatable :: Ki(:,:), Fi(:)

    ! integrate coh crack
    call integrate (elem%cohCrack, nodes(elem%cohCrack_nodes%array), &
    & coh_mat, Ki, Fi, istat, emsg, nofail)
    if (istat == STAT_FAILURE) then
      call clean_up (Ki, Fi)
      return
    end if

    ! assemble the coh crack K and F
    call assembleKF(K_matrix, F_vector, Ki, Fi, elem%cohCrack_nodes%array, &
    & NDIM, istat, emsg)
    if (istat == STAT_FAILURE) then
      call clean_up (Ki, Fi)
      return
    end if

    call clean_up (Ki, Fi)
  end subroutine integrate_assemble_cohCrack




end subroutine integrate_assemble_subelems




end module fBrick_element_module
