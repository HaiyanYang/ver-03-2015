module fBrickPly_elem_module
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
!  topological definition of INTACT element, type brickPly_elem:
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
use parameter_module,       only : NDIM, DP, ZERO, INT_ALLOC_ARRAY, NST_STANDARD, NST_COHESIVE
use global_clock_module,    only : program_clock
use brickPly_elem_module,   only : brickPly_elem
use abstPly_elem_module,    only : abstPly_elem
use coh8Crack_elem_module,  only : coh8Crack_elem


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

type, public :: fBrickPly_elem
  private

  integer :: curr_status                  = 0
  integer :: edge_status_lcl(NEDGE_SURF)  = 0
  logical :: newpartition                 = .false.
  type(program_clock)                :: local_clock
  type(brickPly_elem),   allocatable :: intact_elem
  type(abstPly_elem),    allocatable :: subBulks(:)
  type(coh8Crack_elem),  allocatable :: cohCrack
  type(INT_ALLOC_ARRAY), allocatable :: subBulks_nodes(:)
  type(INT_ALLOC_ARRAY), allocatable :: cohCrack_nodes

end type fBrickPly_elem

interface set
    module procedure set_fBrickPly_elem
end interface

interface integrate
    module procedure integrate_fBrickPly_elem
end interface

interface extract
    module procedure extract_fBrickPly_elem
end interface


public :: set, integrate, extract



contains



pure subroutine set_fBrickPly_elem (elem)
  type(fBrickPly_elem), intent(inout)  :: elem
  
  allocate(elem%intact_elem)

end subroutine set_fBrickPly_elem



pure subroutine extract_fBrickPly_elem (elem, curr_status, edge_status_lcl, &
& bulks_nodes, crack_nodes, bulks_stress, bulks_strain, bulks_df, crack_tau, crack_delta, crack_dm)
! used for output of this elem
use brickPly_elem_module,   only : extract
use abstPly_elem_module,    only : extract
use coh8Crack_elem_module,  only : extract

  type(fBrickPly_elem), intent(in)  :: elem
  integer,    optional, intent(out) :: curr_status
  integer,    optional, intent(out) :: edge_status_lcl(NEDGE_SURF)
  
  integer, allocatable, optional, intent(out) :: bulks_nodes(:,:)
  integer, allocatable, optional, intent(out) :: crack_nodes(:)
  real(DP),allocatable, optional, intent(out) :: bulks_stress(:,:)
  real(DP),allocatable, optional, intent(out) :: bulks_strain(:,:)
  real(DP),allocatable, optional, intent(out) :: bulks_df(:)
  real(DP),             optional, intent(out) :: crack_tau(NST_COHESIVE)
  real(DP),             optional, intent(out) :: crack_delta(NST_COHESIVE)
  real(DP),             optional, intent(out) :: crack_dm
  
  integer :: nsub, i, subnnd
  
  if(present(curr_status))     curr_status     = elem%curr_status

  if(present(edge_status_lcl)) edge_status_lcl = elem%edge_status_lcl
  
  if(present(bulks_nodes)) then
    ! if subBulks are present, then assign their nodes to bulks_nodes array
    if (allocated(elem%subBulks)) then
      nsub = size(elem%subBulks)
      ! allocate 8 nodes per sub elem
      allocate(bulks_nodes(8,nsub))
      do i = 1, nsub
        subnnd = size(elem%subBulks_nodes(i)%array)
        ! copy the node no. of subelem to arg. array
        bulks_nodes(1:subnnd,i) = elem%subBulks_nodes(i)%array(:)
        ! if less than 8 nodes, then fill the rest nodes with the last node no.
        if (subnnd < 8) bulks_nodes(subnnd+1:8,i) = bulks_nodes(subnnd,i)
      end do
    ! if not, allocate and assign intact elem nodes to bulks_nodes array
    else
      ! store the intact elem nodes
      allocate(bulks_nodes(8,1))
      bulks_nodes(:,1) = INTACT_ELEM_NODES
    end if
  end if
  
  if(present(crack_nodes)) then
    if (allocated(elem%cohCrack_nodes)) then
      allocate(crack_nodes(8))
      crack_nodes = elem%cohCrack_nodes%array
    end if
  end if
  
  if(present(bulks_stress)) then
    if (allocated(elem%subBulks)) then
      nsub = size(elem%subBulks)
      allocate(bulks_stress(NST_STANDARD,nsub))
      do i = 1, nsub
        call extract(elem%subBulks(i), stress=bulks_stress(:,i))
      end do
    else
      allocate(bulks_stress(NST_STANDARD,1))
      call extract(elem%intact_elem, stress=bulks_stress(:,1))
    end if
  end if
  
  if(present(bulks_strain)) then
    if (allocated(elem%subBulks)) then
      nsub = size(elem%subBulks)
      allocate(bulks_strain(NST_STANDARD,nsub))
      do i = 1, nsub
        call extract(elem%subBulks(i), strain=bulks_strain(:,i))
      end do
    else
      allocate(bulks_strain(NST_STANDARD,1))
      call extract(elem%intact_elem, strain=bulks_strain(:,1))
    end if
  end if
  
  if(present(bulks_df)) then
    if (allocated(elem%subBulks)) then
      nsub = size(elem%subBulks)
      allocate(bulks_df(nsub))
      do i = 1, nsub
        call extract(elem%subBulks(i), df=bulks_df(i))
      end do
    else
      allocate(bulks_df(1))
      call extract(elem%intact_elem, df=bulks_df(1))
    end if
  end if
  
  if(present(crack_tau)) then
    if (allocated(elem%cohCrack)) then
      call extract(elem%cohCrack, traction=crack_tau)
    end if
  end if
  
  if(present(crack_delta)) then
    if (allocated(elem%cohCrack)) then
      call extract(elem%cohCrack, separation=crack_delta)
    end if
  end if
  
  if(present(crack_dm)) then
    if (allocated(elem%cohCrack)) then
      call extract(elem%cohCrack, dm=crack_dm)
    end if
  end if

end subroutine extract_fBrickPly_elem



pure subroutine integrate_fBrickPly_elem (elem, nodes, ply_angle, lam_mat, coh_mat,&
& K_matrix, F_vector, istat, emsg)
use parameter_module,         only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,&
                              & ZERO,    INTACT, TRANSITION_ELEM,              &
                              & REFINEMENT_ELEM,  CRACK_TIP_ELEM,              &
                              & CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM
use fnode_module,             only : fnode
use lamina_material_module,   only : lamina_material
use cohesive_material_module, only : cohesive_material
use global_clock_module,      only : GLOBAL_CLOCK, clock_in_sync

  type(fBrickPly_elem),     intent(inout) :: elem
  type(fnode),              intent(inout) :: nodes(NNODE)
  real(DP),                 intent(in)    :: ply_angle
  type(lamina_material),    intent(in)    :: lam_mat
  type(cohesive_material),  intent(in)    :: coh_mat
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  !:::: local variables ::::
  ! local copy of intent inout variables
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
  elstatus        = 0
  nofailure       = .false.
  msgloc          = ', integrate, fBrickPly_elem module'

  ! check if last iteration has converged; if so, the elem's current partition
  ! has lead to the convergence of the last increment, and hence it is NOT a
  ! new partition but a converged partition, elem%newpartition = .false.
  if (.not. clock_in_sync(GLOBAL_CLOCK, elem%local_clock)) then
    elem%local_clock  = GLOBAL_CLOCK
    elem%newpartition = .false.
  end if


  !---------------------------------------------------------------------!
  !********** MAIN CALCULATIONS **********
  !---------------------------------------------------------------------!

  elstatuscase: select case (elem%curr_status)
  !
  ! if elem has not yet reached final partition, then in each ITERATION,
  ! check for EDGE_STATUS_PARTITION of the elem first, and update newpartition:
  ! -> if elem partition is UPDATED by edge status, elem newpartition is set TRUE
  ! -> if elem partition is UNCHANGED by edge status, elem newpartition is UNCHANGED
  !
  ! after the edge status partition, calculate based on elem partition status:
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
      ! save existing elem status for comparison purpose later
      elstatus = elem%curr_status
      ! update elem status according to edge status values
      ! and update the edge status values
      ! and partition the elem
      call edge_status_partition (elem, nodes, ply_angle, istat, emsg)
      if (istat == STAT_FAILURE) then
        emsg = trim(emsg)//trim(msgloc)
        call clean_up (K_matrix, F_vector)
        return
      end if

      !***** select what to do based on elem partition status *****
      newPartition: select case (elem%newpartition)
      ! if the current elem partition is a new partition (not converged),
      ! no failure in subelems and no failure_criterion_partition is done
      ! just integrate and assemble subelems
      case (.true.) newPartition
          ! suppress subelem failure for increment of new elem partition
          nofailure = .true.
          ! integrate the subelems and assemle their K and F
          call integrate_assemble_subelems (elem, nodes, ply_angle, lam_mat, coh_mat, &
          & K_matrix, F_vector, istat, emsg, nofailure)
          if (istat == STAT_FAILURE) then
            emsg = trim(emsg)//trim(msgloc)//'-mark 1'
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
          call integrate_assemble_subelems (elem, nodes, ply_angle, lam_mat, coh_mat, &
          & K_matrix, F_vector, istat, emsg, nofailure)
          if (istat == STAT_FAILURE) then
            emsg = trim(emsg)//trim(msgloc)//'-mark 2'
            call clean_up (K_matrix, F_vector)
            return
          end if
          !***** check failure criterion *****
          ! failure criterion partitions elem of any status directly into
          ! MATRIX_CRACK_ELEM partition if the failure criterion judges
          ! any subelem reaches MATRIX/FIBRE failure onset
          call failure_criterion_partition (elem, nodes, ply_angle, istat, emsg)
          if (istat == STAT_FAILURE) then
            emsg = trim(emsg)//trim(msgloc)
            call clean_up (K_matrix, F_vector)
            return
          end if
          ! check to see if the elem partition is updated by failure criterion
          ! NOTE: any update goes straight to MATRIX_CRACK_ELEM status
          ! if updated, newpartition becomes TRUE and
          ! subelems need to be re-integrated with failure suppressed
          if (elem%newpartition) then
            ! suppress failure for subelem during new partition increment
            nofailure = .true.
            ! re-integrate the subelems and re-assemle their K and F
            call integrate_assemble_subelems (elem, nodes, ply_angle, lam_mat, coh_mat, &
            & K_matrix, F_vector, istat, emsg, nofailure)
            if (istat == STAT_FAILURE) then
              emsg = trim(emsg)//trim(msgloc)
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
      nofailure = elem%newpartition

      ! integrate and assemble sub elems
      call integrate_assemble_subelems (elem, nodes, ply_angle, lam_mat, coh_mat, &
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

  return


  contains


  pure subroutine clean_up (K_matrix, F_vector)
    real(DP), intent(inout) :: K_matrix(:,:), F_vector(:)
    K_matrix = ZERO
    F_vector = ZERO
  end subroutine clean_up

end subroutine integrate_fBrickPly_elem



pure subroutine edge_status_partition (elem, nodes, ply_angle, istat, emsg)
use parameter_module,      only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,  &
                          & REAL_ALLOC_ARRAY, ZERO, INTACT,                   &
                          & TRANSITION_EDGE, REFINEMENT_EDGE,                 &
                          & CRACK_TIP_EDGE,  COH_CRACK_EDGE,                  &
                          & TRANSITION_ELEM, REFINEMENT_ELEM,                 &
                          & CRACK_TIP_ELEM,  CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM
use fnode_module,          only : fnode, extract, update
use global_toolkit_module, only : crack_elem_cracktip2d

  ! passed-in variables
  type(fBrickPly_elem),     intent(inout) :: elem
  type(fnode),              intent(inout) :: nodes(NNODE)
  real(DP),                 intent(in)    :: ply_angle
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local variable
  character(len=MSGLENGTH)  :: msgloc
  integer                   :: edge_status(NEDGE)
  integer                   :: elstatus, eledgestatus_lcl(NEDGE_SURF)
  type(REAL_ALLOC_ARRAY)    :: coords(NNODE)
  integer                   :: nfailedge
  integer                   :: ifailedge(NEDGE_SURF)
  real(DP)                  :: crackpoint1(2), crackpoint2(2)
  real(DP)                  :: botsurf_coords(2,NEDGE_SURF)
  real(DP)                  :: ztop, zbot
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
  edge_status       = 0
  elstatus          = 0
  eledgestatus_lcl  = 0
  nfailedge         = 0
  ifailedge         = 0
  crackpoint1       = ZERO
  crackpoint2       = ZERO
  botsurf_coords    = ZERO
  ztop              = ZERO
  zbot              = ZERO
  jbe1              = 0
  jbe2              = 0
  jnode             = 0
  loc               = 0
  i = 0

  ! extract edge_status from nodes
  do i = 1, NEDGE
    call extract (nodes(NODES_ON_EDGES(3,i)), nstat=edge_status(i))
  end do

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
  
  ! find ztop and zbot
  ! NOTE: this element assumes that top surf nodes have the same z coord: ztop
  !                            and  bot surf nodes have the same z coord: zbot
  ztop = coords( NODES_ON_TOP_EDGES(1,1) )%array(3)
  zbot = coords( NODES_ON_BOT_EDGES(1,1) )%array(3)


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
      &                              crack_angle = ply_angle,                 &
      &                                   coords = botsurf_coords,            &
      &                           nodes_on_edges = ENDNODES_ON_BOT_EDGES,     &
      &                                    istat = istat,                     &
      &                                     emsg = emsg,                      &
      &                         edge_crack_point = crackpoint2,               &
      &                         crack_edge_index = jbe2)
      if (istat == STAT_FAILURE) then
        emsg = trim(emsg)//trim(msgloc)
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

    ! update edge jbe2 fl. nodes coords when elem was previously intact
    ! or transition elem (on both surfs)
    if (elem%curr_status == INTACT .or. &
    &   elem%curr_status == TRANSITION_ELEM) then
      ! update the coords of two fl. nodes on edge jbe2
      ! of both top and bot surfaces, assuming a perpendicular matrix crack
      jnode = NODES_ON_BOT_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      coords(jnode)%array(3)  =zbot
      call update(nodes(jnode), x=coords(jnode)%array)

      jnode = NODES_ON_BOT_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      coords(jnode)%array(3)  =zbot
      call update(nodes(jnode), x=coords(jnode)%array)

      jnode = NODES_ON_TOP_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      coords(jnode)%array(3)  =ztop
      call update(nodes(jnode), x=coords(jnode)%array)

      jnode = NODES_ON_TOP_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      coords(jnode)%array(3)  =ztop
      call update(nodes(jnode), x=coords(jnode)%array)
    end if
    
    ! update bot surf fl. nodes nstat
    call update(nodes(NODES_ON_EDGES(3,jbe1)), nstat=eledgestatus_lcl(jbe1))
    call update(nodes(NODES_ON_EDGES(3,jbe2)), nstat=eledgestatus_lcl(jbe2))
    call update(nodes(NODES_ON_EDGES(4,jbe1)), nstat=eledgestatus_lcl(jbe1))
    call update(nodes(NODES_ON_EDGES(4,jbe2)), nstat=eledgestatus_lcl(jbe2))
    ! update top surf fl. nodes nstat
    call update(nodes(NODES_ON_EDGES(3,jbe1+NEDGE_SURF)), nstat=eledgestatus_lcl(jbe1))
    call update(nodes(NODES_ON_EDGES(3,jbe2+NEDGE_SURF)), nstat=eledgestatus_lcl(jbe2))
    call update(nodes(NODES_ON_EDGES(4,jbe1+NEDGE_SURF)), nstat=eledgestatus_lcl(jbe1))
    call update(nodes(NODES_ON_EDGES(4,jbe2+NEDGE_SURF)), nstat=eledgestatus_lcl(jbe2))

    ! update elem partition when changing from intact or trans. elem
    if (elem%curr_status == INTACT .or. &
    &   elem%curr_status == TRANSITION_ELEM) then
      ! update to elem components
      elem%curr_status     = elstatus
      elem%edge_status_lcl = eledgestatus_lcl
      ! update elem partition
      call partition_elem (elem, istat, emsg)
      if (istat == STAT_FAILURE) then
        emsg = trim(emsg)//trim(msgloc)
        return
      end if
      ! set newpartition logical variable to be true
      elem%newpartition = .true.
    else
      ! update to elem components
      elem%curr_status     = elstatus
      elem%edge_status_lcl = eledgestatus_lcl
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



pure subroutine failure_criterion_partition (elem, nodes, ply_angle, istat, emsg)
! Purpose:
! this subroutine updates elem status & partition to MATRIX_CRACK_ELEM
! if any sub elem is nolonger INTACT
use parameter_module,       only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                          & REAL_ALLOC_ARRAY, ZERO, INTACT,                   &
                          & TRANSITION_EDGE, COH_CRACK_EDGE,                  &
                          & TRANSITION_ELEM, REFINEMENT_ELEM,                 &
                          & CRACK_TIP_ELEM,  CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM
use fnode_module,           only : fnode, extract, update
use brickPly_elem_module,   only : extract
use abstPly_elem_module,    only : extract
use coh8Crack_elem_module,  only : extract
use global_toolkit_module,  only : crack_elem_centroid2d, crack_elem_cracktip2d

  ! passed-in variables
  type(fBrickPly_elem),     intent(inout) :: elem
  type(fnode),              intent(inout) :: nodes(NNODE)
  real(DP),                 intent(in)    :: ply_angle
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
  real(DP)                  :: ztop, zbot
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
  ztop              = ZERO
  zbot              = ZERO
  crackedges        = 0
  jbe1              = 0
  jbe2              = 0
  jnode             = 0
  loc               = 0
  failed            = .false.
  i = 0

  ! check input validity

  ! extract nodal coords from passed in nodes array
  do i = 1, NNODE
    call extract (nodes(i), x=coords(i)%array)
  end do
  
  ! find ztop and zbot
  ! NOTE: this element assumes that top surf nodes have the same z coord: ztop
  !                            and  bot surf nodes have the same z coord: zbot
  ztop = coords( NODES_ON_TOP_EDGES(1,1) )%array(3)
  zbot = coords( NODES_ON_BOT_EDGES(1,1) )%array(3)

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
  ! - if elem is REFINEMENT/CRACK TIP/CRACK WAKE ELEM, check subBulks and cohCrack fstat:
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
        &                     crack_angle = ply_angle,              &
        &                          coords = botsurf_coords,         &
        &                  nodes_on_edges = ENDNODES_ON_BOT_EDGES,  &
        &                           istat = istat,                  &
        &                            emsg = emsg,                   &
        &               edge_crack_points = crackpoints,            &
        &              crack_edge_indices = crackedges)
        if (istat == STAT_FAILURE) then
          emsg = trim(emsg)//trim(msgloc)
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
        &                              crack_angle = ply_angle,                &
        &                                   coords = botsurf_coords,           &
        &                           nodes_on_edges = ENDNODES_ON_BOT_EDGES,    &
        &                                    istat = istat,                    &
        &                                     emsg = emsg,                     &
        &                         edge_crack_point = crackpoint2,              &
        &                         crack_edge_index = jbe2)
        if (istat == STAT_FAILURE) then
          emsg = trim(emsg)//trim(msgloc)
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
      ! check if coh crack reaches failure onset
      call extract(elem%cohCrack, fstat=subelstatus)
      if (subelstatus /= INTACT) then
        ! update failed to true
        failed = .true.
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

    !:::: update passed-in fl. nodes coords ::::!

    ! update edges jbe1 and jbe2 fl. nodes when elem was previously intact
    ! (on both surfs), coords update to crackpoints
    if (elem%curr_status == INTACT) then
      ! update the coords of two fl. nodes on edge jbe1
      ! of both top and bot surfaces, now assuming a perpendicular matrix crack
      jnode = NODES_ON_BOT_EDGES(3,jbe1)
      coords(jnode)%array(1:2)=crackpoints(:,1)
      coords(jnode)%array(3)  =zbot
      call update(nodes(jnode), x=coords(jnode)%array)
      
      jnode = NODES_ON_BOT_EDGES(4,jbe1)
      coords(jnode)%array(1:2)=crackpoints(:,1)
      coords(jnode)%array(3)  =zbot
      call update(nodes(jnode), x=coords(jnode)%array)

      jnode = NODES_ON_TOP_EDGES(3,jbe1)
      coords(jnode)%array(1:2)=crackpoints(:,1)
      coords(jnode)%array(3)  =ztop
      call update(nodes(jnode), x=coords(jnode)%array)

      jnode = NODES_ON_TOP_EDGES(4,jbe1)
      coords(jnode)%array(1:2)=crackpoints(:,1)
      coords(jnode)%array(3)  =ztop
      call update(nodes(jnode), x=coords(jnode)%array)

      ! update the coords of two fl. nodes on edge jbe2
      ! of both top and bot surfaces, now assuming a perpendicular matrix crack
      jnode = NODES_ON_BOT_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoints(:,2)
      coords(jnode)%array(3)  =zbot
      call update(nodes(jnode), x=coords(jnode)%array)
      
      jnode = NODES_ON_BOT_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoints(:,2)
      coords(jnode)%array(3)  =zbot
      call update(nodes(jnode), x=coords(jnode)%array)
      
      jnode = NODES_ON_TOP_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoints(:,2)
      coords(jnode)%array(3)  =ztop
      call update(nodes(jnode), x=coords(jnode)%array)
      
      jnode = NODES_ON_TOP_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoints(:,2)
      coords(jnode)%array(3)  =ztop
      call update(nodes(jnode), x=coords(jnode)%array)
      
    ! update edge jbe2 fl. nodes when elem was previously transition elem
    ! (on both surfs), coords update to crackpoint2
    else if (elem%curr_status == TRANSITION_ELEM) then
      ! update the coords of two fl. nodes on edge jbe2
      ! of both top and bot surfaces, now assuming a perpendicular matrix crack
      jnode = NODES_ON_BOT_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      coords(jnode)%array(3)  =zbot
      call update(nodes(jnode), x=coords(jnode)%array)
      
      jnode = NODES_ON_BOT_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      coords(jnode)%array(3)  =zbot
      call update(nodes(jnode), x=coords(jnode)%array)
      
      jnode = NODES_ON_TOP_EDGES(3,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      coords(jnode)%array(3)  =ztop
      call update(nodes(jnode), x=coords(jnode)%array)
      
      jnode = NODES_ON_TOP_EDGES(4,jbe2)
      coords(jnode)%array(1:2)=crackpoint2(:)
      coords(jnode)%array(3)  =ztop
      call update(nodes(jnode), x=coords(jnode)%array)
      
    end if

    !:::: update passed-in fl. nodes nstat (both surfs) ::::!
    
    ! update bot surf fl. nodes nstat
    call update(nodes(NODES_ON_EDGES(3,jbe1)), nstat=COH_CRACK_EDGE)
    call update(nodes(NODES_ON_EDGES(3,jbe2)), nstat=COH_CRACK_EDGE)
    call update(nodes(NODES_ON_EDGES(4,jbe1)), nstat=COH_CRACK_EDGE)
    call update(nodes(NODES_ON_EDGES(4,jbe2)), nstat=COH_CRACK_EDGE)
    ! update top surf fl. nodes nstat
    call update(nodes(NODES_ON_EDGES(3,jbe1+NEDGE_SURF)), nstat=COH_CRACK_EDGE)
    call update(nodes(NODES_ON_EDGES(3,jbe2+NEDGE_SURF)), nstat=COH_CRACK_EDGE)
    call update(nodes(NODES_ON_EDGES(4,jbe1+NEDGE_SURF)), nstat=COH_CRACK_EDGE)
    call update(nodes(NODES_ON_EDGES(4,jbe2+NEDGE_SURF)), nstat=COH_CRACK_EDGE)

    !:::: update elem partition ::::!
    ! only updates partition when changing from intact and trans elem
    if (elem%curr_status == INTACT .or. &
    &   elem%curr_status == TRANSITION_ELEM) then
      !:::: update elem components ::::!
      elem%edge_status_lcl(jbe1)    = COH_CRACK_EDGE
      elem%edge_status_lcl(jbe2)    = COH_CRACK_EDGE
      elem%curr_status              = MATRIX_CRACK_ELEM
      ! update elem partition
      call partition_elem (elem, istat, emsg)
      if (istat == STAT_FAILURE) then
        emsg = trim(emsg)//trim(msgloc)
        return
      end if
      ! set newpartition logical variable to be true
      elem%newpartition = .true.
    else
      !:::: update elem components ::::!
      elem%edge_status_lcl(jbe1)    = COH_CRACK_EDGE
      elem%edge_status_lcl(jbe2)    = COH_CRACK_EDGE
      elem%curr_status              = MATRIX_CRACK_ELEM
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

  ! passed-in variables
  type(fBrickPly_elem),     intent(inout) :: el
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
  integer                            :: nsub, subnnode
  character(len=ELTYPELENGTH)        :: subeltype
  integer                            :: cohcrack_nds_top(4)
  integer                            :: cohcrack_nds_bot(4)
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
  
  !********** check correctness based on elem's curr status **********
  select case (el%curr_status)
  case (TRANSITION_ELEM)
  ! partition into subBulks only, for the purpose of mesh refinement
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
  case (REFINEMENT_ELEM, CRACK_TIP_ELEM)
  ! partition into subBulks only, for the purpose of mesh refinement
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
  case (CRACK_WAKE_ELEM)
  ! partition into TWO subBulks and ONE cohCrack, with one crack edge being
  ! the CRACK_TIP_EDGE and the other the COH_CRACK_EDGE
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
  case (MATRIX_CRACK_ELEM)
  ! partition into TWO subBulks and ONE cohCrack
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

  ! copy the top & bot edges topology to local arrays for changes later
  nds_botsurf = NODES_ON_BOT_EDGES
  nds_topsurf = NODES_ON_TOP_EDGES

  !********** partition the top and bottom surfaces first **********
  ! call the partition quad elem subroutine from global toolkit module
  ! pass the TOP surface topology nds_topsurf into the subroutine,
  ! this will partition the element top surface into 2D subBulk elems,
  ! whose nodes (in lcl indices) will be stored in subBulks_nds_top;
  ! nodes of the 2D coh crack will be stored in cohcrack_nds_top.
  call partition_quad_elem (nds_topsurf, ifailedge, subBulks_nds_top, &
  & istat, emsg, cohcrack_nds_top)
  if (istat == STAT_FAILURE) then
    emsg = trim(emsg)//trim(msgloc)
    call clean_up (subBulks_nds_top, subBulks_nds_bot)
    return
  end if
  ! pass the BOTTOM surface topology nds_botsurf into the subroutine,
  ! this will partition the element bottom surface into 2D subBulk elems,
  ! whose nodes (in lcl indices) will be stored in subBulks_nds_bot;
  ! nodes of the 2D coh crack will be stored in cohcrack_nds_bot.
  call partition_quad_elem (nds_botsurf, ifailedge, subBulks_nds_bot, &
  & istat, emsg, cohcrack_nds_bot)
  if (istat == STAT_FAILURE) then
    emsg = trim(emsg)//trim(msgloc)
    call clean_up (subBulks_nds_top, subBulks_nds_bot)
    return
  end if
  ! ** NOTE **
  ! the ifailedge array is applicable to both top and bottom edges; so
  ! the top and bottom surfaces are partitioned in exactly the same way, i.e.,
  ! same number of sub elems, same sub elem types.


  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! based on the top/bot surface partitions in subBulks_nds_top/bot,
  ! allocate and populate the el%subBulks_nodes array
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  ! extract no. of sub elems
  nsub = size(subBulks_nds_top)

  ! allocate nsub no. of subBulks, subBulks_nodes and
  ! subBulks_glb_connec
  if(allocated(el%subBulks))         deallocate(el%subBulks)
  if(allocated(el%subBulks_nodes))   deallocate(el%subBulks_nodes)
  allocate(el%subBulks(nsub))
  allocate(el%subBulks_nodes(nsub))

  ! allocate & define subBulks
  do_subbulks: do j = 1, nsub
      ! determine no. of nodes in sub elem j
      subnnode = 2 * size(subBulks_nds_top(j)%array)
      ! allocate connec for sub elems
      allocate(el%subBulks_nodes(j)%array(subnnode))
      ! initialize these arrays
      el%subBulks_nodes(j)%array   = 0
      ! populate lcl connec
      ! copy top surf. lcl connec
      el%subBulks_nodes(j)%array( subnnode/2 + 1 : subnnode   ) = &
      &  subBulks_nds_top(j)%array(:)
      ! copy bot surf. lcl connec
      el%subBulks_nodes(j)%array(              1 : subnnode/2 ) = &
      &  subBulks_nds_bot(j)%array(:)
  end do do_subbulks

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! based on the el's current partition status,
  ! allocate and populate the el%cohCrack_nodes array and set el%cohCrack
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  if (el%curr_status /= TRANSITION_ELEM ) then
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
      el%cohCrack_nodes%array(1) = cohcrack_nds_top(1)
      el%cohCrack_nodes%array(2) = cohcrack_nds_top(2)
      el%cohCrack_nodes%array(3) = cohcrack_nds_bot(2)
      el%cohCrack_nodes%array(4) = cohcrack_nds_bot(1)
      el%cohCrack_nodes%array(5) = cohcrack_nds_top(4)
      el%cohCrack_nodes%array(6) = cohcrack_nds_top(3)
      el%cohCrack_nodes%array(7) = cohcrack_nds_bot(3)
      el%cohCrack_nodes%array(8) = cohcrack_nds_bot(4)   
  end if


  ! deallocate local alloc arrays before successful return
  call clean_up (subBulks_nds_top, subBulks_nds_bot)

  return

  contains

  pure subroutine clean_up (subBulks_nds_top, subBulks_nds_bot)
    type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subBulks_nds_top(:)
    type(INT_ALLOC_ARRAY), allocatable, intent(inout) :: subBulks_nds_bot(:)

    if (allocated(subBulks_nds_top))    deallocate(subBulks_nds_top)
    if (allocated(subBulks_nds_bot))    deallocate(subBulks_nds_bot)

  end subroutine clean_up

end subroutine partition_elem



pure subroutine integrate_assemble_subelems (elem, nodes, ply_angle, lam_mat, coh_mat, &
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
use brickPly_elem_module,     only : integrate
use abstPly_elem_module,      only : integrate
use coh8Crack_elem_module,    only : integrate
use global_toolkit_module,    only : assembleKF

  ! - passed in variables
  type(fBrickPly_elem),     intent(inout) :: elem
  type(fnode),              intent(in)    :: nodes(NNODE)
  real(DP),                 intent(in)    :: ply_angle
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
      call integrate_assemble_intact_elem (elem, nodes, ply_angle, lam_mat, &
      & K_matrix, F_vector, istat, emsg, nofail)
      if (istat == STAT_FAILURE) then
        emsg = trim(emsg)//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if

  case (TRANSITION_ELEM)
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
      call integrate_assemble_subBulks (elem, nodes, ply_angle, lam_mat, &
      & K_matrix, F_vector, istat, emsg, nofail)
      if (istat == STAT_FAILURE) then
        emsg = trim(emsg)//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if

  case (REFINEMENT_ELEM, CRACK_TIP_ELEM, CRACK_WAKE_ELEM, MATRIX_CRACK_ELEM)
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
      call integrate_assemble_subBulks (elem, nodes, ply_angle, lam_mat, &
      & K_matrix, F_vector, istat, emsg, nofail)
      if (istat == STAT_FAILURE) then
        emsg = trim(emsg)//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
      end if
      ! integrate and assmeble cohCrack
      call integrate_assemble_cohCrack (elem, nodes, coh_mat, &
      & K_matrix, F_vector, istat, emsg, nofail)
      if (istat == STAT_FAILURE) then
        emsg = trim(emsg)//trim(msgloc)
        call zeroKF (K_matrix, F_vector)
        return
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


  pure subroutine integrate_assemble_intact_elem (elem, nodes, ply_angle, lam_mat, &
  & K_matrix, F_vector, istat, emsg, nofail)
    type(fBrickPly_elem),     intent(inout) :: elem
    type(fnode),              intent(in)    :: nodes(NNODE)
    real(DP),                 intent(in)    :: ply_angle
    type(lamina_material),    intent(in)    :: lam_mat
    real(DP),                 intent(out)   :: K_matrix(NDOF,NDOF)
    real(DP),                 intent(out)   :: F_vector(NDOF)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    logical,                  intent(in)    :: nofail

    real(DP), allocatable :: Ki(:,:), Fi(:)

    ! only integrate the intact elem
    call integrate(elem%intact_elem, nodes(INTACT_ELEM_NODES), ply_angle, lam_mat, &
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


  pure subroutine integrate_assemble_subBulks (elem, nodes, ply_angle, lam_mat, &
  & K_matrix, F_vector, istat, emsg, nofail)
    type(fBrickPly_elem),     intent(inout) :: elem
    type(fnode),              intent(in)    :: nodes(NNODE)
    real(DP),                 intent(in)    :: ply_angle
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
        call integrate (elem%subBulks(isub), nodes(elem%subBulks_nodes(isub)%array), &
        & ply_angle, lam_mat, Ki, Fi, istat, emsg, nofail)
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
    type(fBrickPly_elem),     intent(inout) :: elem
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




end module fBrickPly_elem_module
