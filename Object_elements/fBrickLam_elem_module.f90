module fBrickLam_elem_module
!
!  Purpose:
!    define a floating-node brick laminate element
!    
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    29/06/15  B. Y. Chen            Original code
!
use parameter_module,       only : DP, INT_ALLOC_ARRAY, ZERO
use fBrickPly_elem_module,  only : fBrickPly_elem
use fCoh8Delam_elem_module, only : fCoh8Delam_elem

implicit none
private

! parameters for ply-block elem (fbrick) and interface elem (fDelam8)
integer, parameter :: NNODE_PLYBLK = 24
integer, parameter :: NNODE_INTERF = 32
integer, parameter :: NEDGE = 8, NNDRL = 8, NNDFL = 16, NNDIN_INTERF = 8

type, public :: plyblock_layup
  real(DP) :: angle  = ZERO
  integer  :: nplies = 0
end type plyblock_layup

type, public :: fBrickLam_elem
  private
  integer :: curr_status = 0
  integer :: NPLYBLKS    = 0
  integer,               allocatable :: node_connec(:)
  type(plyblock_layup),  allocatable :: layup(:)
  type(fBrickPly_elem),  allocatable :: plyblks(:)
  type(fCoh8Delam_elem), allocatable :: interfs(:)
  type(INT_ALLOC_ARRAY), allocatable :: plyblks_nodes(:)
  type(INT_ALLOC_ARRAY), allocatable :: interfs_nodes(:)
  type(INT_ALLOC_ARRAY), allocatable :: plyblks_edges(:)
  type(INT_ALLOC_ARRAY), allocatable :: interfs_edges(:)
end type fBrickLam_elem

interface set
  module procedure set_fBrickLam_elem
end interface

interface integrate
  module procedure integrate_fBrickLam_elem
end interface

interface extract
  module procedure extract_fBrickLam_elem
end interface


public :: set, integrate, extract




contains




pure subroutine extract_fBrickLam_elem (elem, curr_status, layup, plyblks, &
& interfs)
use fBrickPly_elem_module,  only : fBrickPly_elem
use fCoh8Delam_elem_module, only : fCoh8Delam_elem

  type(fBrickLam_elem),                     intent(in)  :: elem
  integer,                           optional, intent(out) :: curr_status
  type(plyblock_layup), allocatable, optional, intent(out) :: layup(:)
  type(fBrickPly_elem), allocatable, optional, intent(out) :: plyblks(:)
  type(fCoh8Delam_elem),allocatable, optional, intent(out) :: interfs(:)

  if(present(curr_status)) curr_status = elem%curr_status

  if(present(layup)) then
      if(allocated(elem%layup)) then
          allocate(layup(size(elem%layup)))
          layup = elem%layup
      end if
  end if

  if(present(plyblks)) then
      if(allocated(elem%plyblks)) then
          allocate(plyblks(size(elem%plyblks)))
          plyblks = elem%plyblks
      end if
  end if

  if(present(interfs)) then
      if(allocated(elem%interfs)) then
          allocate(interfs(size(elem%interfs)))
          interfs = elem%interfs
      end if
  end if

end subroutine extract_fBrickLam_elem



pure subroutine set_fBrickLam_elem (elem, NPLYBLKS, node_connec, layup)
use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE
use fBrickPly_elem_module,  only : set
use fCoh8Delam_elem_module, only : set

  type(fBrickLam_elem),  intent(inout) :: elem
  integer,  intent(in)  :: NPLYBLKS
  integer,  intent(in)  :: node_connec( NPLYBLKS   * NNODE_PLYBLK + &
                                     & (NPLYBLKS-1)* NNDIN_INTERF )
  type(plyblock_layup), intent(in) :: layup(NPLYBLKS)

  integer                  :: istat
  character(len=MSGLENGTH) :: emsg, msgloc
  integer                  :: nndtotal, jstart, jend
  integer                  :: i, j

  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ' set, fBrickLam element module'
  nndtotal = 0
  jstart = 0
  jend   = 0
  i = 0; j = 0

  ! check input validity: leave it for preproc

  ! set local copy el first, copy it to elem in the end

  ! total no. of nodes in this element
  nndtotal = NPLYBLKS * NNODE_PLYBLK + (NPLYBLKS-1)* NNDIN_INTERF

  ! set NPLYBLKS
  elem%NPLYBLKS = NPLYBLKS

  ! set node connec
  allocate(elem%node_connec(nndtotal))
  elem%node_connec = node_connec

  ! set layup
  allocate(elem%layup(NPLYBLKS))
  elem%layup = layup

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !***** allocate nodes and edges to plyblks and set plyblks elements *****
  allocate(elem%plyblks(NPLYBLKS))
  allocate(elem%plyblks_nodes(NPLYBLKS))
  allocate(elem%plyblks_edges(NPLYBLKS))
  do i = 1, NPLYBLKS
    allocate(elem%plyblks_nodes(i)%array(NNODE_PLYBLK))
    allocate(elem%plyblks_edges(i)%array(NEDGE))
    elem%plyblks_nodes(i)%array = 0
    elem%plyblks_edges(i)%array = 0
    elem%plyblks_nodes(i)%array = [(j, j = (i-1)*NNODE_PLYBLK+1, i*NNODE_PLYBLK)]
    elem%plyblks_edges(i)%array = [(j, j = (i-1)*NEDGE       +1, i*NEDGE       )]
    ! set this plyblk elem
    call set (elem%plyblks(i), ply_angle=layup(i)%angle,    &
    & node_connec=node_connec(elem%plyblks_nodes(i)%array), &
    & istat=istat, emsg=emsg)
    if (istat == STAT_FAILURE) then
      emsg = emsg//trim(msgloc)
      return
    end if
  end do
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! if more than one ply, then interf subelem needs to be allocated
  if (NPLYBLKS > 1) then
  
  !***** allocate nodes to interfs and set interfs elements *****
  allocate(elem%interfs(NPLYBLKS-1))
  allocate(elem%interfs_nodes(NPLYBLKS-1))
  allocate(elem%interfs_edges(NPLYBLKS-1))
  do i = 1, NPLYBLKS-1
      ! allocate array dimensions
      allocate(elem%interfs_nodes(i)%array(NNODE_INTERF))
      allocate(elem%interfs_edges(i)%array(NEDGE))
      elem%interfs_nodes(i)%array = 0
      elem%interfs_edges(i)%array = 0

      ! 1st half of interface real nodes comes from
      ! bottom plyblk elem top surface (nodes 5-8)
        elem%interfs_nodes(i)%array(        1 : NNDRL/2) =                      &
      & elem%plyblks_nodes(i)%array(NNDRL/2+1 : NNDRL  )
      ! 2nd half of interface real nodes comes from
      ! top plyblk elem bottom surface (nodes 1-4)
        elem%interfs_nodes(i  )%array(NNDRL/2+1 : NNDRL  ) =                    &
      & elem%plyblks_nodes(i+1)%array(        1 : NNDRL/2)
      ! 1st half of interface flo nodes come from
      ! bottm plyblk elem top surface (nodes 17-24)
        elem%interfs_nodes(i)%array(NNDRL+        1 : NNDRL+NNDFL/2) =          &
      & elem%plyblks_nodes(i)%array(NNDRL+NNDFL/2+1 : NNDRL+NNDFL  )
      ! 2nd half of interface flo nodes come from
      ! top plyblk elem bottom surface (nodes 9-16)
        elem%interfs_nodes(i  )%array(NNDRL+NNDFL/2+1 : NNDRL+NNDFL  ) =        &
      & elem%plyblks_nodes(i+1)%array(NNDRL+        1 : NNDRL+NNDFL/2)
      ! internal nodes of interface element fDelam8
        jstart = NPLYBLKS*NNODE_PLYBLK + (i-1)* NNDIN_INTERF+1
        jend   = NPLYBLKS*NNODE_PLYBLK +  i   * NNDIN_INTERF
        elem%interfs_nodes(i)%array(NNDRL+NNDFL+1 : NNDRL+NNDFL+NNDIN_INTERF) = &
      & [(j, j=jstart,jend)]

      ! edges of this interface element
      ! 1st half of interface edges comes from
      ! bottom plyblk elem top edges (edges 5-8)
        elem%interfs_edges(i)%array(        1 : NEDGE/2) =                      &
      & elem%plyblks_edges(i)%array(NEDGE/2+1 : NEDGE  )
      ! 2nd half of interface edges comes from
      ! top plyblk elem bottom edges (edges 1-4)
        elem%interfs_edges(i  )%array(NEDGE/2+1 : NEDGE  ) =                    &
      & elem%plyblks_edges(i+1)%array(        1 : NEDGE/2)

      ! set this fDelam8 element
      call set (elem%interfs(i), node_connec(elem%interfs_nodes(i)%array),        &
      & istat, emsg)
      if (istat == STAT_FAILURE) then
        emsg = emsg//trim(msgloc)
        return
      end if
  end do
  
  end if ! NPLYBLK > 1
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!


end subroutine set_fBrickLam_elem



pure subroutine integrate_fBrickLam_elem (elem, nodes, edge_status,        &
& plylam_mat, plycoh_mat, interf_mat, K_matrix, F_vector, istat, emsg)

use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, NDIM, &
                      & ZERO, NDIM, MATRIX_CRACK_ELEM
use fnode_module,             only : fnode
use lamina_material_module,   only : lamina_material, lamina_scaled_Gfc
use cohesive_material_module, only : cohesive_material
use fBrickPly_elem_module,    only : extract, integrate
use fCoh8Delam_elem_module,   only : extract, update, integrate
use global_toolkit_module,    only : assembleKF

  type(fBrickLam_elem),     intent(inout) :: elem
  type(fnode),              intent(inout) :: nodes(size(elem%node_connec))
  integer,                  intent(inout) :: edge_status(elem%NPLYBLKS*NEDGE)
  type(lamina_material),    intent(in)    :: plylam_mat
  type(cohesive_material),  intent(in)    :: plycoh_mat
  type(cohesive_material),  intent(in)    :: interf_mat
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg


  ! local variables
  character(len=MSGLENGTH) :: msgloc
  integer                  :: nplyblks, ninterfs, ndof
  type(fnode)              :: plyblknds(NNODE_PLYBLK), interfnds(NNODE_INTERF)
  integer                  :: subegstatus(NEDGE)
  type(lamina_material)    :: plyblklam_mat
  real(DP)                 :: theta1, theta2
  real(DP), allocatable    :: Ki(:,:), Fi(:)
  logical                  :: topsurf_set, botsurf_set
  integer                  :: plyblk_status, plyblk_egstatus_lcl(NEDGE/2)
  ! loop counters
  integer :: i

  ! initialize intent out and local variables
  istat       = STAT_SUCCESS
  emsg        = ''
  msgloc      = ' integrate, fBrickLam element module'
  nplyblks    = 0
  ninterfs    = 0
  ndof        = 0
  subegstatus = 0
  theta1      = ZERO
  theta2      = ZERO
  topsurf_set = .false.
  botsurf_set = .false.
  plyblk_status       = 0
  plyblk_egstatus_lcl = 0
  i=0

  ! check if the elem has been set by checking layup .
  ! note that in set subroutine all attributes  must be allocated together,
  ! so checking one of them (here layup) would suffice.
  if(.not.allocated(elem%layup)) then
    istat = STAT_FAILURE
    emsg  = 'fBrickLam element must be set before integration,'//trim(msgloc)
    return
  end if


  ! extract no. plyblock and no. interfaces from layup, and calculate ndof
  nplyblks = elem%NPLYBLKS
  ninterfs = nplyblks - 1
  ndof     = NDIM * size(elem%node_connec)

  ! initialize K & F
  allocate(K_matrix(ndof,ndof),F_vector(ndof))
  K_matrix = ZERO
  F_vector = ZERO

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !***** integrate plyblock elements and assemble into global matrix *****
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  loop_plyblks: do i = 1, nplyblks
      ! extract nodes and edge status of this plyblk for integration
      plyblknds     = nodes(elem%plyblks_nodes(i)%array)
      subegstatus   = edge_status(elem%plyblks_edges(i)%array)
      ! increase fibre toughness w.r.t no. plies in this plyblock
      plyblklam_mat = lamina_scaled_Gfc(plylam_mat, elem%layup(i)%nplies)
      ! integrate this plyblk elem and update its nodes and edge status
      call integrate (elem%plyblks(i), plyblknds, subegstatus, plyblklam_mat, &
      & plycoh_mat, Ki, Fi, istat, emsg)
      if (istat == STAT_FAILURE) exit loop_plyblks
      ! assemble this elem's K and F
      call assembleKF (K_matrix, F_vector, Ki, Fi, elem%plyblks_nodes(i)%array, &
      & NDIM, istat, emsg)
      if (istat == STAT_FAILURE) exit loop_plyblks
      ! update nodes and edge status
      nodes(elem%plyblks_nodes(i)%array)       = plyblknds
      edge_status(elem%plyblks_edges(i)%array) = subegstatus
  end do loop_plyblks
  if (istat == STAT_FAILURE) then
    emsg = emsg//trim(msgloc)
    call zeroKF(K_matrix, F_vector)
    call clean_up(Ki, Fi)
    return
  end if


  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !***** integrate interface elements and assemble into global matrix *****
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  if (ninterfs > 0) then
  
  loop_interfs: do i = 1, ninterfs
      !----- update edge status -----
      ! extract status of this interface
      call extract (elem%interfs(i), top_subelem_set=topsurf_set, &
      &                              bot_subelem_set=botsurf_set )
      ! if top surf of this interface has not yet partitioned
      if (.not. topsurf_set) then
        ! check if the top plyblk elem has reached final partition
        call extract (elem%plyblks(i+1), curr_status = plyblk_status)
        ! if it has reached final partition, extract its edge status array
        ! and update to this interface elem
        if (plyblk_status == MATRIX_CRACK_ELEM) then
          call extract (elem%plyblks(i+1), edge_status_lcl=plyblk_egstatus_lcl)
          call update  (elem%interfs(i  ), ply_edge_status=plyblk_egstatus_lcl, &
          & top_or_bottom='top', istat=istat, emsg=emsg)
          if (istat == STAT_FAILURE) exit loop_interfs
        end if
      end if
      ! if bot surf of this interface has not yet partitioned
      if (.not. botsurf_set) then
        ! check if the top plyblk elem has reached final partition
        call extract (elem%plyblks(i), curr_status = plyblk_status)
        ! if it has reached final partition, extract its edge status array
        ! and update to this interface elem
        if (plyblk_status == MATRIX_CRACK_ELEM) then
          call extract (elem%plyblks(i), edge_status_lcl=plyblk_egstatus_lcl)
          call update  (elem%interfs(i), ply_edge_status=plyblk_egstatus_lcl, &
          & top_or_bottom='bottom', istat=istat, emsg=emsg)
          if (istat == STAT_FAILURE) exit loop_interfs
        end if
      end if
      !----- integration and assembly -----
      ! extract nodes of this interf for integration
      interfnds = nodes(elem%interfs_nodes(i)%array)
      ! extract top and bot ply angles of this interf
      theta1 = elem%layup(i)%angle
      theta2 = elem%layup(i+1)%angle
      ! integrate this interf elem and update its nodes
      call integrate (elem%interfs(i), interfnds, interf_mat, theta1, theta2,   &
      & Ki, Fi, istat, emsg)
      if (istat == STAT_FAILURE) exit loop_interfs
      ! assemble this elem's K and F
      call assembleKF (K_matrix, F_vector, Ki, Fi, elem%interfs_nodes(i)%array, &
      & NDIM, istat, emsg)
      if (istat == STAT_FAILURE) exit loop_interfs
      ! update nodes
      nodes(elem%interfs_nodes(i)%array) = interfnds
  end do loop_interfs
  if (istat == STAT_FAILURE) then
    emsg = emsg//trim(msgloc)
    call zeroKF(K_matrix, F_vector)
    call clean_up(Ki, Fi)
    return
  end if
  
  end if ! ninterf > 0

  call clean_up(Ki, Fi)

  return


  contains


  pure subroutine zeroKF (K, F)
    real(DP), intent(inout) :: K(:,:), F(:)
    K = ZERO
    F = ZERO
  end subroutine zeroKF

  pure subroutine clean_up (Ki, Fi)
    real(DP), allocatable, intent(inout) :: Ki(:,:), Fi(:)
    if(allocated(Ki)) deallocate(Ki)
    if(allocated(Fi)) deallocate(Fi)
  end subroutine clean_up


end subroutine integrate_fBrickLam_elem




end module fBrickLam_elem_module
