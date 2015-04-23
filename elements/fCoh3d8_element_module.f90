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
!  topological definition of bot sub element (reverse view), type fCoh3d8_subelem:
!
!  1_____9___BE1__10___2
!  |\                  |\
!  | \16               | \11
!  |  \BE4             |  \BE2       bot subelem top edges in this elem bot edges  
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
integer, parameter :: TOP_SUBELEM_NODES(20) = [1,2,3,4,5,6,7,8,  &
                   &  17,18,19,20,21,22,23,24,   25,26,27,28]

! NODAL CONNEC OF BOT SUB ELEM: REAL NODES, FLOATING NODES AND INTERNAL NODES                   
integer, parameter :: BOT_SUBELEM_NODES(20) = [8,7,6,5,4,3,2,1,  &
                   &  14,13,12,11,10, 9,16,15,   31,30,29,32]

type, public :: fCoh3d8_element
    private
    
    integer :: node_connec(NNODE) = 0
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
    emsg  = 'node connec indices must be >=1, set, &
    &fCoh3d8_subelem_module'
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



pure subroutine update_fCoh3d8_element (elem, surf_edge_status, top_or_bottom, &
& istat, emsg)
! Purpose:
! this subroutine is used to update the edge_status array of the element
! passed in from the adjacent ply elements, and then alloc and set the subelem
use parameter_module, only : MSGLENGTH, STAT_FAILURE, STAT_SUCCESS,&
                      & INTACT, TRANSITION_EDGE, REFINEMENT_EDGE,  &
                      & CRACK_TIP_EDGE, WEAK_CRACK_EDGE,           &
                      & COH_CRACK_EDGE, STRONG_CRACK_EDGE
                      
  type(fCoh3d8_element),    intent(inout) :: elem
  integer,                  intent(in)    :: surf_edge_status(NEDGE/2)
  character(len=*),         intent(in)    :: top_or_bottom
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  ! local copy of elem
  type(fCoh3d8_element)    :: el
  character(len=MSGLENGTH) :: msgloc
  integer                  :: n_crackedges

  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ' update, fCoh3d8_element module'
  n_crackedges = 0

  ! check edge status, see if there's any unexpected edge status value
  if ( any( .not. ( surf_edge_status == INTACT          .or.         &
  &                 surf_edge_status == TRANSITION_EDGE .or.         &
  &                 surf_edge_status == REFINEMENT_EDGE .or.         &
  &                 surf_edge_status == CRACK_TIP_EDGE  .or.         &
  &                 surf_edge_status == WEAK_CRACK_EDGE .or.         &
  &                 surf_edge_status == COH_CRACK_EDGE  .or.         &
  &                 surf_edge_status == STRONG_CRACK_EDGE )  )  ) then
    istat = STAT_FAILURE
    emsg  = 'surf edge status value is NOT recognized,'//trim(msgloc)
    return
  end if

  ! check the no. of broken edges; only accepts TWO cracked edges, as this is 
  ! the final partition from the ply element
  n_crackedges = count (surf_edge_status >= COH_CRACK_EDGE)
  if (n_crackedges /= 2) then
    istat = STAT_FAILURE
    emsg  = 'no. of cracked edges must be TWO,'//trim(msgloc)
    return
  end if
  
  ! copy elem to its local copy
  el = elem

  ! update to elem component if checkings are passed
  select case (trim(adjustl(top_or_bottom))
    
    case ('top')
        ! check, flag error if top elem has already been set
        if (allocated(el%top_subelem)) then 
          istat = STAT_FAILURE
          emsg  = 'top sub elem has already been defined,'//trim(msgloc)
          return
        end if
        ! proceed if no error encountered
        
        ! deallocate intact elem
        if (allocated(el%intact_elem)) deallocate(el%intact_elem)
        
        ! allocate top sub elem
        allocate(el%top_subelem)
        
        ! set top sub elem; note that the top sub elem's top edges are just
        ! this elem's top edges, without any change of order 
        call set (el%top_subelem, node_connec=el%node_connec(TOP_SUBELEM_NODES),&
        & top_edge_status=surf_edge_status, istat=istat, emsg=emsg)
        if (istat == STAT_FAILURE) then
          emsg = emsg//trim(msgloc)
          return
        end if
    
    case ('bottom')
        ! check, flag error if bot elem has already been set
        if (allocated(el%bot_subelem)) then 
          istat = STAT_FAILURE
          emsg  = 'bot sub elem has already been defined,'//trim(msgloc)
          return
        end if
        ! proceed if no error encountered
        
        ! deallocate intact elem
        if (allocated(el%intact_elem)) deallocate(el%intact_elem)
        
        ! allocate top sub elem
        allocate(el%bot_subelem)
        
        ! set top sub elem; note that the bot sub elem's top edges are
        ! this elem's bot edges but with permutated order (see illustration at top
        ! of this module, in comment)
        call set (el%bot_subelem, node_connec=el%node_connec(BOT_SUBELEM_NODES),&
        & top_edge_status=surf_edge_status([3, 2, 1, 4]), istat=istat, &
        & emsg=emsg)
        if (istat == STAT_FAILURE) then
          emsg = emsg//trim(msgloc)
          return
        end if
    
    case default
        istat = STAT_FAILURE
        emsg  = "unsupported top_or_bottom value, input either 'top' or &
        &'bottom',"//trim(msgloc)
        return
  
  end select
  
  ! copy definition to input arg. before successful return
  elem = el
    
end subroutine update_fCoh3d8_element



pure subroutine extract_fCoh3d8_element (elem, intact_elem, top_subelem, bot_subelem)
use coh3d8_element_module,  only : coh3d8_element
use fCoh3d8_subelem_module, only : fCoh3d8_subelem

    type(fCoh3d8_element),                        intent(in)   :: elem
    type(coh3d8_element),  allocatable, optional, intent(out)  :: intact_elem
    type(fCoh3d8_subelem), allocatable, optional, intent(out)  :: top_subelem
    type(fCoh3d8_subelem), allocatable, optional, intent(out)  :: bot_subelem

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



pure subroutine integrate_fCoh3d8_element (elem, K_matrix, F_vector)

    type(fCoh3d8_element),intent(inout)       :: elem 
    real(kind=dp),allocatable,intent(out)   :: K_matrix(:,:), F_vector(:)


    ! local variables
    type(int_alloc_array), allocatable  :: mainglbcnc(:)
    
    integer :: i,j,l, elstat, mainelstat

    logical :: nofailure
    
      

    ! initialize K & F
    allocate(K_matrix(ndof,ndof),F_vector(ndof))
    K_matrix=zero; F_vector=zero
    
    ! initialize local variables
    i=0; j=0; l=0
    elstat=0; mainelstat=0
    
    nofailure=.false.

    ! extract current status value
    elstat=elem%pstat  
    
    ! if elem is intact
    if(elstat==intact) then
    
        ! assign 1 coh3d8 elem as the main elem before failure, if not yet done
        if(.not.allocated(elem%intelem)) then 
            allocate(elem%intelem(1))
            allocate(elem%intelem_lcl_connec(1))
            allocate(elem%intelem_lcl_connec(1)%array(nndrl))   ! coh3d8 elem
            allocate(mainglbcnc(1))
            allocate(mainglbcnc(1)%array(nndrl))
            ! main elm 1 connec
            elem%intelem_lcl_connec(1)%array=[(i, i=1,nndrl)]
            mainglbcnc(1)%array(:)=elem%node_connec(elem%intelem_lcl_connec(1)%array(:))
            ! create sub elements
            call set(elem%intelem(1),key=0,connec=mainglbcnc(1)%array,matkey=elem%matkey)   
        end if   
        
        
        ! check if elem has started to fail; if so, no more edge status partitioning later
        call extract(elem%intelem(1),pstat=mainelstat)   
        
        if(mainelstat>intact) then
        ! if elem has reached failure onset, then update curr status
            elstat=elfail1
            elem%pstat=elstat
        end if

        call partition(elem)
        
        ! no damage/failure modelling at the iteration of new partition
        if(elem%pstat/=elstat) nofailure=.true.
        
        call integrate_assemble(elem,K_matrix,F_vector,nofailure)

    else if(elstat==elfail1) then
    ! if main elem coh3d8 has started to fail
    
        call partition(elem)
        
        ! no damage/failure modelling at the iteration of new partition
        if(elem%pstat/=elstat) nofailure=.true.
        
        call integrate_assemble(elem,K_matrix,F_vector,nofailure)   
    
    else if(elstat==elfail2) then
    ! if elem has already been partitioned into 2 subxcoh elems,
    ! update their lcl_ID_crack_edges arrays and elem curr status 
        
        call update_edgestatus(elem)
        
        ! no damage/failure modelling at the iteration of new partition
        if(elem%pstat/=elstat) nofailure=.true.
        
        call integrate_assemble(elem,K_matrix,F_vector,nofailure) 
    
    else if(elstat==elfail3) then
    
        call integrate_assemble(elem,K_matrix,F_vector,nofailure) 
    
    else
  write(msg_file,*)'unsupported elstat value in xcoh elem module'
  call exit_function
    end if      
        

    !---------------------------------------------------------------------!
    !               deallocate local arrays 
    !---------------------------------------------------------------------!

    if(allocated(mainglbcnc)) deallocate(mainglbcnc)


end subroutine integrate_fCoh3d8_element

 
 
pure subroutine partition (elem)
! check edge status, and partition into 2 subxcoh elems if any edge status is not intact
    
    type(fCoh3d8_element), intent(inout) :: elem
    
    type(int_alloc_array), allocatable  :: subglbcnc(:),subedge_connec(:)
    
    integer :: i
    
    i=0
    
    
    if(elem%pstat==intact .or. elem%pstat==elfail1) then
    
        if(maxval(elem%lcl_ID_crack_edges)>0) then

            elem%pstat=elfail2
            
            ! deallocate original elem
            deallocate(elem%intelem)
            deallocate(elem%intelem_lcl_connec)
            
            ! allocate two subxcoh elems
            allocate(elem%subelem(2))
            allocate(elem%subelem_lcl_connec(2))
            allocate(subglbcnc(2))
            allocate(subedge_connec(2))
            
            do i=1, 2
                allocate(elem%subelem_lcl_connec(i)%array(16))
                allocate(subglbcnc(i)%array(16))
                allocate(subedge_connec(i)%array(4))
                elem%subelem_lcl_connec(i)%array=0
                subglbcnc(i)%array=0
                subedge_connec(i)%array=0
            end do
            
            ! local connec of two subxcoh elems
            elem%subelem_lcl_connec(1)%array=[1,2,3,4,5,6,7,8,17,18,19,20,21,22,23,24]
            elem%subelem_lcl_connec(2)%array=[6,5,8,7,2,1,4,3,10,9,16,15,14,13,12,11]
            
            ! glb connec of two subxcoh elems
            subglbcnc(1)%array(:)=elem%node_connec(elem%subelem_lcl_connec(1)%array(:))
            subglbcnc(2)%array(:)=elem%node_connec(elem%subelem_lcl_connec(2)%array(:))
            
            ! glb edge cnc of two subxcoh elems
            subedge_connec(1)%array=elem%edge_connec([5,6,7,8])
            subedge_connec(2)%array=elem%edge_connec([1,4,3,2])
            
            ! set two subxcoh elems
            call set(elem%subelem(1),key=0,matkey=elem%matkey,&
            & node_connec=subglbcnc(1)%array,edge_connec=subedge_connec(1)%array)
            call set(elem%subelem(2),key=0,matkey=elem%matkey,&
            & node_connec=subglbcnc(2)%array,edge_connec=subedge_connec(2)%array)
            
            ! update edge status array of two sub elems
            call update_edgestatus(elem)
            
        end if
        
    else
        write(msg_file,*) 'unsupported elstat in xcoh partition!'
        call exit_function
    end if
    
    if(allocated(subglbcnc)) deallocate(subglbcnc)
    if(allocated(subedge_connec)) deallocate(subedge_connec)

end subroutine partition
 


pure subroutine update_edgestatus (elem)

    type(fCoh3d8_element), intent(inout) :: elem
    
    
    integer :: subelstat1, subelstat2, nfe1, nfe2
    integer, allocatable :: lcl_ID_crack_edges1(:), lcl_ID_crack_edges2(:)
    integer :: i
    
    subelstat1=0; subelstat2=0; nfe1=0; nfe2=0
    i=0
    
    
    if(elem%pstat==elfail2) then

        call extract(elem%subelem(1),pstat=subelstat1)
        call extract(elem%subelem(2),pstat=subelstat2)
        
        
        if(subelstat1==elfail2 .and. subelstat2==elfail2) then
        ! if both subxcoh elems have reached final partition state, 
        ! then update elstat to elfail3 and no lcl_ID_crack_edges update is needed;
            elem%pstat=elfail3
            
        else
        ! update lcl_ID_crack_edges arrays of two subxcoh elems
            
            allocate(lcl_ID_crack_edges1(4)); lcl_ID_crack_edges1=0
            allocate(lcl_ID_crack_edges2(4)); lcl_ID_crack_edges2=0
            
            nfe1=0
            nfe2=0
            do i=1, nedge
                select case (elem%lcl_ID_crack_edges(i))
                    case(1) ! edge 1 here is edge 1 of subxcoh2
                        nfe2=nfe2+1
                        lcl_ID_crack_edges2(nfe2)=1
                    case(2) ! edge 2 here is edge 4 of subxcoh2
                        nfe2=nfe2+1
                        lcl_ID_crack_edges2(nfe2)=4
                    case(3) ! edge 3 here is edge 3 of subxcoh2
                        nfe2=nfe2+1
                        lcl_ID_crack_edges2(nfe2)=3
                    case(4) ! edge 4 here is edge 2 of subxcoh2
                        nfe2=nfe2+1
                        lcl_ID_crack_edges2(nfe2)=2
                    case(5:8) ! top 4 edges are the 4 edges of subxcoh1
                        nfe1=nfe1+1
                        lcl_ID_crack_edges1(nfe1)=elem%lcl_ID_crack_edges(i)-4
                    case(0)
                    ! do nothing
                        continue
                    case default
                        write(msg_file,*)'sth wrong in xcoh integration subxcoh lcl_ID_crack_edges update!'
                        call exit_function            
                end select   
            end do
            
            
            if(subelstat1<elfail2) call update(elem%subelem(1),lcl_ID_crack_edges=lcl_ID_crack_edges1)
            if(subelstat2<elfail2) call update(elem%subelem(2),lcl_ID_crack_edges=lcl_ID_crack_edges2)  

        end if

    else
        write(msg_file,*)'unsupported elstat value in xcoh update edge status!'
        call exit_function
    end if
    
    if(allocated(lcl_ID_crack_edges1)) deallocate(lcl_ID_crack_edges1)
    if(allocated(lcl_ID_crack_edges2)) deallocate(lcl_ID_crack_edges2)

end subroutine update_edgestatus



pure subroutine integrate_assemble (elem, K_matrix, F_vector, nofailure) 
!---------------------------------------------------------------------!
!       integrate and assemble sub element system arrays
!---------------------------------------------------------------------!        
  ! - passed in variables   
  type(fCoh3d8_element), intent(inout)	:: elem
  real(kind=dp), 	intent(inout)			:: K_matrix(:,:), F_vector(:)
    logical, intent(in)                 :: nofailure
  ! - local variables
  real(kind=dp),	allocatable           	:: Ki(:,:), Fi(:)   ! sub_elem K matrix and F vector
  integer, 		allocatable 			:: dofcnc(:)
    integer :: i,j,l
    
    i=0;j=0;l=0
    
    ! empty K and F for reuse
    K_matrix=zero; F_vector=zero  
 
    ! integrate sub elements and assemble into global matrix
    ! if elem partition just changed in this iteration, no failure modelling at this iteration
    if(allocated(elem%intelem)) then
    
        call integrate(elem%intelem(1),Ki,Fi)
        if(allocated(dofcnc)) deallocate(dofcnc)
        allocate(dofcnc(size(Fi))); dofcnc=0
        do j=1, size(elem%intelem_lcl_connec(1)%array) ! no. of nodes in sub elem i
            do l=1, ndim
                ! dof indices of the jth node of sub elem i 
                dofcnc((j-1)*ndim+l)=(elem%intelem_lcl_connec(1)%array(j)-1)*ndim+l
            end do
        end do
        call assembleKF(K_matrix,F_vector,Ki,Fi,dofcnc)
        deallocate(Ki)
        deallocate(Fi)
        deallocate(dofcnc)
    
    else if(allocated(elem%subelem)) then
        
        do i=1, size(elem%subelem)
            call integrate(elem%subelem(i),Ki,Fi,nofailure)
            if(allocated(dofcnc)) deallocate(dofcnc)
            allocate(dofcnc(size(Fi))); dofcnc=0
            do j=1, size(elem%subelem_lcl_connec(i)%array) ! no. of nodes in sub elem i
                do l=1, ndim
                    ! dof indices of the jth node of sub elem i 
                    dofcnc((j-1)*ndim+l)=(elem%subelem_lcl_connec(i)%array(j)-1)*ndim+l
                end do
            end do
            ! each subxcoh elem contributes to half(in weight) of the system matrices
            call assembleKF(K_matrix,F_vector,half*Ki,half*Fi,dofcnc)
            deallocate(Ki)
            deallocate(Fi)
            deallocate(dofcnc)
        end do

    else
        write(msg_file,*)'elem not allocated in xcoh element module!'
        call exit_function
    end if
    
    if(allocated(Ki)) deallocate(Ki)
    if(allocated(Fi)) deallocate(Fi)
    if(allocated(dofcnc)) deallocate(dofcnc)
            
end subroutine integrate_assemble 




end module fCoh3d8_element_module
