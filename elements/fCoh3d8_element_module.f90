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
!  |   \24             |   \19      E5, E6, E7, E8
!  |    \              |    \
!  |     \5___17__(29)_|_18__\6
!  |______|____________|      |
! 4\   14 | (27)  13  3\      |
!   \     |             \     |
!  15\    |            12\    |
! (28)\   |           (26)\   |     Bottom edges (anti-clock wise from front):
!    16\  |              11\  |     E1, E2, E3, E4
!       \ |                 \ |
!        \|__________________\|
!         1    9  (25)  10    2
!
!  :::: bottom surface definition ::::
!  real     nodes on bot edges :  1,  2,  3,  4
!  floating nodes on bot edges :  9, 10, 11, 12, 13, 14, 15, 16
!  internal nodes on bot edges : 25, 26, 27, 28
!  nodes of bot E1 ::    <end nodes: 1, 2>  <fl. nodes:  9, 10>  <in. node: 25>
!  nodes of bot E2 ::    <end nodes: 2, 3>  <fl. nodes: 11, 12>  <in. node: 26>
!  nodes of bot E3 ::    <end nodes: 3, 4>  <fl. nodes: 13, 14>  <in. node: 27>
!  nodes of bot E4 ::    <end nodes: 4, 1>  <fl. nodes: 15, 16>  <in. node: 28>
!
!  :::: top surface definition ::::
!  real     nodes on top edges :  5,  6,  7,  8
!  floating nodes on top edges : 17, 18, 19, 20, 21, 22, 23, 24
!  internal nodes on top edges : 29, 30, 31, 32
!  nodes of top E1 ::    <end nodes: 5, 6>  <fl. nodes: 17, 18>  <in. node: 29>
!  nodes of top E2 ::    <end nodes: 6, 7>  <fl. nodes: 19, 20>  <in. node: 30>
!  nodes of top E3 ::    <end nodes: 7, 8>  <fl. nodes: 21, 22>  <in. node: 31>
!  nodes of top E4 ::    <end nodes: 8, 5>  <fl. nodes: 23, 24>  <in. node: 32>
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
!  8____22________21___7
!  |\                  |\
!  | \23               | \20
!  |  \                |  \         Top edges : E5, E6, E7, E8
!  |   \24             |   \19
!  |    \              |    \
!  |     \5___17_______|_18__\6
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
!  1_____9________10___2
!  |\                  |\
!  | \16               | \11
!  |  \                |  \         Top edges : E3, E2, E1, E4
!  |   \15             |   \12
!  |    \              |    \
!  |     \4___14_______|_13__\3
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
                   
! NODAL CONNEC OF INTACT ELEM: REAL NODES ONLY
integer, parameter :: INTACTELEM_NODES(8) = [1,2,3,4,5,6,7,8]

! NODAL CONNEC OF 2 SUB ELEMS: REAL NODES, FLOATING NODES AND INTERNAL NODES
integer, parameter :: TOPSUBELEM_NODES(20) = [1,2,3,4,5,6,7,8,  &
                   &  17,18,19,20,21,22,23,24,   25,26,27,28]
integer, parameter :: BOTSUBELEM_NODES(20) = [8,7,6,5,4,3,2,1,  &
                   &  14,13,12,11,10, 9,16,15,   31,30,29,32]

! EDGE CONNEC OF 2 SUB ELEMS:
integer, parameter :: TOPSUBELEM_TOPEDGES(4) = [5, 6, 7, 8]
integer, parameter :: BOTSUBELEM_TOPEDGES(4) = [3, 2, 1, 4]

type, public :: fCoh3d8_element
    private
    
    integer :: pstat                      = 0
    integer :: node_connec(NNODE)         = 0
    integer :: lcl_ID_crack_edges(NEDGE)  = 0
    type(coh3d8_element),  allocatable   :: intactElem
    type(fCoh3d8_subelem), allocatable   :: topSubElem
    type(fCoh3d8_subelem), allocatable   :: botSubElem
    
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



! this subroutine is used to set the connectivity and material lib index of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
pure subroutine set_fCoh3d8_element (elem, node_connec, edge_connec)

    type(fCoh3d8_element),  intent(inout)   :: elem
    integer,                intent(in)      :: node_connec(nnode)
    integer,                intent(in)      :: edge_connec(nedge)


    elem%node_connec=node_connec
    elem%edge_connec=edge_connec

end subroutine set_fCoh3d8_element



! this subroutine is used to update the lcl_ID_crack_edges array of the element
pure subroutine update_fCoh3d8_element (elem, lcl_ID_crack_edges)

    type(fCoh3d8_element),    intent(inout)    :: elem
    integer,                  intent(in)    :: lcl_ID_crack_edges(nedge)

    elem%lcl_ID_crack_edges=lcl_ID_crack_edges

end subroutine update_fCoh3d8_element



pure subroutine extract_fCoh3d8_element (elem, pstat, key, matkey, node_connec, edge_connec, &
& lcl_ID_crack_edges,intelem,subelem,intelem_lcl_connec,subelem_lcl_connec,sdv)

    type(fCoh3d8_element),                      intent(in)  :: elem
    integer,                        optional, intent(out) :: pstat
    integer,                        optional, intent(out) :: key
    integer,                        optional, intent(out) :: matkey
    integer,            allocatable,optional, intent(out) :: node_connec(:)
    integer,            allocatable,optional, intent(out) :: edge_connec(:)
    integer,            allocatable,optional, intent(out) :: lcl_ID_crack_edges(:)
    type(coh3d8_element),allocatable,optional,intent(out) :: intelem(:)
    type(fCoh3d8_subelem),allocatable,optional, intent(out) :: subelem(:)
    type(int_alloc_array),allocatable,optional,intent(out):: intelem_lcl_connec(:),subelem_lcl_connec(:)
    type(sdv_array),    allocatable,optional, intent(out) :: sdv(:)

    if(present(pstat)) pstat=elem%pstat
    if(present(key)) key=elem%key 
    if(present(matkey)) matkey=elem%matkey
    
    if(present(node_connec)) then 
        allocate(node_connec(nnode))
        node_connec=elem%node_connec
    end if
    
    if(present(edge_connec)) then 
        allocate(edge_connec(nedge))
        edge_connec=elem%edge_connec
    end if

    if(present(lcl_ID_crack_edges)) then 
        allocate(lcl_ID_crack_edges(nedge))
        lcl_ID_crack_edges=elem%lcl_ID_crack_edges
    end if
    
    if(present(intelem)) then
        if(allocated(elem%intelem)) then
            allocate(intelem(size(elem%intelem)))
            intelem=elem%intelem
        end if
    end if
    
    if(present(subelem)) then
        if(allocated(elem%subelem)) then
            allocate(subelem(size(elem%subelem)))
            subelem=elem%subelem
        end if
    end if
    
    
    if(present(intelem_lcl_connec)) then
        if(allocated(elem%intelem_lcl_connec)) then
            allocate(intelem_lcl_connec(size(elem%intelem_lcl_connec)))
            intelem_lcl_connec=elem%intelem_lcl_connec
        end if
    end if
    
    
    if(present(subelem_lcl_connec)) then
        if(allocated(elem%subelem_lcl_connec)) then
            allocate(subelem_lcl_connec(size(elem%subelem_lcl_connec)))
            subelem_lcl_connec=elem%subelem_lcl_connec
        end if
    end if
    
    if(present(sdv)) then        
        if(allocated(elem%sdv)) then
            allocate(sdv(size(elem%sdv)))
            sdv=elem%sdv
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
