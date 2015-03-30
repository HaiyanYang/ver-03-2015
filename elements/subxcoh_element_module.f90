module subxcoh_element_module
    use parameter_module
    use glb_clock_module
    use toolkit_module                  ! global tools for element integration
    use lib_edge_module                 ! global edge library
    use lib_node_module                 ! global node library
    use lib_mat_module                  ! global material library
    use coh3d6_element_module
    use coh3d8_element_module
    use sub3d_element_module
    use integration_point_module
  
  
    implicit none
    private

    integer,parameter :: ndim=3, nndrl=8, nedge=4, nndfl=2*nedge, nnode=nndrl+nndfl, ndof=ndim*nnode
    ! Topology: nodes on each edge; 4 nodes per edge, 1-2 are end nodes, 3-4 are fl. nodes; assigned in lcl node numbers
    integer,parameter :: topo(4,nedge)=reshape([5,6,9,10,6,7,11,12,7,8,13,14,8,5,15,16],[4,nedge])
                                              
    

    type, public :: subxcoh_element             ! breakable brick
        private
        
        integer :: curr_status=0        ! 0 means intact
        integer :: key=0 
        integer :: matkey=0
        
        integer :: nodecnc(nnode)=0     ! cnc to glb node arrays for accessing nodal variables (x, u, du, v, dof ...)
        integer :: edgecnc(nedge)=0     ! cnc to glb edge arrays for accessing edge variables (failure status)
        integer :: ifailedge(nedge)=0   ! indices of failed edges
        
        logical :: newpartition=.false. ! true when elem is changing partition, no failure should be considered then
        integer :: nstep=0, ninc=0      ! to store curr step and increment no.
        
        type(sub3d_element), allocatable :: subelem(:)
        type(int_alloc_array), allocatable :: subcnc(:)      ! sub_elem connec to parent elem nodes
        
        type(sdv_array), allocatable :: sdv(:)
        
    end type subxcoh_element
  
    interface empty
        module procedure empty_subxcoh_element
    end interface
  
    interface prepare
        module procedure prepare_subxcoh_element
    end interface
    
    interface update
        module procedure update_subxcoh_element
    end interface
    
    !~interface precrack
    !~    module procedure precrack_subxcoh_element
    !~end interface
    
    interface integrate
        module procedure integrate_subxcoh_element
    end interface
    
    interface extract
        module procedure extract_subxcoh_element
    end interface




    public :: empty,prepare,update,integrate,extract



    contains




    ! empty a breakable quadrilateral
    subroutine empty_subxcoh_element(elem)
  
        type(subxcoh_element),intent(out) :: elem
        
        elem%curr_status=0
        elem%key=0 
        elem%matkey=0
        
        elem%nodecnc=0
        elem%edgecnc=0
        elem%ifailedge=0
        
        elem%newpartition=.false.
        elem%nstep=0
        elem%ninc=0
        
        if(allocated(elem%subelem)) deallocate(elem%subelem)
        if(allocated(elem%subcnc))  deallocate(elem%subcnc)
        if(allocated(elem%sdv)) deallocate(elem%sdv)

    end subroutine empty_subxcoh_element
  
  
  
  
  
    ! this subroutine is used to prepare the connectivity and material lib index of the element
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine prepare_subxcoh_element(elem,key,matkey,nodecnc,edgecnc)
    
        type(subxcoh_element),    intent(inout)   :: elem
        integer,                intent(in)      :: key
        integer,                intent(in)      :: matkey
        integer,                intent(in)      :: nodecnc(nnode)
        integer,                intent(in)      :: edgecnc(nedge)

        elem%key=key 
        elem%matkey=matkey
        elem%nodecnc=nodecnc
        elem%edgecnc=edgecnc
    
    end subroutine prepare_subxcoh_element
 



    ! this subroutine is used to update the ifailedge array of the element
    subroutine update_subxcoh_element(elem,ifailedge)
    
        type(subxcoh_element),    intent(inout)   :: elem
        integer,                  intent(in)      :: ifailedge(nedge)

        elem%ifailedge=ifailedge
    
    end subroutine update_subxcoh_element


   
    
    subroutine extract_subxcoh_element(elem,curr_status,key,matkey,nodecnc,edgecnc, &
    & ifailedge,newpartition,nstep,ninc,subelem,subcnc,sdv)
    
        type(subxcoh_element),                      intent(in)  :: elem
        integer,                        optional, intent(out) :: curr_status
        integer,                        optional, intent(out) :: key
        integer,                        optional, intent(out) :: matkey
        integer,            allocatable,optional, intent(out) :: nodecnc(:)
        integer,            allocatable,optional, intent(out) :: edgecnc(:)
        integer,            allocatable,optional, intent(out) :: ifailedge(:)
        logical,                        optional, intent(out) :: newpartition
        integer,                        optional, intent(out) :: nstep, ninc
        type(sub3d_element),allocatable,optional, intent(out) :: subelem(:)
        type(int_alloc_array),allocatable,optional,intent(out):: subcnc(:)
        type(sdv_array),    allocatable,optional, intent(out) :: sdv(:)

        if(present(curr_status)) curr_status=elem%curr_status
        if(present(key)) key=elem%key 
        if(present(matkey)) matkey=elem%matkey
        
        if(present(nodecnc)) then 
            allocate(nodecnc(nnode))
            nodecnc=elem%nodecnc
        end if
        
        if(present(edgecnc)) then 
            allocate(edgecnc(nedge))
            edgecnc=elem%edgecnc
        end if

        if(present(ifailedge)) then 
            allocate(ifailedge(nedge))
            ifailedge=elem%ifailedge
        end if
        
        if(present(newpartition)) newpartition=elem%newpartition
        
        if(present(nstep)) nstep=elem%nstep
        
        if(present(ninc)) ninc=elem%ninc
        
        if(present(subelem)) then
            if(allocated(elem%subelem)) then
                allocate(subelem(size(elem%subelem)))
                subelem=elem%subelem
            end if
        end if
        
        if(present(subcnc)) then
            if(allocated(elem%subcnc)) then
                allocate(subcnc(size(elem%subcnc)))
                subcnc=elem%subcnc
            end if
        end if
        
        if(present(sdv)) then        
            if(allocated(elem%sdv)) then
                allocate(sdv(size(elem%sdv)))
                sdv=elem%sdv
            end if
        end if
    
    end subroutine extract_subxcoh_element






    subroutine integrate_subxcoh_element(elem, K_matrix, F_vector, nofailure)
    
        type(subxcoh_element),intent(inout)       :: elem 
        real(kind=dp),allocatable,intent(out)   :: K_matrix(:,:), F_vector(:)
        logical, optional, intent(in)           :: nofailure
    
    
        ! local variables
        type(int_alloc_array), allocatable  :: subglbcnc(:)     ! glb cnc of sub element, used when elem is intact
        
        integer :: i,j,l, elstat, subelstat
        
           
        ! - glb clock step and increment no. extracted from glb clock module
        integer :: curr_step, curr_inc
        logical :: last_converged               ! true if last iteration has converged: a new increment/step has started
        
        ! control parameter to prevent damage modelling if true
        logical :: nofail
        
    
        ! initialize K & F
        allocate(K_matrix(ndof,ndof),F_vector(ndof))
        K_matrix=zero; F_vector=zero
        
        ! initialize local variables
        i=0; j=0; l=0
        elstat=0; subelstat=0
        curr_step=0; curr_inc=0
        last_converged=.false.
        
        ! by default, damage modelling is allowed, unless specified
        nofail=.false.; if(present(nofailure)) nofail=nofailure

        ! assign 1 coh3d8 subelem if not yet done
        if(.not.allocated(elem%subelem)) then 
            allocate(elem%subelem(1))
            allocate(elem%subcnc(1))
            allocate(elem%subcnc(1)%array(nndrl))   ! coh3d8 elem
            allocate(subglbcnc(1))
            allocate(subglbcnc(1)%array(nndrl))
            ! sub elm 1 connec
            elem%subcnc(1)%array=[(i, i=1,nndrl)]
            subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
            ! create sub elements
            call prepare(elem%subelem(1),eltype='coh3d8', matkey=elem%matkey, &
            & glbcnc=subglbcnc(1)%array)
        end if
            

        ! - extract curr step and inc values from glb clock module
        call extract_glb_clock(kstep=curr_step,kinc=curr_inc)
        
        ! - check if last iteration has converged, and update the current step & increment no.
        if(elem%nstep.ne.curr_step .or. elem%ninc.ne.curr_inc) then
            last_converged=.true.
            elem%nstep = curr_step
            elem%ninc = curr_inc
        end if
        
        

        !---------------------------------------------------------------------!
        !       update elem partition using edge status variable
        !---------------------------------------------------------------------!

        ! extract current status value
        elstat=elem%curr_status     
        
        
        ! if elem is intact
        if(elstat==intact) then
        
            ! check if elem has started to fail
            call extract(elem%subelem(1),curr_status=subelstat)
            
            if(subelstat>intact) then 
            ! if elem reached delamination onset
                elstat=elfail1
                elem%curr_status=elstat
            end if
            
            call edge_status_partition(elem)
            call integrate_assemble(elem,K_matrix,F_vector,nofail) 
            
    	else if(elstat==elfail1) then
    	! elem already delaminated
            call edge_status_partition(elem)
    		call integrate_assemble(elem,K_matrix,F_vector,nofail) 
    	
    	else if(elstat==elfail2) then
    	! elem already edge partitioned; need to update mnodes, then integrate and assemble
    	
    		! update mnodes
    		call update_mnode(elem)
    		! integrate and assemble
    		call integrate_assemble(elem,K_matrix,F_vector,nofail) 
    	
    	else
         	write(msg_file,*)'unsupported elstat value in subxcoh elem module'
			call exit_function
        end if
        
        

        !---------------------------------------------------------------------!
        !               update elem sdv for output 
        !---------------------------------------------------------------------!
        if(last_converged) call update_sdv(elem)
  
        
        !---------------------------------------------------------------------!
        !               deallocate local arrays 
        !---------------------------------------------------------------------!
        
        if(allocated(subglbcnc)) deallocate(subglbcnc)
    
 
    end subroutine integrate_subxcoh_element
    
    



	subroutine integrate_assemble(elem,K_matrix,F_vector,nofail)
	!---------------------------------------------------------------------!
    !       integrate and assemble sub element system arrays
    !---------------------------------------------------------------------!     
    	! - passed in variables   
    	type(subxcoh_element), intent(inout)	:: elem
    	real(kind=dp), 	intent(inout)			:: K_matrix(:,:), F_vector(:)
        logical, intent(in)                     :: nofail
    	! - local variables
    	real(kind=dp),	allocatable           	:: Ki(:,:), Fi(:)   ! sub_elem K matrix and F vector
    	integer, 		allocatable 			:: dofcnc(:)
        integer :: i,j,l
        
        i=0;j=0;l=0
        
        ! empty K and F for reuse
        K_matrix=zero; F_vector=zero  
     
        ! integrate sub elements and assemble into global matrix
        do i=1, size(elem%subelem)
            call integrate(elem%subelem(i),Ki,Fi,nofail)
            if(allocated(dofcnc)) deallocate(dofcnc)
            allocate(dofcnc(size(Fi))); dofcnc=0
            do j=1, size(elem%subcnc(i)%array) ! no. of nodes in sub elem i
                do l=1, ndim
                    ! dof indices of the jth node of sub elem i 
                    dofcnc((j-1)*ndim+l)=(elem%subcnc(i)%array(j)-1)*ndim+l
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


    
    
    
    
    subroutine edge_status_partition(elem)

    
    ! passed-in variables
    type(subxcoh_element), intent(inout) :: elem


    ! extracted variables, from glb libraries
    integer :: edgstat(nedge)               ! status variable array of element edges
    type(real_alloc_array) :: coord(nnode)  ! nodal coord arrays to store the coords of elem nodes extracted from glb node lib
    
    
    
    ! local variable
    
    integer :: nfailedge        ! no. of failed edges in the element
    integer :: ifedg(nedge)     ! index of failed edges in the element
    integer :: elstat           ! local copy of elem curr status
    integer :: jbe1,jbe2,jbe3,jbe4, jnode ! indices of broken edges, and a variable to hold a glb node index
    integer :: iscross          ! indicator of line intersection; >0 if two lines intersect
    integer :: i, j, l, k       ! counters
      
    real(dp) :: xp1, yp1, xp2, yp2      ! (x,y) of point 1 and point 2 on the crack line
    real(dp) :: x1, y1, z1, x2, y2, z2  ! (x,y) of node 1 and node 2 of an element edge
    real(dp) :: xct, yct, zct           ! (x,y) of a crack tip on an edge, i.e., intersection of crack line & edge
    real(dp) :: detlc,a1,b1,a2,b2,xmid,ymid

    
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
    
    
        ! initialize local variables
        
        edgstat=0; 
        nfailedge=0; ifedg=0; elstat=0
        jbe1=0; jbe2=0; jbe3=0; jbe4=0; jnode=0
        iscross=0
        i=0; j=0; l=0; k=0
        
        xp1=zero; yp1=zero
        xp2=zero; yp2=zero

        x1=zero; y1=zero; z1=zero
        x2=zero; y2=zero; z2=zero
        
        xct=zero; yct=zero; zct=zero

        detlc=zero; a1=zero; b1=zero; a2=zero; b2=zero; xmid=zero; ymid=zero




!-----------------------------------------------------------------------!
!               EXTRACTION INTERFACE
!           extract variables from elem components and global libraries
!-----------------------------------------------------------------------!
        ! extract elem status variable
        elstat=elem%curr_status
        
        ! no need to update partition if elem is already in failed partition (2 broken edges)
        if(elstat==elfail2) return

        ! extract edge status variables from glb edge library
        edgstat(:)=lib_edge(elem%edgecnc(:))
        
        
        ! extract nodal coords from glb node library
        do i=1, nnode
            call extract(lib_node(elem%nodecnc(i)),x=coord(i)%array)
        end do

        ! extract elem failed edges' indices 
        ! this info is passed from adj. ply elem in the xcoh module
        ifedg=elem%ifailedge

!-----------------------------------------------------------------------!
!       procedure calculations (pure)
!-----------------------------------------------------------------------!
    
!       find the no. of broken edges
        do i=1,size(ifedg)
            if(ifedg(i)>0) nfailedge=nfailedge+1
        end do
    
!       update elstat w.r.t no. of failed edges and edge status
        if(nfailedge==0) then
        ! adj. ply elem remains intact, do nothing
            if(.not.(elstat==intact .or. elstat==elfail1)) then
                write(msg_file,*)'inconsistency btw elstat and nfailedge in subxcoh elem!'
                call exit_function
            end if
            
        else if(nfailedge==1) then 
        ! adj. ply elem is in transition partition; edge status should be egtrans
        
            if(edgstat(ifedg(1))==egtrans) then
              ! edge marks the refinement end and the trans elem start 
              ! elem is a trans elem, only this edge needs to be partitioned
                elstat=eltrans
                
            else
              ! ifailedge is not correct / edgstat not correct
              write(msg_file,*)'wrong edge status for nfailedge=1 in subxcoh'
              call exit_function
        
            endif
            
        else if(nfailedge==2) then
        ! adj. ply elem could be cracked, wake, tip, refinement elem
        ! only update the elstat, not the edge status, nor the crack tip coords      
            
            jbe1=ifedg(1)
            jbe2=ifedg(2)
            
            if(edgstat(jbe1)<=egref .and. edgstat(jbe2)<=egref) then
            ! refinement elem
                elstat=elref
            else if((edgstat(jbe1)<=egtip .and. edgstat(jbe2)==egtip).or. &
            & (edgstat(jbe2)<=egtip .and. edgstat(jbe1)==egtip)) then
            ! tip elem
                elstat=eltip
            else if((edgstat(jbe1)<cohcrack .and. edgstat(jbe2)>=cohcrack).or. &
            & (edgstat(jbe2)<cohcrack .and. edgstat(jbe1)>=cohcrack)) then
            ! wake elem, cohesive/stress-free crack
                elstat=elwake
            else if(edgstat(jbe1)>=cohcrack .and. edgstat(jbe2)>=cohcrack) then
            ! cracked elem, cohesive/stress-free crack
                elstat=elfail2
            else ! unknown combination
                write(msg_file,*)'unknown combination of 2 edge status in subxcoh!'
                call exit_function
            end if
            
        else
            write(msg_file,*)'unsupported nfailedge value for edge and el stat update in subxcoh edge stat partition!'
            call exit_function 
        end if     

!-----------------------------------------------------------------------!
!                   UPDATE INTERFACE
!               update global libraries
!-----------------------------------------------------------------------!
            
!       update element curr_status and sub-element cnc matrices

        !if(elstat>elem%curr_status) then
        if (elstat==elfail2) then
        ! only update elem curr status and partitions into sub elems when it reaches elfail2 partition
        ! otherwise, elem curr status remains intact

            elem%curr_status=elstat                    
            
            call update_subcnc(elem,edgstat,ifedg,nfailedge)
        
        end if





!       deallocate local dynamic arrays

    
    
    end subroutine edge_status_partition
    
  



       
    subroutine update_subcnc(elem,edgstat,ifedg,nfailedge)
    
    ! passed-in variables
    type(subxcoh_element),    intent(inout)   :: elem
    integer,                intent(in)      :: edgstat(:), ifedg(:), nfailedge



    ! local variables
    type(int_alloc_array), allocatable :: subglbcnc(:)  ! glb cnc of sub elements
    integer :: i, j, l                                  ! counters
    integer :: ibe, ibe1, ibe2                          ! indices of broken edges
    integer :: e1,e2,e3,e4                              ! edge indices, used for partitioning element
    integer :: nsub, nbulk                              ! no. of sub elements, and no. of bulk partitions
    integer :: jnode, jnode1, jnode2                    ! node index variables

    real(dp), allocatable   :: x1(:), x2(:), xc(:)      ! coords of broken edge end nodes (x1 and x2) and crack tip (xc)
    real(dp)                :: tratio, tratio1, tratio2 ! tratio=|xc-x1|/|x2-x1|
    type(xnode),allocatable :: mnode(:)                 ! material nodes of cohesive elem
    real(dp), allocatable   :: Tmatrix(:,:)             ! interpolation matrix btw bottom num nodes and mat nodes


!       initialize local variables

        i=0; j=0; l=0
        e1=0; e2=0; e3=0; e4=0
        ibe=0; ibe1=0; ibe2=0
        nsub=0; nbulk=0
        jnode=0; jnode1=0; jnode2=0
        tratio=zero; tratio1=zero; tratio2=zero
        
        

10      select case (nfailedge)
        case (0) !- no cracked edge, do nothing
            continue
            
            
        case (1) !- one edge cracked, trans partition
            ! find the index of the broken edge
            ibe=ifedg(1)

            ! ibe1 must be between 1 to 4
            if(ibe1<1 .or. ibe1>4) then
                write(msg_file,*) 'something wrong in subxcoh update subcnc case nfailedge=1'
                call exit_function
            end if

            ! verify its status variable value
            if(edgstat(ibe)/=egtrans) then
                write(msg_file,*)'transition partition only accepts edgstat=egtrans!'
                call exit_function
            end if
            
            ! allocate sub element arrays; in this case, 3 coh3d6 sub3d elements
            nsub=3
            if(allocated(elem%subelem)) deallocate(elem%subelem)
            if(allocated(elem%subcnc)) deallocate(elem%subcnc)
            if(allocated(subglbcnc)) deallocate(subglbcnc)
            allocate(elem%subelem(nsub))
            allocate(elem%subcnc(nsub))
            allocate(subglbcnc(nsub)) 
            ! these coh3d6 elems have 7 num nodes and 6 mat nodes
            do j=1, nsub
                allocate(elem%subcnc(j)%array(7))
                allocate(subglbcnc(j)%array(7))
                elem%subcnc(j)%array=0
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
                    write(msg_file,*)'wrong broken edge in tip partition, subcnc'
                    call exit_function
            end select
            
            ! find the smaller glb node on the broken edge
            if(elem%nodecnc(topo(3,ibe))<elem%nodecnc(topo(4,ibe))) then
                jnode=topo(3,ibe)
            else
                jnode=topo(4,ibe)
            end if
            
            ! find the relative position of crack tip on this edge
            call extract(lib_node(elem%nodecnc(topo(1,ibe))),x=x1)
            call extract(lib_node(elem%nodecnc(topo(2,ibe))),x=x2)
            call extract(lib_node(elem%nodecnc(jnode)),x=xc)
            tratio=distance(x1,xc)/distance(x1,x2)
            
            
            
            !*** sub elm 1 connec; 7 numerical nodes, 6 material nodes     
            
            elem%subcnc(1)%array(1)=topo(1,e1)-nndrl/2
            elem%subcnc(1)%array(2)=topo(2,e1)-nndrl/2
            elem%subcnc(1)%array(3)=topo(1,e3)-nndrl/2
            elem%subcnc(1)%array(4)=topo(2,e3)-nndrl/2
            elem%subcnc(1)%array(5)=topo(1,e1)
            elem%subcnc(1)%array(6)=topo(2,e1)
            elem%subcnc(1)%array(7)=jnode
            
            subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
            
            ! update the mnode array
            ! first two mat nodes of bottom surf are the same as num nodes
            mnode(1)=lib_node(subglbcnc(1)%array(1))
            mnode(2)=lib_node(subglbcnc(1)%array(2))
            mnode(3)=tratio*lib_node(subglbcnc(1)%array(4))+(one-tratio)*lib_node(subglbcnc(1)%array(1))
            mnode(4)=lib_node(subglbcnc(1)%array(5))
            mnode(5)=lib_node(subglbcnc(1)%array(6))
            mnode(6)=lib_node(subglbcnc(1)%array(7))
            
            Tmatrix=zero
            Tmatrix(1,1)=one
            Tmatrix(2,2)=one  
            Tmatrix(3,1)=one-tratio
            Tmatrix(3,4)=tratio
            
            call prepare(elem%subelem(1),eltype='coh3d6',matkey=elem%matkey,glbcnc=subglbcnc(1)%array &
            & ,Tmatrix=Tmatrix,mnode=mnode)
            
            
            
            
            !*** sub elm 2 connec; 7 numerical nodes; 6 material nodes
            
            elem%subcnc(2)%array(1)=topo(1,e2)-nndrl/2
            elem%subcnc(2)%array(2)=topo(2,e2)-nndrl/2
            elem%subcnc(2)%array(3)=topo(1,ibe)-nndrl/2
            elem%subcnc(2)%array(4)=topo(2,ibe)-nndrl/2
            elem%subcnc(2)%array(5)=topo(1,e2)
            elem%subcnc(2)%array(6)=topo(2,e2)
            elem%subcnc(2)%array(7)=jnode

            subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))
            
            ! update the mnode array
            ! first two mat nodes of bottom surf are the same as num nodes
            mnode(1)=lib_node(subglbcnc(2)%array(1))
            mnode(2)=lib_node(subglbcnc(2)%array(2))
            mnode(3)=tratio*lib_node(subglbcnc(2)%array(3))+(one-tratio)*lib_node(subglbcnc(2)%array(4))
            mnode(4)=lib_node(subglbcnc(2)%array(5))
            mnode(5)=lib_node(subglbcnc(2)%array(6))
            mnode(6)=lib_node(subglbcnc(2)%array(7))
            
            Tmatrix=zero
            Tmatrix(1,1)=one
            Tmatrix(2,2)=one  ! first two mat nodes of bottom surf are the same as num nodes
            Tmatrix(3,3)=tratio
            Tmatrix(3,4)=one-tratio
            
            call prepare(elem%subelem(2),eltype='coh3d6',matkey=elem%matkey,glbcnc=subglbcnc(2)%array &
            & ,Tmatrix=Tmatrix,mnode=mnode)
            
            
            
            
            !*** sub elm 3 connec; 7 numerical nodes, 6 material nodes
            
            elem%subcnc(3)%array(1)=topo(1,e3)-nndrl/2
            elem%subcnc(3)%array(2)=topo(2,e3)-nndrl/2
            elem%subcnc(3)%array(3)=topo(1,e1)-nndrl/2
            elem%subcnc(3)%array(4)=topo(2,e1)-nndrl/2 
            elem%subcnc(3)%array(5)=topo(1,e3)
            elem%subcnc(3)%array(6)=topo(2,e3)
            elem%subcnc(3)%array(7)=jnode  
            
            subglbcnc(3)%array(:)=elem%nodecnc(elem%subcnc(3)%array(:))
            
            ! update the mnode array
            ! first two mat nodes of bottom surf are the same as num nodes
            mnode(1)=lib_node(subglbcnc(3)%array(1))
            mnode(2)=lib_node(subglbcnc(3)%array(2))
            mnode(3)=tratio*lib_node(subglbcnc(3)%array(2))+(one-tratio)*lib_node(subglbcnc(3)%array(3))
            mnode(4)=lib_node(subglbcnc(3)%array(5))
            mnode(5)=lib_node(subglbcnc(3)%array(6))
            mnode(6)=lib_node(subglbcnc(3)%array(7))
            
            Tmatrix=zero
            Tmatrix(1,1)=one
            Tmatrix(2,2)=one  ! first two mat nodes of bottom surf are the same as num nodes
            Tmatrix(3,2)=tratio
            Tmatrix(3,3)=one-tratio
            
            call prepare(elem%subelem(3),eltype='coh3d6',matkey=elem%matkey,glbcnc=subglbcnc(3)%array &
            & ,Tmatrix=Tmatrix,mnode=mnode)
         

 
        case (2) !- two edges cracked
         
            ibe1=min(ifedg(1),ifedg(2))   ! local edge index of 1st broken edge
            ibe2=max(ifedg(1),ifedg(2)) 
            
            
            ! ibe1 must be between 1 to 3, and ibe2 between 2 to 4, with ibe2 > ibe1
            if(ibe1<1 .or. ibe1>3 .or. ibe2<2 .or. ibe2>4 .or. ibe2<=ibe1) then
                write(msg_file,*) 'something wrong in subxcoh update subcnc case nfailedge=2'
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
                            write(msg_file,*)'wrong 2nd broken edge in update subcnc subxcoh'
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
                            write(msg_file,*)'wrong 2nd broken edge in update subcnc subxcoh'
                            call exit_function
                    end select
                case(3)
                    if(ibe2==4) then
                        nbulk=4
                        e1=3; e2=4; e3=1; e4=2
                    else
                        write(msg_file,*)'wrong 2nd broken edge in update subcnc subxcoh'
                        call exit_function                    
                    end if    
                case default
                    write(msg_file,*)'wrong broken edge in update subcnc subxcoh'
                    call exit_function
            end select
            
            
            select case(nbulk)
                case(2)
                ! two quad subdomains

                    nsub=2  ! only two coh3d8 sub3d elems

                    if(allocated(elem%subelem)) deallocate(elem%subelem)
                    if(allocated(elem%subcnc)) deallocate(elem%subcnc)
                    if(allocated(subglbcnc)) deallocate(subglbcnc)
                    allocate(elem%subelem(nsub))
                    allocate(elem%subcnc(nsub))
                    allocate(subglbcnc(nsub))
                    
                    ! allocate 8 numerical nodes for sub elems
                    do j=1, nsub
                        allocate(elem%subcnc(j)%array(8))
                        allocate(subglbcnc(j)%array(8))
                        elem%subcnc(j)%array=0
                        subglbcnc(j)%array=0
                    end do
                    
                    ! allocate 8 material nodes for sub elems
                    if(allocated(mnode)) deallocate(mnode)
                    allocate(mnode(8))
                    
                    ! allocate Tmatrix for bottom surface (4 mat nodes interpolated by 4 num nodes)
                    if(allocated(Tmatrix)) deallocate(Tmatrix)
                    allocate(Tmatrix(4,4))
                    
                    ! find the smaller glb fl. node on the 1st broken edge
                    if(elem%nodecnc(topo(3,e1))<elem%nodecnc(topo(4,e1))) then
                        jnode1=topo(3,e1)
                    else
                        jnode1=topo(4,e1)
                    end if
                    
                    ! find the relative position of crack tip on this edge
                    call extract(lib_node(elem%nodecnc(topo(1,e1))),x=x1)
                    call extract(lib_node(elem%nodecnc(topo(2,e1))),x=x2)
                    call extract(lib_node(elem%nodecnc(jnode1)),x=xc)
                    tratio1=distance(x1,xc)/distance(x1,x2)
                    
                    ! find the smaller glb fl. node on the 2nd broken edge
                    if(elem%nodecnc(topo(3,e3))<elem%nodecnc(topo(4,e3))) then
                        jnode2=topo(3,e3)
                    else
                        jnode2=topo(4,e3)
                    end if
                    
                    ! find the relative position of crack tip on this edge
                    call extract(lib_node(elem%nodecnc(topo(1,e3))),x=x1)
                    call extract(lib_node(elem%nodecnc(topo(2,e3))),x=x2)
                    call extract(lib_node(elem%nodecnc(jnode2)),x=xc)
                    tratio2=distance(x1,xc)/distance(x1,x2)
                    
                    
                    !*** sub elm 1 connec
                    elem%subcnc(1)%array(1)=topo(1,e1)-nndrl/2
                    elem%subcnc(1)%array(2)=topo(2,e1)-nndrl/2
                    elem%subcnc(1)%array(3)=topo(1,e3)-nndrl/2
                    elem%subcnc(1)%array(4)=topo(2,e3)-nndrl/2
                    elem%subcnc(1)%array(5)=topo(1,e1)
                    elem%subcnc(1)%array(6)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(1)%array(6)=jnode1
                    elem%subcnc(1)%array(7)=topo(4,e3); if(edgstat(e3)<cohcrack) elem%subcnc(1)%array(7)=jnode2
                    elem%subcnc(1)%array(8)=topo(2,e3)
                    
                    subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
                    
                    ! update the mnode array
                    mnode(1)=lib_node(subglbcnc(1)%array(1))
                    mnode(2)=tratio1*lib_node(subglbcnc(1)%array(1))+(one-tratio1)*lib_node(subglbcnc(1)%array(2))
                    mnode(3)=tratio2*lib_node(subglbcnc(1)%array(3))+(one-tratio2)*lib_node(subglbcnc(1)%array(4))
                    mnode(4)=lib_node(subglbcnc(1)%array(4))
                    mnode(5)=lib_node(subglbcnc(1)%array(5))
                    mnode(6)=lib_node(subglbcnc(1)%array(6))
                    mnode(7)=lib_node(subglbcnc(1)%array(7))
                    mnode(8)=lib_node(subglbcnc(1)%array(8))
                    
                    Tmatrix=zero
                    Tmatrix(1,1)=one
                    Tmatrix(2,1)=tratio1
                    Tmatrix(2,2)=one-tratio1
                    Tmatrix(3,3)=tratio2
                    Tmatrix(3,4)=one-tratio2
                    Tmatrix(4,4)=one
                    
                    call prepare(elem%subelem(1),eltype='coh3d8',matkey=elem%matkey,glbcnc=subglbcnc(1)%array &
                    & ,Tmatrix=Tmatrix,mnode=mnode)
                    
                    
                    
                    
                    
                    !*** sub elm 2 connec
                    elem%subcnc(2)%array(1)=topo(1,e1)-nndrl/2
                    elem%subcnc(2)%array(2)=topo(2,e1)-nndrl/2
                    elem%subcnc(2)%array(3)=topo(1,e3)-nndrl/2
                    elem%subcnc(2)%array(4)=topo(2,e3)-nndrl/2
                    elem%subcnc(2)%array(5)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(2)%array(5)=jnode1
                    elem%subcnc(2)%array(6)=topo(2,e1)
                    elem%subcnc(2)%array(7)=topo(1,e3)
                    elem%subcnc(2)%array(8)=topo(3,e3); if(edgstat(e3)<cohcrack) elem%subcnc(2)%array(8)=jnode2

                    subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))
                    
                    ! update the mnode array
                    mnode(1)=tratio1*lib_node(subglbcnc(2)%array(1))+(one-tratio1)*lib_node(subglbcnc(2)%array(2))
                    mnode(2)=lib_node(subglbcnc(2)%array(2))
                    mnode(3)=lib_node(subglbcnc(2)%array(3))
                    mnode(4)=tratio2*lib_node(subglbcnc(2)%array(3))+(one-tratio2)*lib_node(subglbcnc(2)%array(4))
                    mnode(5)=lib_node(subglbcnc(2)%array(5))
                    mnode(6)=lib_node(subglbcnc(2)%array(6))
                    mnode(7)=lib_node(subglbcnc(2)%array(7))
                    mnode(8)=lib_node(subglbcnc(2)%array(8))
                    
                    Tmatrix=zero
                    Tmatrix(1,1)=tratio1
                    Tmatrix(1,2)=one-tratio1
                    Tmatrix(2,2)=one
                    Tmatrix(3,3)=one
                    Tmatrix(4,3)=tratio2
                    Tmatrix(4,4)=one-tratio2
                    
                    call prepare(elem%subelem(2),eltype='coh3d8',matkey=elem%matkey,glbcnc=subglbcnc(2)%array &
                    & ,Tmatrix=Tmatrix,mnode=mnode)
                    
    
                    
                    
                    
                    
                case(4)
                ! four triangular subdomains

                    nsub=4

                    if(allocated(elem%subelem)) deallocate(elem%subelem)
                    if(allocated(elem%subcnc)) deallocate(elem%subcnc)
                    if(allocated(subglbcnc)) deallocate(subglbcnc)
                    allocate(elem%subelem(nsub))
                    allocate(elem%subcnc(nsub))
                    allocate(subglbcnc(nsub))
                    do j=1, nsub   ! coh3d6 elem cnc
                        allocate(elem%subcnc(j)%array(7))
                        allocate(subglbcnc(j)%array(7))
                        elem%subcnc(j)%array=0
                        subglbcnc(j)%array=0
                    end do

                    ! allocate 6 material nodes for sub elems
                    if(allocated(mnode)) deallocate(mnode)
                    allocate(mnode(6))
                    
                    ! allocate Tmatrix for bottom surface (3 mat nodes interpolated by 4 num nodes)
                    if(allocated(Tmatrix)) deallocate(Tmatrix)
                    allocate(Tmatrix(3,4))

                    
                    ! find the smaller glb fl. node on the 1st broken edge
                    if(elem%nodecnc(topo(3,e1))<elem%nodecnc(topo(4,e1))) then
                        jnode1=topo(3,e1)
                    else
                        jnode1=topo(4,e1)
                    end if
                    
                    ! find the relative position of crack tip on this edge
                    call extract(lib_node(elem%nodecnc(topo(1,e1))),x=x1)
                    call extract(lib_node(elem%nodecnc(topo(2,e1))),x=x2)
                    call extract(lib_node(elem%nodecnc(jnode1)),x=xc)
                    tratio1=distance(x1,xc)/distance(x1,x2)
                    
                    ! find the smaller glb fl. node on the 2nd broken edge
                    if(elem%nodecnc(topo(3,e2))<elem%nodecnc(topo(4,e2))) then
                        jnode2=topo(3,e2)
                    else
                        jnode2=topo(4,e2)
                    end if
                    
                    ! find the relative position of crack tip on this edge
                    call extract(lib_node(elem%nodecnc(topo(1,e2))),x=x1)
                    call extract(lib_node(elem%nodecnc(topo(2,e2))),x=x2)
                    call extract(lib_node(elem%nodecnc(jnode2)),x=xc)
                    tratio2=distance(x1,xc)/distance(x1,x2)
                    
                    
                    
                    !*** sub elm 1 connec
                    elem%subcnc(1)%array(1)=topo(1,e1)-nndrl/2 ! upper surf nodes
                    elem%subcnc(1)%array(2)=topo(2,e1)-nndrl/2
                    elem%subcnc(1)%array(3)=topo(1,e3)-nndrl/2
                    elem%subcnc(1)%array(4)=topo(2,e3)-nndrl/2
                    elem%subcnc(1)%array(5)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(1)%array(5)=jnode1
                    elem%subcnc(1)%array(6)=topo(1,e2)
                    elem%subcnc(1)%array(7)=topo(3,e2); if(edgstat(e2)<cohcrack) elem%subcnc(1)%array(7)=jnode2
                    
                    subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
                    
                    mnode(1)=tratio1*lib_node(subglbcnc(1)%array(1))+(one-tratio1)*lib_node(subglbcnc(1)%array(2))
                    mnode(2)=lib_node(subglbcnc(1)%array(2))
                    mnode(3)=tratio2*lib_node(subglbcnc(1)%array(2))+(one-tratio2)*lib_node(subglbcnc(1)%array(3))
                    mnode(4)=lib_node(subglbcnc(1)%array(5))
                    mnode(5)=lib_node(subglbcnc(1)%array(6))
                    mnode(6)=lib_node(subglbcnc(1)%array(7))
                    
                    Tmatrix=zero
                    Tmatrix(1,1)=tratio1
                    Tmatrix(1,2)=one-tratio1
                    Tmatrix(2,2)=one
                    Tmatrix(3,2)=tratio2
                    Tmatrix(3,3)=one-tratio2
                    
                    call prepare(elem%subelem(1),eltype='coh3d6',matkey=elem%matkey,glbcnc=subglbcnc(1)%array &
                    & ,Tmatrix=Tmatrix,mnode=mnode)
 
                    
                    !*** sub elm 2 connec
                    elem%subcnc(2)%array(1)=topo(1,e2)-nndrl/2 ! upper surf nodes
                    elem%subcnc(2)%array(2)=topo(2,e2)-nndrl/2
                    elem%subcnc(2)%array(3)=topo(1,e4)-nndrl/2
                    elem%subcnc(2)%array(4)=topo(2,e4)-nndrl/2
                    elem%subcnc(2)%array(5)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(2)%array(5)=jnode2
                    elem%subcnc(2)%array(6)=topo(1,e3)
                    elem%subcnc(2)%array(7)=topo(2,e3) 
                    
                    subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))
                    
                    mnode(1)=tratio2*lib_node(subglbcnc(2)%array(1))+(one-tratio2)*lib_node(subglbcnc(2)%array(2))
                    mnode(2)=lib_node(subglbcnc(2)%array(2))
                    mnode(3)=lib_node(subglbcnc(2)%array(3))
                    mnode(4)=lib_node(subglbcnc(2)%array(5))
                    mnode(5)=lib_node(subglbcnc(2)%array(6))
                    mnode(6)=lib_node(subglbcnc(2)%array(7))
                    
                    Tmatrix=zero
                    Tmatrix(1,1)=tratio2
                    Tmatrix(1,2)=one-tratio2
                    Tmatrix(2,2)=one
                    Tmatrix(3,3)=one
                    
                    call prepare(elem%subelem(2),eltype='coh3d6',matkey=elem%matkey,glbcnc=subglbcnc(2)%array &
                    & ,Tmatrix=Tmatrix,mnode=mnode)
                    
                    
                    !*** sub elm 3 connec
                    elem%subcnc(3)%array(1)=topo(1,e4)-nndrl/2 ! upper surf nodes
                    elem%subcnc(3)%array(2)=topo(2,e4)-nndrl/2
                    elem%subcnc(3)%array(3)=topo(1,e2)-nndrl/2
                    elem%subcnc(3)%array(4)=topo(2,e2)-nndrl/2
                    elem%subcnc(3)%array(5)=topo(1,e4)
                    elem%subcnc(3)%array(6)=topo(2,e4)
                    elem%subcnc(3)%array(7)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(7)=jnode1
                    
                    subglbcnc(3)%array(:)=elem%nodecnc(elem%subcnc(3)%array(:))
                    
                    mnode(1)=lib_node(subglbcnc(3)%array(1))
                    mnode(2)=lib_node(subglbcnc(3)%array(2))
                    mnode(3)=tratio1*lib_node(subglbcnc(3)%array(2))+(one-tratio1)*lib_node(subglbcnc(3)%array(3))
                    mnode(4)=lib_node(subglbcnc(3)%array(5))
                    mnode(5)=lib_node(subglbcnc(3)%array(6))
                    mnode(6)=lib_node(subglbcnc(3)%array(7))
                    
                    Tmatrix=zero
                    Tmatrix(1,1)=one
                    Tmatrix(2,2)=one
                    Tmatrix(3,2)=tratio1
                    Tmatrix(3,3)=one-tratio1
                    
                    call prepare(elem%subelem(3),eltype='coh3d6',matkey=elem%matkey,glbcnc=subglbcnc(3)%array &
                    & ,Tmatrix=Tmatrix,mnode=mnode)
                    
                    
                    !*** sub elm 4 connec
                    elem%subcnc(4)%array(1)=topo(1,e1)-nndrl/2 ! upper surf nodes
                    elem%subcnc(4)%array(2)=topo(2,e1)-nndrl/2
                    elem%subcnc(4)%array(3)=topo(1,e3)-nndrl/2
                    elem%subcnc(4)%array(4)=topo(2,e3)-nndrl/2
                    elem%subcnc(4)%array(5)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(4)%array(5)=jnode1
                    elem%subcnc(4)%array(6)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(4)%array(6)=jnode2
                    elem%subcnc(4)%array(7)=topo(2,e3)             
                    
                    subglbcnc(4)%array(:)=elem%nodecnc(elem%subcnc(4)%array(:))
                    
                    mnode(1)=tratio1*lib_node(subglbcnc(4)%array(1))+(one-tratio1)*lib_node(subglbcnc(4)%array(2))
                    mnode(2)=tratio2*lib_node(subglbcnc(4)%array(2))+(one-tratio2)*lib_node(subglbcnc(4)%array(3))
                    mnode(3)=lib_node(subglbcnc(4)%array(4))
                    mnode(4)=lib_node(subglbcnc(4)%array(5))
                    mnode(5)=lib_node(subglbcnc(4)%array(6))
                    mnode(6)=lib_node(subglbcnc(4)%array(7))
                    
                    Tmatrix=zero
                    Tmatrix(1,1)=tratio1
                    Tmatrix(1,2)=one-tratio1
                    Tmatrix(2,2)=tratio2
                    Tmatrix(2,3)=one-tratio2
                    Tmatrix(3,4)=one
                    
                    call prepare(elem%subelem(4),eltype='coh3d6',matkey=elem%matkey,glbcnc=subglbcnc(4)%array &
                    & ,Tmatrix=Tmatrix,mnode=mnode)
                    
                    
                    
                    
                case default
                    write(msg_file,*)'wrong nbulk in update subcnc subxcoh'
                    call exit_function
            end select

          
!        case(6)
        ! do not update partition
            

!        case(8)


           
        case default
           write(msg_file,*) 'WARNING: subxcoh update subcnc case selection default!'
           
        end select
        
        
        ! deallocate local array
        
        if(allocated(subglbcnc)) deallocate(subglbcnc)
        if(allocated(mnode)) deallocate(mnode)
        if(allocated(Tmatrix)) deallocate(Tmatrix)
        if(allocated(x1)) deallocate(x1)
        if(allocated(x2)) deallocate(x2)
        if(allocated(xc)) deallocate(xc)


    end subroutine update_subcnc


  
  
  
    subroutine update_mnode(elem)
    
    	! passed-in variables
    	type(subxcoh_element),    intent(inout)   :: elem
    	
    	! local variables
    	type(int_alloc_array), allocatable :: subglbcnc(:)  ! glb cnc of sub elements
    	integer :: i, j, l                                  ! counters
    	integer :: nsub                             		! no. of sub elements


    	real(dp)                :: tratio1, tratio2 		! tratio=|xc-x1|/|x2-x1|
    	type(xnode),allocatable :: mnode(:)                 ! material nodes of cohesive elem
    	real(dp), allocatable   :: Tmatrix(:,:)             ! interpolation matrix btw bottom num nodes and mat nodes

!       initialize local variables

        i=0; j=0; l=0
        nsub=0
        tratio1=zero; tratio2=zero
    
    	if(.not.allocated(elem%subelem)) then
    		write(msg_file,*)'sub elems not allocated in subxcoh update mnode!'
    		call exit_function
    	end if
 
        ! currently only support elem partition elfail2
        select case(elem%curr_status)
        
            case(elfail2)
   	
                nsub=size(elem%subelem)
                
                select case(nsub)
                        case(2)
                        ! two coh3d8 sub elems

                            allocate(subglbcnc(nsub))
                            
                            ! allocate 8 numerical nodes for sub elems
                            do j=1, nsub
                                allocate(subglbcnc(j)%array(8))
                                subglbcnc(j)%array=0
                                subglbcnc(j)%array(:)=elem%nodecnc(elem%subcnc(j)%array(:))
                            end do
                            
                            !*** sub elem 1 mnode update
                            
                            call extract(elem%subelem(1),Tmatrix=Tmatrix,mnode=mnode)                                       
                            
                            ! extract ratio values from Tmatrix
                            tratio1=Tmatrix(2,1)
                            tratio2=Tmatrix(3,3)
                            
                            ! update the mnode array
                            mnode(1)=lib_node(subglbcnc(1)%array(1))
                            mnode(2)=tratio1*lib_node(subglbcnc(1)%array(1))+(one-tratio1)*lib_node(subglbcnc(1)%array(2))
                            mnode(3)=tratio2*lib_node(subglbcnc(1)%array(3))+(one-tratio2)*lib_node(subglbcnc(1)%array(4))
                            mnode(4)=lib_node(subglbcnc(1)%array(4))
                            mnode(5)=lib_node(subglbcnc(1)%array(5))
                            mnode(6)=lib_node(subglbcnc(1)%array(6))
                            mnode(7)=lib_node(subglbcnc(1)%array(7))
                            mnode(8)=lib_node(subglbcnc(1)%array(8))                  
                            
                            call update(elem%subelem(1),mnode=mnode)
                                             
                            
                            !*** sub elm 2 mnode update
                            
                            call extract(elem%subelem(2),Tmatrix=Tmatrix,mnode=mnode) 
                            
                            tratio1=Tmatrix(1,1)
                            tratio2=Tmatrix(4,3)
                            
                            ! update the mnode array
                            mnode(1)=tratio1*lib_node(subglbcnc(2)%array(1))+(one-tratio1)*lib_node(subglbcnc(2)%array(2))
                            mnode(2)=lib_node(subglbcnc(2)%array(2))
                            mnode(3)=lib_node(subglbcnc(2)%array(3))
                            mnode(4)=tratio2*lib_node(subglbcnc(2)%array(3))+(one-tratio2)*lib_node(subglbcnc(2)%array(4))
                            mnode(5)=lib_node(subglbcnc(2)%array(5))
                            mnode(6)=lib_node(subglbcnc(2)%array(6))
                            mnode(7)=lib_node(subglbcnc(2)%array(7))
                            mnode(8)=lib_node(subglbcnc(2)%array(8))
                            
                            call update(elem%subelem(2),mnode=mnode)                
                            
                            
                            
                        case(4)
                        ! four coh3d6 sub elems


                            if(allocated(subglbcnc)) deallocate(subglbcnc)
                            allocate(subglbcnc(nsub))
                            do j=1, nsub   ! coh3d6 elem cnc
                                allocate(subglbcnc(j)%array(7))
                                subglbcnc(j)%array=0
                                subglbcnc(j)%array(:)=elem%nodecnc(elem%subcnc(j)%array(:))
                            end do

                                        
                            
                            !*** sub elm 1 connec
                            
                            call extract(elem%subelem(1),Tmatrix=Tmatrix,mnode=mnode)
                            
                            tratio1=Tmatrix(1,1)
                            tratio2=Tmatrix(3,2)
                            
                            mnode(1)=tratio1*lib_node(subglbcnc(1)%array(1))+(one-tratio1)*lib_node(subglbcnc(1)%array(2))
                            mnode(2)=lib_node(subglbcnc(1)%array(2))
                            mnode(3)=tratio2*lib_node(subglbcnc(1)%array(2))+(one-tratio2)*lib_node(subglbcnc(1)%array(3))
                            mnode(4)=lib_node(subglbcnc(1)%array(5))
                            mnode(5)=lib_node(subglbcnc(1)%array(6))
                            mnode(6)=lib_node(subglbcnc(1)%array(7))
                            
                            call update(elem%subelem(1),mnode=mnode)
                            
         
                            
                            !*** sub elm 2 connec
                            
                            call extract(elem%subelem(2),Tmatrix=Tmatrix,mnode=mnode)
                            
                            tratio2=Tmatrix(1,1)
                            
                            mnode(1)=tratio2*lib_node(subglbcnc(2)%array(1))+(one-tratio2)*lib_node(subglbcnc(2)%array(2))
                            mnode(2)=lib_node(subglbcnc(2)%array(2))
                            mnode(3)=lib_node(subglbcnc(2)%array(3))
                            mnode(4)=lib_node(subglbcnc(2)%array(5))
                            mnode(5)=lib_node(subglbcnc(2)%array(6))
                            mnode(6)=lib_node(subglbcnc(2)%array(7))        
                            
                            call update(elem%subelem(2),mnode=mnode)

                            
                            
                            
                            !*** sub elm 3 connec
                            
                            call extract(elem%subelem(3),Tmatrix=Tmatrix,mnode=mnode)
                            
                            tratio1=Tmatrix(3,2)
                            
                            mnode(1)=lib_node(subglbcnc(3)%array(1))
                            mnode(2)=lib_node(subglbcnc(3)%array(2))
                            mnode(3)=tratio1*lib_node(subglbcnc(3)%array(2))+(one-tratio1)*lib_node(subglbcnc(3)%array(3))
                            mnode(4)=lib_node(subglbcnc(3)%array(5))
                            mnode(5)=lib_node(subglbcnc(3)%array(6))
                            mnode(6)=lib_node(subglbcnc(3)%array(7))
              
                            call update(elem%subelem(3),mnode=mnode)
                            
                            
                            
                            !*** sub elm 4 connec
                            
                            call extract(elem%subelem(4),Tmatrix=Tmatrix,mnode=mnode)
                            
                            tratio1=Tmatrix(1,1)
                            tratio2=Tmatrix(2,2)
                            
                            mnode(1)=tratio1*lib_node(subglbcnc(4)%array(1))+(one-tratio1)*lib_node(subglbcnc(4)%array(2))
                            mnode(2)=tratio2*lib_node(subglbcnc(4)%array(2))+(one-tratio2)*lib_node(subglbcnc(4)%array(3))
                            mnode(3)=lib_node(subglbcnc(4)%array(4))
                            mnode(4)=lib_node(subglbcnc(4)%array(5))
                            mnode(5)=lib_node(subglbcnc(4)%array(6))
                            mnode(6)=lib_node(subglbcnc(4)%array(7))          
                            
                            call update(elem%subelem(4),mnode=mnode)                   
                            
                            
                        case default
                            write(msg_file,*)'wrong nbulk in update subcnc subxcoh'
                            call exit_function
                end select

            case default
                write(msg_file,*)'partition status not supported in subxcoh update mnode!'
                call exit_function
                
        end select
            
        ! deallocate local array
        
        if(allocated(subglbcnc)) deallocate(subglbcnc)
        if(allocated(mnode)) deallocate(mnode)
        if(allocated(Tmatrix)) deallocate(Tmatrix)
  
  
  	end subroutine update_mnode
  
  
  
  
  
    subroutine update_sdv(elem)
    
    	! passed-in variables
    	type(subxcoh_element),    intent(inout)   :: elem
    	
    	! local variables   
        character(len=eltypelength) ::  subeltype
    
        ! real and integer failure variables
        real(dp) :: rfvar
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
        rfvar=zero; ifvar=0
        i=0; j=0; l=0
        
    
        ! update subxcoh elem sdv for output
        if(.not.allocated(elem%sdv)) then
            allocate(elem%sdv(1))
            allocate(elem%sdv(1)%i(1))  ! fstat
            allocate(elem%sdv(1)%r(1))  ! dm
            elem%sdv(1)%i(1)=0
            elem%sdv(1)%r(1)=0
        end if
                    
                    
        do i=1,size(elem%subelem)           
        
            ifvar=0
            rfvar=zero
            
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
                    write(msg_file,*)'subelem type not supported in subxcoh elem!'
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
            
            ! update this fvar to subxcoh sdv
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
  
  
end module subxcoh_element_module
