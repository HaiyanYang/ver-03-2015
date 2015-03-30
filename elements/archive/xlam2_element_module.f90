module xlam_element_module
    use parameter_module
    use toolkit_module                  ! global tools for element integration
    use lib_edge_module                 ! global edge library
    use lib_node_module                 ! global node library
    use lib_mat_module                  ! global material library
    use lib_bcd_module                  ! global bcd library
    use xbrick_element_module
    use xcoh_element_module

    implicit none
    private
    
    ! parameters for no. nodes of ply-block element (xbrick) and interface element (coh3d8)
    integer, parameter :: ndim=3, nndplyblk=24, nedgplyblk=8, nndinterf=24, nedginterf=8, ncorner=8
    

    type, public :: xlam_element
        private
        
        integer :: curr_status=0        ! 0 means intact
        integer :: key=0 
        integer :: bulkmat=0
        integer :: cohmat=0
        integer :: interfmat=0
        
        integer,allocatable :: nodecnc(:)    ! cnc to glb node arrays for accessing nodal variables (x, u, du, v, dof ...)
        integer,allocatable :: edgecnc(:)
        real(dp),allocatable:: layup(:,:)   ! angle and relative thickness of all plies
    
        type(xbrick_element),allocatable :: plyblk(:)
        type(xcoh_element),allocatable :: interf(:)
        
        type(int_alloc_array), allocatable :: plyblknodecnc(:), plyblkedgecnc(:)      ! plyblk_elem connec to parent elem nodes
        type(int_alloc_array), allocatable :: interfnodecnc(:), interfedgecnc(:)      ! interf_elem connec to parent elem nodes
    
    end type xlam_element
    
    interface empty
        module procedure empty_xlam_element
    end interface
  
    interface prepare
        module procedure prepare_xlam_element
    end interface
    
    interface integrate
        module procedure integrate_xlam_element
    end interface
    
    interface extract
        module procedure extract_xlam_element
    end interface




    public :: empty,prepare,integrate,extract
    
    
    
    
    contains
    
    
    
    
    
    subroutine empty_xlam_element(elem)
  
        type(xlam_element),intent(out) :: elem
        
        elem%curr_status=0
        elem%key=0 
        elem%bulkmat=0
        elem%cohmat=0
        elem%interfmat=0
        
        if(allocated(elem%nodecnc)) deallocate(elem%nodecnc)
        if(allocated(elem%edgecnc)) deallocate(elem%edgecnc)
        if(allocated(elem%layup)) deallocate(elem%layup)
        if(allocated(elem%plyblk)) deallocate(elem%plyblk)
        if(allocated(elem%interf)) deallocate(elem%interf)
        if(allocated(elem%plyblknodecnc)) deallocate(elem%plyblknodecnc)
        if(allocated(elem%plyblkedgecnc)) deallocate(elem%plyblkedgecnc)
        if(allocated(elem%interfnodecnc)) deallocate(elem%interfnodecnc)
        if(allocated(elem%interfedgecnc)) deallocate(elem%interfedgecnc)

    end subroutine empty_xlam_element
  
    
    
    
    
    
    subroutine extract_xlam_element(elem,curr_status,key,bulkmat,cohmat,interfmat,nodecnc,edgecnc,layup &
    & ,plyblk,interf,plyblknodecnc,plyblkedgecnc,interfnodecnc,interfedgecnc)
    
        type(xlam_element),intent(in)  :: elem
        
        integer, optional, intent(out) :: curr_status
        integer, optional, intent(out) :: key 
        integer, optional, intent(out) :: bulkmat
        integer, optional, intent(out) :: cohmat
        integer, optional, intent(out) :: interfmat
        
        integer,allocatable, optional, intent(out) :: nodecnc(:) 
        integer,allocatable, optional, intent(out) :: edgecnc(:)
        real(dp),allocatable, optional, intent(out):: layup(:,:)
    
        type(xbrick_element),allocatable, optional, intent(out) :: plyblk(:)
        type(xcoh_element),allocatable, optional, intent(out) :: interf(:)
        
        type(int_alloc_array), allocatable, optional, intent(out) :: plyblknodecnc(:), plyblkedgecnc(:)
        type(int_alloc_array), allocatable, optional, intent(out) :: interfnodecnc(:), interfedgecnc(:)
        
        
        
        if(present(curr_status)) curr_status=elem%curr_status
        if(present(key)) key=elem%key 
        if(present(bulkmat)) bulkmat=elem%bulkmat
        if(present(cohmat)) cohmat=elem%cohmat
        if(present(interfmat)) interfmat=elem%interfmat
        
        if(present(nodecnc)) then 
            if(allocated(elem%nodecnc)) then
                allocate(nodecnc(size(elem%nodecnc)))
                nodecnc=elem%nodecnc
            end if
        end if
        
        if(present(edgecnc)) then 
            if(allocated(elem%edgecnc)) then
                allocate(edgecnc(size(elem%edgecnc)))
                edgecnc=elem%edgecnc
            end if
        end if
        
        if(present(layup)) then
            if(allocated(elem%layup)) then
                allocate(layup(size(elem%layup(:,1)),size(elem%layup(1,:))))
                layup=elem%layup
            end if
        end if
        
        if(present(plyblk)) then
            if(allocated(elem%plyblk)) then
                allocate(plyblk(size(elem%plyblk)))
                plyblk=elem%plyblk
            end if
        end if
        
        if(present(interf)) then
            if(allocated(elem%interf)) then
                allocate(interf(size(elem%interf)))
                interf=elem%interf
            end if
        end if
        
        if(present(plyblknodecnc)) then
            if(allocated(elem%plyblknodecnc)) then
                allocate( plyblknodecnc(size(elem%plyblknodecnc)))
                plyblknodecnc=elem%plyblknodecnc
            end if
        end if
        
        if(present(plyblkedgecnc)) then
            if(allocated(elem%plyblkedgecnc)) then
                allocate( plyblkedgecnc(size(elem%plyblkedgecnc)))
                plyblkedgecnc=elem%plyblkedgecnc
            end if
        end if
        
        if(present(interfnodecnc)) then
            if(allocated(elem%interfnodecnc)) then
                allocate( interfnodecnc(size(elem%interfnodecnc)))
                interfnodecnc=elem%interfnodecnc
            end if
        end if
        
        if(present(interfedgecnc)) then
            if(allocated(elem%interfedgecnc)) then
                allocate( interfedgecnc(size(elem%interfedgecnc)))
                interfedgecnc=elem%interfedgecnc
            end if
        end if
    
    end subroutine extract_xlam_element
    
    
    
    
    
    subroutine prepare_xlam_element(elem,key,bulkmat,cohmat,interfmat,nodecnc,edgecnc,layup)
    
        type(xlam_element),    intent(inout)    :: elem
        integer,                intent(in)      :: key
        integer,                intent(in)      :: bulkmat, cohmat, interfmat
        integer,                intent(in)      :: nodecnc(:)
        integer,                intent(in)      :: edgecnc(:)
        real(dp),               intent(in)      :: layup(:,:)

        elem%key=key 
        elem%bulkmat=bulkmat
        elem%cohmat=cohmat
        elem%interfmat=interfmat
        
        if (size(nodecnc)/=nndplyblk*size(layup(1,:))) then
            write(msg_file,*)'layup and no. nodes do not match in prepare_xlam!'
            call exit_function
        end if
        
        if (size(edgecnc)/=nedgplyblk*size(layup(1,:))) then
            write(msg_file,*)'layup and no. edges do not match in prepare_xlam!'
            call exit_function
        end if
        
        allocate(elem%nodecnc(size(nodecnc)))
        allocate(elem%edgecnc(size(edgecnc)))
        allocate(elem%layup(size(layup(:,1)),size(layup(1,:))))
        elem%nodecnc=nodecnc
        elem%edgecnc=edgecnc
        elem%layup=layup
    
    end subroutine prepare_xlam_element
    
    
    


    subroutine integrate_xlam_element(elem, K_matrix, F_vector)
    
        type(xlam_element),intent(inout)       :: elem 
        real(kind=dp),allocatable,intent(out)   :: K_matrix(:,:), F_vector(:)
    
    
        ! local variables
        
        ! sub_elem K matrix and F vector
        real(kind=dp),allocatable           :: Ki(:,:), Fi(:)
        
        ! loop counters
        integer :: i,j,l  

        ! no. dof, no. plyblocks and no. interfaces in this elem
        integer :: ndof, nplyblk, ninterf
        
        ! lcl arrays to store temporarily the glb cnc of plyblock nodes, edges, and interface nodes
        integer :: plyblknode(nndplyblk), plyblkedge(nedgplyblk), interfnode(nndinterf), interfedge(nedginterf) 
        
        ! lcl arrays to store temporarily the cnc of dofs of each sub elem to parent elem dofs
        integer, allocatable :: dofcnc(:)
        
        ! temp. coordinates of corner nodes (used to update the z-coord of corner nodes of each plyblk)
        type(real_alloc_array) :: coord(ncorner)
        
        ! corner edge thickness of this laminate (no. of corner edges = half * no. of corner nodes)
        real(dp) :: shellthickness(ncorner/2)
        
        ! penalty stiffness to enforce boundary conditions
        real(dp) :: Kpn
        
        ! variables for imposing penalty constrains on bcd nodes
        integer :: nd1, nd2, dof0, dof1, dof2
        real(dp), allocatable :: u0(:), u1(:), u2(:)
        
        ! interface status variable
        integer :: interfstat
        
        ! interface no. failed edges and indices of failed edges
        integer :: nfe, ifailedge(nedginterf)
        
        ! ifailedge array of bottom and top plyblks of an interface
        integer, allocatable :: ifedg1(:), ifedg2(:)
        
        
        ! initialize local variables
        i=0; j=0; l=0
        ndof=0; nplyblk=0; ninterf=0
        plyblknode=0; plyblkedge=0; interfnode=0; interfedge=0
        shellthickness=zero
        Kpn=zero
        nd1=0; nd2=0; dof0=0; dof1=0; dof2=0
        interfstat=0; nfe=0; ifailedge=0
        
        ! penalty stiffness = 1GPa
        Kpn=1000000._dp
        
        
        ! check if the elem has been prepared (by checking layup; note that in prepare subroutine all attributes 
        ! must be allocated together, so checking one of them (here layup) would suffice)
        if(.not.allocated(elem%layup)) then
            write(msg_file,*)'xlam element must be prepared before integration'
            call exit_function
        end if

        ! extract no. plyblock and no. interfaces from layup, and calculate ndof
        nplyblk=size(elem%layup(1,:))
        ninterf=nplyblk-1
        ndof=ndim*size(elem%nodecnc)

     

        !---------------------------------------------------------------------!
        !       prepare sub elements
        !---------------------------------------------------------------------! 
        if (.not.allocated(elem%plyblk)) then
        ! if plyblock elems (& consequently interface elems) not yet allocated (beginning of analysis), 
        ! allocate the plyblock and interface elems first
            
            ! allocate plyblk elems and interface elems
            allocate(elem%plyblk(nplyblk))
            allocate(elem%interf(ninterf))
            
            ! allocate plyblk node and edge cnc arrays...
            allocate(elem%plyblknodecnc(nplyblk))
            allocate(elem%plyblkedgecnc(nplyblk))

            
            do i=1, nplyblk
                allocate(elem%plyblknodecnc(i)%array(nndplyblk))
                allocate(elem%plyblkedgecnc(i)%array(nedgplyblk))
                elem%plyblknodecnc(i)%array=[( j, j=(i-1)*nndplyblk+1 , i*nndplyblk )]
                elem%plyblkedgecnc(i)%array=[( j, j=(i-1)*nedgplyblk+1 , i*nedgplyblk )]
            end do
            
            ! ...and interface node cnc arrays
            allocate(elem%interfnodecnc(ninterf))
            allocate(elem%interfedgecnc(ninterf))
            
            do i=1, ninterf
                allocate(elem%interfnodecnc(i)%array(nndinterf))  
                allocate(elem%interfedgecnc(i)%array(nedginterf))  
                
                ! following needs to be filled (real nodes first, then flo nodes)
                           
                ! 1st half of interface real nodes comes from bottom plyblk elem top surface (nodes 5-8)
                !elem%interfnodecnc(i)%array(1 : ncorner/2)=&
                !& [( j, j=(i-1)*nndplyblk+ncorner/2+1 , (i-1)*nndplyblk+ncorner )]
                elem%interfnodecnc(i)%array(1 : ncorner/2)=&
                & elem%plyblknodecnc(i)%array(ncorner/2+1 : ncorner)
                
                ! 2nd half of interface real nodes comes from top plyblk elem bottom surface (nodes 1-4)
                !elem%interfnodecnc(i)%array(ncorner/2+1 : ncorner)=&
                !& [( j, j=i*nndplyblk+1 , i*nndplyblk+ncorner/2 )]    
                elem%interfnodecnc(i)%array(ncorner/2+1 : ncorner)=&
                & elem%plyblknodecnc(i+1)%array(1 : ncorner/2)
                
                ! 1st half of interface flo nodes come from bottm plyblk elem top surface (nodes 17-24)
                !elem%interfnodecnc(i)%array(ncorner+1 : ncorner+(nndinterf-ncorner)/2)=&
                !& [( j, j=(i-1)*nndplyblk+ncorner+(nndinterf-ncorner)/2+1 , (i-1)*nndplyblk+nndinterf)]
                elem%interfnodecnc(i)%array(ncorner+1 : ncorner+(nndinterf-ncorner)/2)=&
                & elem%plyblknodecnc(i)%array(ncorner+(nndinterf-ncorner)/2+1 : nndinterf)
                
                ! 2nd half of interface flo nodes come from top plyblk elem bottom surface (nodes 9-16)
                !elem%interfnodecnc(i)%array(ncorner+(nndinterf-ncorner)/2+1 : nndinterf)=&
                !& [( j, j=i*nndplyblk+ncorner+1 , i*nndplyblk+ncorner+(nndinterf-ncorner)/2 )] 
                elem%interfnodecnc(i)%array(ncorner+(nndinterf-ncorner)/2+1 : nndinterf)=&
                & elem%plyblknodecnc(i+1)%array(ncorner+1 : ncorner+(nndinterf-ncorner)/2)
                
                ! interface bottom edges are from bottom plyblk elem top edges 
                elem%interfedgecnc(i)%array(1 : nedginterf/2)=elem%plyblkedgecnc(i)%array(nedginterf/2+1 : nedginterf)
                
                ! interface top edges are from top plyblk elem bottom edges
                elem%interfedgecnc(i)%array(nedginterf/2+1 : nedginterf)=elem%plyblkedgecnc(i+1)%array(1 : nedginterf/2)
                
                       
            end do
            
            
            
            ! extract elem corner nodes' coords from glb node library
            do j=1, ncorner ! first ncorner nodes of this elem
                call extract(lib_node(elem%nodecnc(j)),x=coord(j)%array)
                if(.not.allocated(coord(j)%array)) then
                    write(msg_file,*)'error: corner node coords undefined!'
                    call exit_function
                end if
            end do
            ! compute shell thickness at the corner edges (ncorner/2 edges)
            do j=1, ncorner/2
                shellthickness(j)=coord(j+ncorner/2)%array(3)-coord(j)%array(3)
            end do


            
            ! prepare plyblk elems
            do i=1, nplyblk
                ! extract the glb node and edge cnc of this plyblock element from elem glb cnc and plyblk i's local cnc
                plyblknode(:)=elem%nodecnc(elem%plyblknodecnc(i)%array(:))
                plyblkedge(:)=elem%edgecnc(elem%plyblkedgecnc(i)%array(:))
                
                
                ! update corner node coords (z coords) of this plyblk accord. to plyblk thickness ratio
                ! coord(:) is updated to store corner nodes coords of each plyblk
                if (i==1) then
                    do j=1, ncorner/2
                        ! 1st ply top nodes = 1st ply bottom nodes + 1stplyratio* shellthickness
                        coord(j+ncorner/2)%array(3)=coord(j)%array(3)+elem%layup(2,i)*shellthickness(j)
                        ! update 1st ply top node coords into glb lib_node array
                        call update(lib_node(plyblknode(j+ncorner/2)),x=coord(j+ncorner/2)%array)
                    end do
                else
                    do j=1, ncorner/2
                        ! ith ply bottom nodes = (i-1)th ply top nodes (already stored in coord(j+ncorner/2))
                        coord(j)%array(3)=coord(j+ncorner/2)%array(3)
                        ! update ith ply bottom node coords into glb lib_node array
                        call update(lib_node(plyblknode(j)),x=coord(j)%array)
                        ! ith ply top nodes = ith ply bottom nodes + ithplyratio*shellthickness
                        coord(j+ncorner/2)%array(3)=coord(j)%array(3)+elem%layup(2,i)*shellthickness(j)
                        ! update ith ply top node coords into glb lib_node array
                        call update(lib_node(plyblknode(j+ncorner/2)),x=coord(j+ncorner/2)%array)
                    end do                    
                end if
                
                ! prepare each plyblk elem (here xbrick elem type)
                call prepare(elem%plyblk(i),key=0,bulkmat=elem%bulkmat,cohmat=elem%cohmat, &
                & plyangle=elem%layup(1,i), nodecnc=plyblknode,edgecnc=plyblkedge)     
          
            end do
            
            ! prepare interf elems
            do i=1, ninterf
                ! extract the glb node cnc of this interface element from elem glb cnc and interface i's local cnc
                interfnode(:)=elem%nodecnc(elem%interfnodecnc(i)%array(:))
                interfedge(:)=elem%edgecnc(elem%interfedgecnc(i)%array(:))
                
                ! prepare each interface elem (here xcoh elem type)
                !call prepare(elem%interf(i),key=0,connec=interfnode,matkey=elem%interfmat)
                call prepare(elem%interf(i),key=0,matkey=elem%interfmat,nodecnc=interfnode,edgecnc=interfedge)
            end do
        
        end if
        


        !---------------------------------------------------------------------!
        !       integrate and assemble sub element system arrays
        !---------------------------------------------------------------------!

        ! initialize K & F
        allocate(K_matrix(ndof,ndof),F_vector(ndof))
        K_matrix=zero; F_vector=zero 
     
        ! integrate plyblock elements and assemble into global matrix

        do i=1, nplyblk

            call integrate(elem%plyblk(i),Ki,Fi)

            if(allocated(dofcnc)) deallocate(dofcnc)
            allocate(dofcnc(size(Fi))); dofcnc=0
            
            do j=1, nndplyblk ! no. of nodes in sub elem i
                do l=1, ndim
                    ! dof indices of the jth node of sub elem i 
                    dofcnc((j-1)*ndim+l)=(elem%plyblknodecnc(i)%array(j)-1)*ndim+l
                end do
            end do
            call assembleKF(K_matrix,F_vector,Ki,Fi,dofcnc)
            deallocate(Ki)
            deallocate(Fi)
            deallocate(dofcnc)
        end do

        ! integrate cohesive elements and assemble into global matrix

        do i=1, ninterf
        
        	! extract status of this interface
        	call extract(elem%interf(i),curr_status=interfstat)
        	
        	! if this interface elem has not yet reached final partition,
        	! update its ifailedge array before integration
        	if(interfstat<elfail3) then 
        
        		! extract failed edges from bottom and top plyblk elems
        		call extract(elem%plyblk(i),ifailedge=ifedg1)
        		call extract(elem%plyblk(i+1),ifailedge=ifedg2)
        	
        		! no of failed edges in this interface
        		nfe=0
        		ifailedge=0
        	
        		! pass bottom plyblk failed edge info into this interface ifailedge
        		do j=1, size(ifedg1)
        			select case (ifedg1(j))
        				! upper edges failed: pass into interface ifailedge as bottom edges
        				case(5:8)
        					nfe=nfe+1
        					ifailedge(nfe)=ifedg1(j)-4	
        				! bottom edges failed: ignored	
        				case(0:4)
        					continue
        				case default
        					write(msg_file,*)'wrong failed edge index in xlam integration'
        					call exit_function
        			end select
        		end do
        	
        		! pass top plyblk failed edge info into this interface ifailedge
        		do j=1, size(ifedg2)
        			select case (ifedg2(j))
        				! bottom edges failed: pass into interface ifailedge as top edges
        				case(1:4)
        					nfe=nfe+1
        					ifailedge(nfe)=ifedg2(j)+4	
        				! top edges failed: ignored	
        				case(0,5:8)
        					continue
        				case default
        					write(msg_file,*)'wrong failed edge index in xlam integration'
        					call exit_function
        			end select
        		end do
        	
        		! update ifailedge array into this interface elem
        		call update(elem%interf(i),ifailedge=ifailedge)
        	
        	end if
        	
            call integrate(elem%interf(i),Ki,Fi)
            if(allocated(dofcnc)) deallocate(dofcnc)
            allocate(dofcnc(size(Fi))); dofcnc=0
            
            do j=1, nndinterf ! no. of nodes in sub elem i
                do l=1, ndim
                    ! dof indices of the jth node of sub elem i 
                    dofcnc((j-1)*ndim+l)=(elem%interfnodecnc(i)%array(j)-1)*ndim+l
                end do
            end do
            call assembleKF(K_matrix,F_vector,Ki,Fi,dofcnc)
            deallocate(Ki)
            deallocate(Fi)
            deallocate(dofcnc)
        end do


        ! check if nodes are applied bcd
        do i=1, ncorner/2
            ! check if the bcd nodes contain vertical edges of the elem
            if(any(lib_bcdnodes==elem%nodecnc(i)) .and. any(lib_bcdnodes==elem%nodecnc(i+ncorner/2))) then
                ! constrain all nodes along this vertical edge to the disp. of these two nodes
                if(nplyblk > 1) then
                    do j=2, nplyblk
                        ! find corresponding nodes along this vertical edge of all plyblks
                        nd1=i+nndplyblk*(j-1)
                        nd2=i+ncorner/2+nndplyblk*(j-1)   
                        
                        
                        ! constrain node 1 to node i (use penalty stiffness Kpn)
                        dof0=(i-1)*ndim
                        dof1=(nd1-1)*ndim
                        call extract(lib_node(elem%nodecnc(i)),u=u0)
                        call extract(lib_node(elem%nodecnc(nd1)),u=u1)
                        do l=1, ndim
                            K_matrix(dof0+l,dof0+l)=K_matrix(dof0+l,dof0+l)+Kpn
                            K_matrix(dof1+l,dof1+l)=K_matrix(dof1+l,dof1+l)+Kpn
                            K_matrix(dof0+l,dof1+l)=K_matrix(dof0+l,dof1+l)-Kpn
                            K_matrix(dof1+l,dof0+l)=K_matrix(dof1+l,dof0+l)-Kpn
                            F_vector(dof0+l)=F_vector(dof0+l)+Kpn*(u0(l)-u1(l))
                            F_vector(dof1+l)=F_vector(dof1+l)-Kpn*(u0(l)-u1(l))
                        end do
                         
                        ! constrain node 2 to node i+ncorner/2
                        dof0=(i+ncorner/2-1)*ndim
                        dof2=(nd2-1)*ndim
                        call extract(lib_node(elem%nodecnc(i+ncorner/2)),u=u0)
                        call extract(lib_node(elem%nodecnc(nd2)),u=u2)
                        do l=1, ndim
                            K_matrix(dof0+l,dof0+l)=K_matrix(dof0+l,dof0+l)+Kpn
                            K_matrix(dof2+l,dof2+l)=K_matrix(dof2+l,dof2+l)+Kpn
                            K_matrix(dof0+l,dof2+l)=K_matrix(dof0+l,dof2+l)-Kpn
                            K_matrix(dof2+l,dof0+l)=K_matrix(dof2+l,dof0+l)-Kpn
                            F_vector(dof0+l)=F_vector(dof0+l)+Kpn*(u0(l)-u2(l))
                            F_vector(dof2+l)=F_vector(dof2+l)-Kpn*(u0(l)-u2(l))
                        end do   
                    end do
                end if
                ! in the future, interpolate the disp. of floating nodes on bcd edges
            end if
        end do        
        
        
        !---------------------------------------------------------------------!
        !               deallocate local arrays 
        !---------------------------------------------------------------------!
        if(allocated(Ki)) deallocate(Ki)
        if(allocated(Fi)) deallocate(Fi)
        if(allocated(dofcnc)) deallocate(dofcnc)
        if(allocated(u0)) deallocate(u0)
        if(allocated(u1)) deallocate(u1)
        if(allocated(u2)) deallocate(u2)
        if(allocated(ifedg1)) deallocate(ifedg1)
        if(allocated(ifedg2)) deallocate(ifedg2)
    
 
    end subroutine integrate_xlam_element
    
    
    

end module xlam_element_module
