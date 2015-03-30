module xbrick_element_module
    use parameter_module
    use glb_clock_module
    use toolkit_module                  ! global tools for element integration
    use lib_edge_module                 ! global edge library
    use lib_node_module                 ! global node library
    use lib_mat_module                  ! global material library
    use sub3d_element_module
  
  
    implicit none
    private

    integer,parameter :: ndim=3, nndrl=8, nedge=8, nndfl=2*nedge, nnode=nndrl+nndfl, ndof=ndim*nnode
    ! Topology: nodes on each edge; 4 nodes per edge, 1-2 are end nodes, 3-4 are fl. nodes; assigned in lcl node numbers
    integer,parameter :: topo(4,nedge)=reshape([1,2,9,10,2,3,11,12,3,4,13,14,4,1,15,16, &
                                              & 5,6,17,18,6,7,19,20,7,8,21,22,8,5,23,24],[4,nedge])
                                              
    

    type, public :: xbrick_element             ! breakable brick
        private
        
        integer :: curr_status=0        ! 0 means intact
        integer :: key=0 
        integer :: bulkmat=0
        integer :: cohmat=0
        
        real(dp):: plyangle=zero        ! ply angle for composite lamina (rotation around z axis) 
        
        integer :: nodecnc(nnode)=0     ! cnc to glb node arrays for accessing nodal variables (x, u, du, v, dof ...)
        integer :: edgecnc(nedge)=0     ! cnc to glb edge arrays for accessing edge variables (failure status)
        integer :: ifailedge(nedge)=0   ! indices of failed edges
        
        logical :: newpartition=.false. ! true when elem is changing partition, no failure should be considered then
        integer :: nstep=0, ninc=0      ! to store curr step and increment no.
        
        type(sub3d_element), allocatable :: subelem(:)
        type(int_alloc_array), allocatable :: subcnc(:)      ! sub_elem connec to parent elem nodes
        
        type(sdv_array), allocatable :: sdv(:)
        
    end type xbrick_element
  
    interface empty
        module procedure empty_xbrick_element
    end interface
  
    interface prepare
        module procedure prepare_xbrick_element
    end interface
    
    !~interface precrack
    !~    module procedure precrack_xbrick_element
    !~end interface
    
    interface integrate
        module procedure integrate_xbrick_element
    end interface
    
    interface extract
        module procedure extract_xbrick_element
    end interface




    public :: empty,prepare,integrate,extract



    contains




    ! empty a breakable quadrilateral
    subroutine empty_xbrick_element(elem)
  
        type(xbrick_element),intent(out) :: elem
        
        elem%curr_status=0
        elem%key=0 
        elem%bulkmat=0
        elem%cohmat=0
        elem%plyangle=zero
        
        elem%nodecnc=0
        elem%edgecnc=0
        elem%ifailedge=0
        
        elem%newpartition=.false.
        elem%nstep=0
        elem%ninc=0
        
        if(allocated(elem%subelem)) deallocate(elem%subelem)
        if(allocated(elem%subcnc))  deallocate(elem%subcnc)
        if(allocated(elem%sdv)) deallocate(elem%sdv)

    end subroutine empty_xbrick_element
  
  
  
    ! this subroutine is used to prepare the connectivity and material lib index of the element
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine prepare_xbrick_element(elem,key,bulkmat,cohmat,plyangle,nodecnc,edgecnc)
    
        type(xbrick_element),    intent(inout)   :: elem
        integer,                intent(in)      :: key
        integer,                intent(in)      :: bulkmat, cohmat
        real(dp),               intent(in)      :: plyangle
        integer,                intent(in)      :: nodecnc(nnode)
        integer,                intent(in)      :: edgecnc(nedge)

        elem%key=key 
        elem%bulkmat=bulkmat
        elem%cohmat=cohmat
        elem%plyangle=plyangle
        elem%nodecnc=nodecnc
        elem%edgecnc=edgecnc
    
    end subroutine prepare_xbrick_element
    
    
    subroutine extract_xbrick_element(elem,curr_status,key,bulkmat,cohmat,plyangle,nodecnc,edgecnc, &
    & ifailedge,newpartition,nstep,ninc,subelem,subcnc,sdv)
    
        type(xbrick_element),                      intent(in)  :: elem
        integer,                        optional, intent(out) :: curr_status
        integer,                        optional, intent(out) :: key
        integer,                        optional, intent(out) :: bulkmat
        integer,                        optional, intent(out) :: cohmat
        real(dp),                       optional, intent(out) :: plyangle
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
        if(present(bulkmat)) bulkmat=elem%bulkmat
        if(present(cohmat)) cohmat=elem%cohmat
        if(present(plyangle)) plyangle=elem%plyangle
        
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
    
    end subroutine extract_xbrick_element






    subroutine integrate_xbrick_element(elem, K_matrix, F_vector)
    
        type(xbrick_element),intent(inout)       :: elem 
        real(kind=dp),allocatable,intent(out)   :: K_matrix(:,:), F_vector(:)
    
    
        ! local variables
        type(int_alloc_array), allocatable  :: subglbcnc(:)     ! glb cnc of sub element, used when elem is intact
        
        integer :: i,j,l, elstat, subelstat
        
        character(len=eltypelength) ::  subeltype
        
        logical :: nofailure
        
        ! - glb clock step and increment no. extracted from glb clock module
        integer :: curr_step, curr_inc
        logical :: last_converged               ! true if last iteration has converged: a new increment/step has started
    
    
        ! initialize K & F
        allocate(K_matrix(ndof,ndof),F_vector(ndof))
        K_matrix=zero; F_vector=zero
        
        ! initialize local variables
        i=0; j=0; l=0
        elstat=0; subelstat=0; subeltype=''
        
        curr_step=0; curr_inc=0
        last_converged=.false.
        
        nofailure=.false.
        
        if(.not.allocated(elem%subelem)) then 
            allocate(elem%subelem(1))
            allocate(elem%subcnc(1))
            allocate(elem%subcnc(1)%array(nndrl))   ! brick elem
            allocate(subglbcnc(1))
            allocate(subglbcnc(1)%array(nndrl))
            ! sub elm 1 connec
            elem%subcnc(1)%array=[(i, i=1,nndrl)]
            subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
            ! create sub elements
            call prepare(elem%subelem(1),eltype='brick', matkey=elem%bulkmat, &
            & plyangle=elem%plyangle, glbcnc=subglbcnc(1)%array)
        end if


        ! - extract curr step and inc values from glb clock module
        call extract_glb_clock(kstep=curr_step,kinc=curr_inc)
        
        ! - check if last iteration has converged, and update the current step & increment no.
        if(elem%nstep.ne.curr_step .or. elem%ninc.ne.curr_inc) then
            last_converged=.true.
            elem%nstep = curr_step
            elem%ninc = curr_inc
            elem%newpartition=.false.   ! last partition has converged (nolonger 'new')
        end if

        
        

        !---------------------------------------------------------------------!
        !       update elem partition using edge status variable
        !---------------------------------------------------------------------!
     
        ! if elem is not yet failed, check elem edge status variables and update elem status and sub elem cnc
        if(elem%curr_status<elfailm) then 
            ! store current status value
            elstat=elem%curr_status  
            
            ! by default, nofailure is false, i.e., failure criterion will be assessed
            nofailure=.false.
            
            ! partition elem according to edge status values
            call edge_status_partition(elem)  
         
            if(elstat/=elem%curr_status) then
                elem%newpartition=.true.    ! new partition is true for this increment
                nofailure=.true.            ! no material degradation/failure criterion partition for 1st iteration of new partition
            end if 
            

            ! if elem matrix is not yet failed after the edge status partition
            ! integrate and check the failure criterion partition
            if(elem%curr_status<elfailm) then
            
                call integrate_assemble(elem,K_matrix,F_vector,nofailure)

                if(nofailure) then
                ! 1st iteration of new partition, no failure criterion partition for stabilization purpose
                    continue
                else
                ! 2nd and later iteration of new partition
                    
                    !***** check failure criterion *****
                    ! failure criterion partitions elem into elfailm/elfailf partition if sub elem matrix fails/fibre fails
                    call failure_criterion_partition(elem)
                    
                    if(elem%curr_status>=elfailm) then
                    ! elem partition is updated by failure criterion partition, i.e., new partiton is true
                        elem%newpartition=.true.
                        nofailure=.true.            ! no material degradation/failure criterion partition for 1st iteration of new partition
                    end if
                    
                end if
                
            end if
               
        end if
        
        
       
        if(elem%curr_status==elfailm) then
        ! element matrix is already failed
        ! element is already partitioned into 2 bulks and 1 coh, integrate and assemble subelems
        ! bulk sub elems may undergo fibre damage and failure (only after coh sub elem starts failing)
            
            
            ! by default, no material degradation for fibre
            ! until cohesive sub elem starts to fail
            nofailure=.true.
            
            ! during the increment of new partition, no fibre material degradation allowed
            if(elem%newpartition) then
                continue    ! nofailure remains true
            else
            ! after that increment, fibre failure is considered for matrix fail partition only after coh sub elem starts failing
                do i=1, size(elem%subelem)
                    call extract(elem%subelem(i),eltype=subeltype,curr_status=subelstat)
                    if(subeltype=='coh3d6' .or. subeltype=='coh3d8') then
                        if(subelstat > intact) nofailure=.false.
                    end if
                end do
            end if
             
            
            
            ! integrate sub elems
            call integrate_assemble(elem,K_matrix,F_vector,nofailure)
            
            ! check if sub elems have reached fibre failure onset (only after it's allowed);
            ! if so, update curr status to fibre failure status elfailf (no need to update partition)
            if(.not.nofailure) then
                do i=1, size(elem%subelem)
                    call extract(elem%subelem(i),curr_status=subelstat)
                    if(subelstat>=fibre_onset) then
                        elem%curr_status=elfailf
                        goto 10
                        exit
                    end if 
                end do
            end if
        
        end if
        


        if(elem%curr_status==elfailf) then
        ! element fibre is already failed, integrate and assemble subelem
            
            ! during the increment of new partition, no fibre failure allowed
            if(elem%newpartition) then
                nofailure=.true.  ! no fibre failure modelling
            else
            ! after that increment, fibre failure is considered for fibre fail partition
                nofailure=.false.
            end if    
            
            call integrate_assemble(elem,K_matrix,F_vector,nofailure)

        end if
       

        
        !---------------------------------------------------------------------!
        !               deallocate local arrays 
        !---------------------------------------------------------------------!
10        if(allocated(subglbcnc)) deallocate(subglbcnc)

    
 
    end subroutine integrate_xbrick_element
    
    



	subroutine integrate_assemble(elem,K_matrix,F_vector,nofailure)
	!---------------------------------------------------------------------!
    !       integrate and assemble sub element system arrays
    !---------------------------------------------------------------------!     
    	! - passed in variables   
    	type(xbrick_element), intent(inout)	    :: elem
    	real(kind=dp), 	intent(inout)			:: K_matrix(:,:), F_vector(:)
        logical, intent(in)                     :: nofailure
    	! - local variables
    	real(kind=dp),	allocatable           	:: Ki(:,:), Fi(:)   ! sub_elem K matrix and F vector
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
    type(xbrick_element), intent(inout) :: elem


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
    real(dp) :: theta
    
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

        theta=zero


!-----------------------------------------------------------------------!
!               EXTRACTION INTERFACE
!           extract variables from global libraries
!-----------------------------------------------------------------------!
        ! extract elem status variable
        elstat=elem%curr_status

        ! extract edge status variables from glb edge library
        edgstat(:)=lib_edge(elem%edgecnc(:))
        
        ! extract nodal coords from glb node library
        do i=1, nnode
            call extract(lib_node(elem%nodecnc(i)),x=coord(i)%array)
        end do

        ! extract material orientation (fibre angle)
        theta=elem%plyangle

        ! extract elem failed edges' indices
        ifedg=elem%ifailedge

!-----------------------------------------------------------------------!
!       procedure calculations (pure)
!-----------------------------------------------------------------------!
    
!       find and store the broken edges' variables
        if(elstat==intact) then
            do i=1,nedge
                if(edgstat(i)/=intact) then
                    nfailedge=nfailedge+1   ! update total no. of damaged edges
                    ifedg(nfailedge)=i  ! update the indices of damaged edges
                end if
            end do
        else
        ! once partitioned, only check the status of stored failed edges; newly failed edges are ignored
            do i=1,size(ifedg)
                if(ifedg(i)>0) nfailedge=nfailedge+1
            end do
        end if

!       calculate elstat value from edge status variables

!****** if elem is intact, then sort out the most damaged edges first, before updating edge status and doing partitions later
!****** and elem is always partitioned from the most damaged edge, then find the other broken edge along the crack direction
!****** and stores the damaged edge indices of these two edges only; the other damaged edges are ignored; once elem is nolonger
!****** intact (already partitioned along these two edges), no need to go through this sorting again

        if(elstat==intact) then

            if(nfailedge==0) then
                ! no edge failed/damaged, do nothing
                continue
                
            else if(nfailedge==2) then
            ! one pair of edges broken, could be trans elem, not could repartition into ref, tip, wake and failm (matrix failed) elem
                
                ! the failed edge indices must be one pair of corresponding lower and upper edges
                if(.not.((ifedg(2)-ifedg(1))==nedge/2)) then
                    write(msg_file,*)'wrong failed edge indices for nfailedge=2 in xbrick edge status partition'
                    call exit_function
                end if
                    
            else if(nfailedge==4) then
            ! two pairs of edges broken, could be failm (matrix cracked), wake, tip, refinement elem
            
                ! the failed edge indices must be two pairs, with 1 and 2 being the failed lower edges, and 
                ! 3 and 4 being the corresponding failed upper edges
                if(.not.((ifedg(3)-ifedg(1))==nedge/2 .and. (ifedg(4)-ifedg(2))==nedge/2)) then
                    write(msg_file,*)'wrong failed edge indices for nfailedge=4 in xbrick edge status partition'
                    call exit_function
                end if
            
                !**** the partitioning of xbrick elem is entirely based on the bottom four edges ****    
         
                ! if elem partition is intact, then partition accord. to the most damaged edge
              
                ! sort two edges according to their damage severity
                if(edgstat(ifedg(2))>edgstat(ifedg(1))) then
                    jbe1=ifedg(2)
                    jbe2=ifedg(1)
                else
                    jbe1=ifedg(1)
                    jbe2=ifedg(2)
                end if            
                nfailedge=2
                ifedg=0
                ifedg(1)=jbe1
                ifedg(2)=jbe1+nedge/2

                
                ! update of elem status done afterwards
       
            else if( nfailedge==6 ) then
            
                if (.not.((ifedg(4)-ifedg(1))==nedge/2 .and. (ifedg(5)-ifedg(2))==nedge/2 .and. &
                & (ifedg(6)-ifedg(3))==nedge/2)) then
                    write(msg_file,*)'wrong failed edge indices for nfailedge=6 in xbrick edge status partition'
                    call exit_function
                end if
                
                ! if elem partition is intact, then partition accord. to the most damaged edge

                ! find the most damaged edge
                jbe1=ifedg(1)
                do i=2,3
                    if(edgstat(ifedg(i))>edgstat(jbe1)) then
                        jbe1=ifedg(i)
                    end if
                end do
                nfailedge=2
                ifedg=0
                ifedg(1)=jbe1
                ifedg(2)=jbe1+nedge/2


            else if( nfailedge==8 ) then
            
                if(.not. ((ifedg(5)-ifedg(1))==nedge/2 .and. (ifedg(6)-ifedg(2))==nedge/2 .and. &
                & (ifedg(7)-ifedg(3))==nedge/2 .and. (ifedg(8)-ifedg(4))==nedge/2 )) then
                    write(msg_file,*)'wrong failed edge indices for nfailedge=8 in xbrick edge status partition'
                    call exit_function
                end if

                ! if elem partition is intact, then partition accord. to the most damaged edge
                
                ! find the most damaged edge
                jbe1=ifedg(1)
                do i=2,4
                    if(edgstat(ifedg(i))>edgstat(jbe1)) then
                        jbe1=ifedg(i)
                    end if
                end do   
                nfailedge=2
                ifedg=0
                ifedg(1)=jbe1
                ifedg(2)=jbe1+nedge/2
                   
              
            else
            
                write(msg_file,*)'unsupported edge status partition in xbrick!'
                write(msg_file,*) nfailedge, ifedg
                call exit_function
              
            end if

        end if
        
        
!****** if elem is intact, after the above sorting/updating, nfailedge should be 0 or 2 only and they should be a pair of upper/lower edges
!****** if elem is already partitioned, then nfailedge can be 2 or 4 only, and they should be two pairs of upper/lower edges
        if(nfailedge==0 .and. elstat==intact) then
        ! remains intact, do nothing
            continue
        else if(nfailedge==2 .and. (ifedg(2)-ifedg(1))==nedge/2) then 
        ! could be 1st time wake elm, tip elm, ref elm and trans elm
        ! update the edge status, crack tip coords and elstat accordingly
        
            if(edgstat(ifedg(1))==egtrans) then
              ! edge marks the refinement end and the trans elem start 
              ! elem is a trans elem, only this edge needs to be partitioned
                elstat=eltrans
                
            else
              ! another edge must be partitioned to form a ref/tip/wake elem
              ! find the other edge to be partitioned
                
                ! first, find the index of the lower broken edge
                jbe1=ifedg(1)
                
                ! move the corresponding upper broken edge to ifedg(3), s.t. ifedg(2) is free to take the 2nd broken lower edge
                ifedg(3)=ifedg(2)
                
                
                ! find the first (or the second also can) fl. node on the bottom broken edge 
                jnode=topo(3,jbe1)
                
                ! store x,y values of this node in xp1, yp1 (legacy format)
                xp1=coord(jnode)%array(1)
                yp1=coord(jnode)%array(2)

                ! from theta (local fibre dir.), calculate another point on the crack line (could be any point along the line)
                xp2=xp1+cos(theta/halfcirc*pi)
                yp2=yp1+sin(theta/halfcirc*pi)
                
                ! next, find the other edge crossed by the crack line
                do i=1,nedge/2
                    if (i==jbe1) cycle ! the already broken edge, go to next edge
                    
                    ! extract end node 1 coords of edge i
                    jnode=topo(1,i)  
                    x1=coord(jnode)%array(1)
                    y1=coord(jnode)%array(2)
                    z1=coord(jnode)%array(3)
                    
                    ! extract end node 2 coords of edge i
                    jnode=topo(2,i)
                    x2=coord(jnode)%array(1)
                    y2=coord(jnode)%array(2)
                    z2=coord(jnode)%array(3)

                    ! zero cross status and intersection coords for reuse
                    iscross=0
                    xct=zero
                    yct=zero
                    
                    ! check intersection of the crack line and the edge
                    call klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct)
                    if (iscross>0) then
                        nfailedge=nfailedge+1
                        ifedg(2)=i          ! index of the 2nd lower broken edge is i
                        zct=half*(z1+z2)    ! z1 should be the same as z2
                        coord(topo(3,i))%array=[xct,yct,zct]
                        coord(topo(4,i))%array=[xct,yct,zct]
                    end if
                    if (nfailedge==3) exit ! found the edge, no need to proceed
                end do
                
                ! badly shaped element may have large angles; e.g.: 3 or all edges almost parallel to crack, then numerical error may
                ! prevent the algorithm from finding any broken edge or only one broken edge

                if(nfailedge == 2) then
                ! use trial lines: connecting the existing crack tip (xp1,yp1) to the midpoints of the other 3 edges
                    ! crack line equation constants
                    a2=sin(theta/halfcirc*pi)
                    b2=-cos(theta/halfcirc*pi)
                    ! find the midpoint which forms the most parrallel-to-crack line with (xp1,yp1)
                    do i=1,nedge/2 
                        if (i==jbe1) cycle
                        ! tip coords of edge i
                        x1=coord(topo(1,i))%array(1)
                        y1=coord(topo(1,i))%array(2)
                        z1=coord(topo(1,i))%array(3)
                        x2=coord(topo(2,i))%array(1)
                        y2=coord(topo(2,i))%array(2)
                        z2=coord(topo(2,i))%array(3)
                        ! mid point of edge i
                        xmid=half*(x1+x2)
                        ymid=half*(y1+y2)
                        ! line equation constants of midpoint-(xp1,yp1)
                        a1=ymid-yp1
                        b1=xp1-xmid
                        ! initialize detlc and intersection info
                        if(i==1 .or. (jbe1==1 .and. i==2)) then
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
                        end if
                    end do
                    coord(topo(3,ifedg(2)))%array=[xct,yct,zct]
                    coord(topo(4,ifedg(2)))%array=[xct,yct,zct]
                end if
                

                jbe2=ifedg(2)   ! store it in jbe2 (safer)

                ! update edge status variables
                if(edgstat(jbe1)==egref) then ! tip elem end, refinement elem. start (1st time)
                    elstat=elref               
                    ! change 2nd broken edge status to 1 (trans elem start)
                    edgstat(jbe2)=egtrans
             
                else if(edgstat(jbe1)==egtip) then ! wake elem end, tip elem start (1st time)
                    elstat=eltip
                    ! change 2nd broken edge status to 2 (refinement start)
                    edgstat(jbe2)=egref
               
                else if(edgstat(jbe1)>=cohcrack) then ! wake elem start (1st time)
                    elstat=elwake !wake elem               
                    ! change 2nd broken edge status to 3 (tip elem start)
                    edgstat(jbe2)=egtip
                    
                else ! unknown edge status
                    write(msg_file,*)'unknown edge status!'
                    call exit_function
                end if

                !***** project var. values of the bottom broken edges to the top ones *****
                
                nfailedge=4            ! no. of broken edge is now 4
                ifedg(4)=jbe2+nedge/2  ! index of the 4th broken edge is i
                
                ! store it in jbe4 (safer)
                jbe4=ifedg(4)
                
                ! update the two fl. node coords on this edge
                jnode=topo(1,jbe4); z1=coord(jnode)%array(3)
                jnode=topo(2,jbe4); z2=coord(jnode)%array(3)
                zct=half*(z1+z2) ! z1 should be the same as z2
                
                jnode=topo(3,jbe4)
                coord(jnode)%array=[xct,yct,zct]
                jnode=topo(4,jbe4)               
                coord(jnode)%array=[xct,yct,zct]

                ! update edge status var. value
                edgstat(jbe4)=edgstat(jbe2)
        
            endif
            
        else if(nfailedge==4 .and. (ifedg(3)-ifedg(1))==nedge/2 .and. (ifedg(4)-ifedg(2))==nedge/2) then
        ! could be cracked, wake, tip, refinement elem
        ! only update the elstat, not the edge status, nor the crack tip coords
        
            !**** the partitioning of xbrick elem is entirely based on the bottom four edges ****         
            
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
                elstat=elfailm
            else ! unknown combination
                write(msg_file,*)'unknown combination of 2 edge status!'
                call exit_function
            end if
            
        else
            write(msg_file,*)'unsupported nfailedge value for edge and el stat update in xbrick edge stat partition!'
            call exit_function 
        end if     

!-----------------------------------------------------------------------!
!                   UPDATE INTERFACE
!               update global libraries
!-----------------------------------------------------------------------!
            
!       update element curr_status and sub-element cnc matrices

        if(elstat>elem%curr_status) then

            elem%curr_status=elstat 
            elem%ifailedge=ifedg                      

            call update_subcnc(elem,edgstat,ifedg,nfailedge)
        

            ! update glb edge array, only the broken edges' status variables
            lib_edge(elem%edgecnc(ifedg(1:nfailedge)))=edgstat(ifedg(1:nfailedge))
            
            ! update glb node array, only the broken edges' fl. node coord
            do i=1, nfailedge
                j=ifedg(i)          ! local index of broken edge
                
                l=topo(3,j)         ! elem lcl index of broken edge fl. node 1
                k=elem%nodecnc(l)   ! global index of broken edge fl. node 1
                call update(lib_node(k),x=coord(l)%array)
                
                l=topo(4,j)         ! elem lcl index of broken edge fl. node 2
                k=elem%nodecnc(l)   ! global index of broken edge fl. node 2
                call update(lib_node(k),x=coord(l)%array)
            end do
        end if





!       deallocate local dynamic arrays

    
    
    end subroutine edge_status_partition
    
  





!********************************************************************************************
!******************* subroutine kplyfail ****************************************************
!*********** quadratic stress failure criteria and element partition criterion **************
!********************************************************************************************
    subroutine failure_criterion_partition(elem)
! this subroutine update elem status & partition to elfailm if any sub elem is nolonger intact

    ! passed in variables
    type(xbrick_element) :: elem
  
    ! extracted variables, from glb libraries
    integer :: edgstat(nedge)               ! status variable array of element edges
    type(real_alloc_array) :: coord(nnode)  ! nodal coord arrays to store the coords of elem nodes extracted from glb node lib
    
    real(dp)    :: xelm(ndim,nnode)         ! local copy of element nodal coords
    
    real(dp)    :: theta
    
    
    
    ! local variable
    
    integer :: nfailedge        ! no. of failed edges in the element
    integer :: ifedg(nedge)     ! index of failed edges in the element
    integer :: elstat           ! local copy of elem curr status
    integer :: jbe1,jbe2, jnode ! indices of 2 broken edges, and a variable to hold a glb node index
    integer :: iscross          ! indicator of line intersection; >0 if two lines intersect
    integer :: i, j, l, k       ! counters
    integer :: subelstat        ! sub elem status
    logical :: failed           ! true if elem is failed (one of the sub elems reached failure onset)
    logical :: fbfail           ! true if elem has fibre failure (intact elem reached fibre failure onset)
    
      
    real(dp) :: xp1, yp1, xp2, yp2      ! (x,y) of point 1 and point 2 on the crack line
    real(dp) :: x1, y1, z1, x2, y2, z2  ! (x,y) of node 1 and node 2 of an element edge
    real(dp) :: xct, yct, zct           ! (x,y) of a crack tip on an edge, i.e., intersection of crack line & edge
    real(dp) :: xo, yo                  ! (x,y) of element centroid
    real(dp) :: detlc,a1,b1,a2,b2,xmid,ymid
    real(dp) :: xmid1,ymid1,zmid1,xmid2,ymid2,zmid2
    real(dp) :: xct1,yct1,zct1,xct2,yct2,zct2
    
    
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
    
        ! subroutine only partitions elstat < elfailm elems
        if(elem%curr_status>=elfailm) return ! elem already failed, no need to proceed
    
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
        edgstat(:)=lib_edge(elem%edgecnc(:))
        
        ! extract nodal coords from glb node library
        do i=1, nnode
            call extract(lib_node(elem%nodecnc(i)),x=coord(i)%array)
            if(allocated(coord(i)%array)) xelm(:,i)=coord(i)%array(:)
        end do
        

        ! extract material orientation (fibre angle)
        theta=elem%plyangle
        
        elstat=elem%curr_status

        ifedg=elem%ifailedge


!-----------------------------------------------------------------------!
!           check failure criterion on all sub elements
!   elem is judged failed a.l.a. 1 sub elem has reached failure onset
!   elem is judged fibre failed if elstat=intact and subelstat=fibre_onset
!-----------------------------------------------------------------------!

        do i=1, size(elem%subelem)
            call extract(elem%subelem(i),curr_status=subelstat)
            if(subelstat>intact) failed=.true.
            if(failed) exit
        end do

        if (elstat==intact .and. subelstat>=fibre_onset) then
            fbfail=.true.
        end if
!-----------------------------------------------------------------------!

       
            
        if(failed) then
        ! element failed, find newly broken edge and update edge status vars
            
        
          ! find xp1, yp1
            if(elstat .eq. intact) then ! element originally intact
                ! find centroid
                xo=quarter*(xelm(1,1)+xelm(1,2)+xelm(1,3)+xelm(1,4))
                yo=quarter*(xelm(2,1)+xelm(2,2)+xelm(2,3)+xelm(2,4))
                xp1=xo
                yp1=yo
                ! find xp2, yp2
                xp2=xp1+cos(theta/halfcirc*pi)
                yp2=yp1+sin(theta/halfcirc*pi)
                do i=1,nedge/2 
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
                    do i=1,nedge/2-1
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
                        do j=i+1,nedge/2
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
                    do i=1,nedge/2 
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
                ifedg(3)=ifedg(1)+nedge/2
                ifedg(4)=ifedg(2)+nedge/2
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
                    do i=1,nedge/2
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
                        do i=1,nedge/2 
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
                    ifedg(4)=ifedg(3)+nedge/2
                    edgstat(ifedg(4))=cohcrack
                    xelm(1:2,topo(3,ifedg(4)))=xelm(1:2,topo(3,ifedg(3)))
                    xelm(1:2,topo(4,ifedg(4)))=xelm(1:2,topo(4,ifedg(3)))
                    xelm(3,topo(3,ifedg(4)))=half*(xelm(3,topo(1,ifedg(4)))+xelm(3,topo(2,ifedg(4))))
                    xelm(3,topo(4,ifedg(4)))=half*(xelm(3,topo(1,ifedg(4)))+xelm(3,topo(2,ifedg(4))))
                    
                else if (elstat.gt.eltrans .and. elstat.lt.elfailm) then 
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
            ! if intact elem reaches fibre failure onset, new elem partition is named elfailf
                elstat=elfailf
            else
                elstat=elfailm
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
            lib_edge(elem%edgecnc(ifedg(1:nfailedge)))=edgstat(ifedg(1:nfailedge))
            
            ! update glb node array, only the broken edges' fl. node coord
            do i=1, nfailedge
                j=ifedg(i)          ! local index of broken edge
                
                l=topo(3,j)          ! elem lcl index of broken edge fl. node 1
                coord(l)%array(:)=xelm(:,l)
                k=elem%nodecnc(l)   ! global index of broken edge fl. node 1
                call update(lib_node(k),x=coord(l)%array)
                
                l=topo(4,j)          ! elem lcl index of broken edge fl. node 2
                coord(l)%array(:)=xelm(:,l)
                k=elem%nodecnc(l)   ! global index of broken edge fl. node 2
                call update(lib_node(k),x=coord(l)%array)
            end do



        end if

        ! deallocate local arrays
        
        
      return  
      end subroutine failure_criterion_partition






       
    subroutine update_subcnc(elem,edgstat,ifedg,nfailedge)
    
    ! passed-in variables
    type(xbrick_element),    intent(inout)   :: elem
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
            if(edgstat(ibe)/=egtrans) then
                write(msg_file,*)'transition partition only accepts edgstat=egtrans!'
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
            if(elem%nodecnc(topo(3,ibe))<elem%nodecnc(topo(4,ibe))) then
                jnode=topo(3,ibe)
            else
                jnode=topo(4,ibe)
            end if
            
            ! sub elm 1 connec
            elem%subcnc(1)%array(1)=topo(1,e1)
            elem%subcnc(1)%array(2)=topo(2,e1)
            elem%subcnc(1)%array(3)=jnode
            
            elem%subcnc(1)%array(4:5)=elem%subcnc(1)%array(1:2)+nndrl/2 ! upper surf nodes
            elem%subcnc(1)%array(6)=elem%subcnc(1)%array(3)+nndfl/2
            
            subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
            
            ! sub elm 2 connec
            elem%subcnc(2)%array(1)=topo(1,e2)
            elem%subcnc(2)%array(2)=topo(2,e2)
            elem%subcnc(2)%array(3)=jnode 

            elem%subcnc(2)%array(4:5)=elem%subcnc(2)%array(1:2)+nndrl/2 ! upper surf nodes
            elem%subcnc(2)%array(6)=elem%subcnc(2)%array(3)+nndfl/2 

            subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))
            
            ! sub elm 3 connec
            elem%subcnc(3)%array(1)=topo(1,e3)
            elem%subcnc(3)%array(2)=topo(2,e3)
            elem%subcnc(3)%array(3)=jnode  
            
            elem%subcnc(3)%array(4:5)=elem%subcnc(3)%array(1:2)+nndrl/2 ! upper surf nodes
            elem%subcnc(3)%array(6)=elem%subcnc(3)%array(3)+nndfl/2

            subglbcnc(3)%array(:)=elem%nodecnc(elem%subcnc(3)%array(:))

            ! create sub elements
            call prepare(elem%subelem(1),eltype='wedge', matkey=elem%bulkmat, plyangle=elem%plyangle, glbcnc=subglbcnc(1)%array)
            call prepare(elem%subelem(2),eltype='wedge', matkey=elem%bulkmat, plyangle=elem%plyangle, glbcnc=subglbcnc(2)%array)
            call prepare(elem%subelem(3),eltype='wedge', matkey=elem%bulkmat, plyangle=elem%plyangle, glbcnc=subglbcnc(3)%array)
         

 
        case (4) !- two edges cracked
         
            ibe1=min(ifedg(1),ifedg(2),ifedg(3),ifedg(4))   ! local edge index of 1st broken edge
            do i=1, 4
                if(ifedg(i)<=nedge/2) ibe2=max(ibe2,ifedg(i))! local edge index of 2nd broken edge
            end do  
            
            

            ! ibe1 must be between 1 to 3, and ibe2 between 2 to 4, with ibe2 > ibe1
            if(ibe1>=nedge/2 .or. ibe1==0 .or. ibe2<=1 .or. ibe2<=ibe1) then
                write(msg_file,*) 'something wrong in xbrick update subcnc case nfailedge=4'
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
                            write(msg_file,*)'wrong 2nd broken edge in update subcnc xbrick'
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
                            write(msg_file,*)'wrong 2nd broken edge in update subcnc xbrick'
                            call exit_function
                    end select
                case(3)
                    if(ibe2==4) then
                        nbulk=4
                        e1=3; e2=4; e3=1; e4=2
                    else
                        write(msg_file,*)'wrong 2nd broken edge in update subcnc xbrick'
                        call exit_function                    
                    end if    
                case default
                    write(msg_file,*)'wrong broken edge in update subcnc xbrick'
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
                    if(elem%nodecnc(topo(3,e1))<elem%nodecnc(topo(4,e1))) then
                        jnode1=topo(3,e1)
                    else
                        jnode1=topo(4,e1)
                    end if
                    
                    ! find the smaller glb fl. node on the broken edge
                    if(elem%nodecnc(topo(3,e3))<elem%nodecnc(topo(4,e3))) then
                        jnode3=topo(3,e3)
                    else
                        jnode3=topo(4,e3)
                    end if
                    
                    ! sub elm 1 connec
                    elem%subcnc(1)%array(1)=topo(1,e1)
                    elem%subcnc(1)%array(2)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(1)%array(2)=jnode1
                    elem%subcnc(1)%array(3)=topo(4,e3); if(edgstat(e3)<cohcrack) elem%subcnc(1)%array(3)=jnode3
                    elem%subcnc(1)%array(4)=topo(2,e3)
                    
                    elem%subcnc(1)%array(5)=elem%subcnc(1)%array(1)+nndrl/2
                    elem%subcnc(1)%array(6)=elem%subcnc(1)%array(2)+nndfl/2
                    elem%subcnc(1)%array(7)=elem%subcnc(1)%array(3)+nndfl/2
                    elem%subcnc(1)%array(8)=elem%subcnc(1)%array(4)+nndrl/2
                    
                    subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
                    
                    ! sub elm 2 connec
                    elem%subcnc(2)%array(1)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(2)%array(1)=jnode1
                    elem%subcnc(2)%array(2)=topo(2,e1)
                    elem%subcnc(2)%array(3)=topo(1,e3)
                    elem%subcnc(2)%array(4)=topo(3,e3); if(edgstat(e3)<cohcrack) elem%subcnc(2)%array(4)=jnode3
                    
                    elem%subcnc(2)%array(5)=elem%subcnc(2)%array(1)+nndfl/2
                    elem%subcnc(2)%array(6)=elem%subcnc(2)%array(2)+nndrl/2
                    elem%subcnc(2)%array(7)=elem%subcnc(2)%array(3)+nndrl/2
                    elem%subcnc(2)%array(8)=elem%subcnc(2)%array(4)+nndfl/2

                    subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))
                    
                    ! create sub bulk elements
                    
                    call prepare(elem%subelem(1),eltype='brick',matkey=elem%bulkmat,plyangle=elem%plyangle,&
                    &glbcnc=subglbcnc(1)%array)
                    call prepare(elem%subelem(2),eltype='brick',matkey=elem%bulkmat,plyangle=elem%plyangle,&
                    &glbcnc=subglbcnc(2)%array)
                    
                    
                    if(iscoh) then
                    
                    ! sub elm 3 connec
                    elem%subcnc(3)%array(2)=topo(4,e3); if(edgstat(e3)<cohcrack) elem%subcnc(3)%array(2)=jnode3
                    elem%subcnc(3)%array(1)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(1)=jnode1
                    
                    elem%subcnc(3)%array(3)=elem%subcnc(3)%array(2)+nndfl/2
                    elem%subcnc(3)%array(4)=elem%subcnc(3)%array(1)+nndfl/2
                    
                    elem%subcnc(3)%array(6)=topo(3,e3); if(edgstat(e3)<cohcrack) elem%subcnc(3)%array(6)=jnode3
                    elem%subcnc(3)%array(5)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(5)=jnode1
                    
                    elem%subcnc(3)%array(7)=elem%subcnc(3)%array(6)+nndfl/2
                    elem%subcnc(3)%array(8)=elem%subcnc(3)%array(5)+nndfl/2
                    
                    subglbcnc(3)%array(:)=elem%nodecnc(elem%subcnc(3)%array(:))
                    
                    
                    call prepare(elem%subelem(3),eltype='coh3d8',matkey=elem%cohmat,glbcnc=subglbcnc(3)%array)
                    
                    
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
                    if(elem%nodecnc(topo(3,e1))<elem%nodecnc(topo(4,e1))) then
                        jnode1=topo(3,e1)
                    else
                        jnode1=topo(4,e1)
                    end if
                    
                    ! find the smaller glb fl. node on the broken edge
                    if(elem%nodecnc(topo(3,e2))<elem%nodecnc(topo(4,e2))) then
                        jnode2=topo(3,e2)
                    else
                        jnode2=topo(4,e2)
                    end if
                    
                    ! sub elm 1 connec
                    elem%subcnc(1)%array(1)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(1)%array(1)=jnode1
                    elem%subcnc(1)%array(2)=topo(1,e2)
                    elem%subcnc(1)%array(3)=topo(3,e2); if(edgstat(e2)<cohcrack) elem%subcnc(1)%array(3)=jnode2

                    elem%subcnc(1)%array(4)=elem%subcnc(1)%array(1)+nndfl/2 ! upper surf nodes
                    elem%subcnc(1)%array(5)=elem%subcnc(1)%array(2)+nndrl/2
                    elem%subcnc(1)%array(6)=elem%subcnc(1)%array(3)+nndfl/2
                    
                    subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
                    
                    ! sub elm 2 connec
                    elem%subcnc(2)%array(1)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(2)%array(1)=jnode2
                    elem%subcnc(2)%array(2)=topo(1,e3)
                    elem%subcnc(2)%array(3)=topo(2,e3)
                    
                    elem%subcnc(2)%array(4)=elem%subcnc(2)%array(1)+nndfl/2 ! upper surf nodes
                    elem%subcnc(2)%array(5)=elem%subcnc(2)%array(2)+nndrl/2
                    elem%subcnc(2)%array(6)=elem%subcnc(2)%array(3)+nndrl/2
                    
                    subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))
                    
                    ! sub elm 3 connec
                    elem%subcnc(3)%array(1)=topo(1,e4)
                    elem%subcnc(3)%array(2)=topo(2,e4)
                    elem%subcnc(3)%array(3)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(3)=jnode1
                    
                    elem%subcnc(3)%array(4)=elem%subcnc(3)%array(1)+nndrl/2 ! upper surf nodes
                    elem%subcnc(3)%array(5)=elem%subcnc(3)%array(2)+nndrl/2
                    elem%subcnc(3)%array(6)=elem%subcnc(3)%array(3)+nndfl/2
                    
                    subglbcnc(3)%array(:)=elem%nodecnc(elem%subcnc(3)%array(:))
                    
                    ! sub elm 4 connec
                    elem%subcnc(4)%array(1)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(4)%array(1)=jnode1
                    elem%subcnc(4)%array(2)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(4)%array(2)=jnode2
                    elem%subcnc(4)%array(3)=topo(2,e3)

                    elem%subcnc(4)%array(4)=elem%subcnc(4)%array(1)+nndfl/2 ! upper surf nodes
                    elem%subcnc(4)%array(5)=elem%subcnc(4)%array(2)+nndfl/2
                    elem%subcnc(4)%array(6)=elem%subcnc(4)%array(3)+nndrl/2
                    
                    subglbcnc(4)%array(:)=elem%nodecnc(elem%subcnc(4)%array(:))
                    
                    ! create sub bulk elements
                    call prepare(elem%subelem(1),eltype='wedge',matkey=elem%bulkmat,&
                    &plyangle=elem%plyangle,glbcnc=subglbcnc(1)%array)
                    call prepare(elem%subelem(2),eltype='wedge',matkey=elem%bulkmat,&
                    &plyangle=elem%plyangle,glbcnc=subglbcnc(2)%array)
                    call prepare(elem%subelem(3),eltype='wedge',matkey=elem%bulkmat,&
                    &plyangle=elem%plyangle,glbcnc=subglbcnc(3)%array)
                    call prepare(elem%subelem(4),eltype='wedge',matkey=elem%bulkmat,&
                    &plyangle=elem%plyangle,glbcnc=subglbcnc(4)%array)
                    
                    if(iscoh) then
                    
                    ! sub elm 5 connec
                    elem%subcnc(5)%array(5)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(5)%array(5)=jnode1
                    elem%subcnc(5)%array(6)=topo(3,e2); if(edgstat(e2)<cohcrack) elem%subcnc(5)%array(6)=jnode2
                    
                    elem%subcnc(5)%array(7)=elem%subcnc(5)%array(6)+nndfl/2
                    elem%subcnc(5)%array(8)=elem%subcnc(5)%array(5)+nndfl/2
                    
                    elem%subcnc(5)%array(1)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(5)%array(1)=jnode1
                    elem%subcnc(5)%array(2)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(5)%array(2)=jnode2   
                    
                    elem%subcnc(5)%array(3)=elem%subcnc(5)%array(2)+nndfl/2
                    elem%subcnc(5)%array(4)=elem%subcnc(5)%array(1)+nndfl/2
                    
                    subglbcnc(5)%array(:)=elem%nodecnc(elem%subcnc(5)%array(:))
                    
                    call prepare(elem%subelem(5),eltype='coh3d8', matkey=elem%cohmat, glbcnc=subglbcnc(5)%array)
                    
                    end if
                    
                case default
                    write(msg_file,*)'wrong nbulk in update subcnc xbrick'
                    call exit_function
            end select

          
!        case(6)
        ! do not update partition
            

!        case(8)


           
        case default
           write(msg_file,*) 'WARNING: xbrick update subcnc case selection default!'
           
        end select
        
        
        ! deallocate local array
        
        if(allocated(subglbcnc)) deallocate(subglbcnc)


    end subroutine update_subcnc


  
  
  
  
  
  
  
  
  
  
end module xbrick_element_module
