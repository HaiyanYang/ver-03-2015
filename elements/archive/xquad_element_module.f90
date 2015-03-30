module xquad_element_module
    use parameter_module
    use toolkit_module                  ! global tools for element integration
    use lib_edge_module                 ! global edge library
    use lib_node_module                 ! global node library
    use lib_mat_module                  ! global material library
    use sub2d_element_module


    implicit none
    private

    integer,parameter :: ndim=2, nndrl=4, nedge=4, nndfl=2*nedge, nnode=nndrl+nndfl, ndof=ndim*nnode
    ! Topology: nodes on each edge; 4 nodes per edge, 1-2 are end nodes, 3-4 are fl. nodes; assigned in lcl node numbers
    integer,parameter :: topo(4,nedge)=reshape([1,2,5,6,2,3,7,8,3,4,9,10,4,1,11,12],[4,nedge])
    
    ! element status variable values
    integer, parameter :: eltrans=1, elref=2, eltip=3, elwake=4, elfail=5
    
    ! edge status variable values
    integer, parameter :: egtrans=1, egref=2, egtip=3, wkcrack=3, cohcrack=4, strgcrack=5

    type, public :: xquad_element             ! breakable quadrilateral
        private

        integer :: curr_status=0        ! 0 means intact
        integer :: key=0
        integer :: bulkmat=0
        integer :: cohmat=0

        integer :: nodecnc(nnode)=0     ! cnc to glb node arrays for accessing nodal variables (x, u, du, v, dof ...)
        integer :: edgecnc(nedge)=0     ! cnc to glb edge arrays for accessing edge variables (failure status)

        type(sub2d_element), allocatable :: subelem(:)
        type(int_alloc_array), allocatable :: subcnc(:)      ! sub_elem connec to parent elem nodes

    end type xquad_element

    interface empty
        module procedure empty_xquad_element
    end interface

    interface prepare
        module procedure prepare_xquad_element
    end interface

    !~interface precrack
    !~    module procedure precrack_xquad_element
    !~end interface

    interface integrate
        module procedure integrate_xquad_element
    end interface

    interface extract
        module procedure extract_xquad_element
    end interface




    public :: empty,prepare,integrate,extract



    contains




    ! empty a breakable quadrilateral
    subroutine empty_xquad_element(elem)

        type(xquad_element),intent(out) :: elem

        elem%curr_status=0
        elem%key=0
        elem%bulkmat=0
        elem%cohmat=0

        elem%nodecnc=0
        elem%edgecnc=0

        if(allocated(elem%subelem)) deallocate(elem%subelem)
        if(allocated(elem%subcnc))  deallocate(elem%subcnc)

    end subroutine empty_xquad_element



    ! this subroutine is used to prepare the connectivity and material lib index of the element
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine prepare_xquad_element(elem,key,bulkmat,cohmat,nodecnc,edgecnc)

        type(xquad_element),    intent(inout)   :: elem
        integer,                intent(in)      :: key
        integer,                intent(in)      :: bulkmat, cohmat
        integer,                intent(in)      :: nodecnc(nnode)
        integer,                intent(in)      :: edgecnc(nedge)

        elem%key=key
        elem%bulkmat=bulkmat
        elem%cohmat=cohmat
        elem%nodecnc=nodecnc
        elem%edgecnc=edgecnc

    end subroutine prepare_xquad_element


    subroutine extract_xquad_element(elem,curr_status,key,bulkmat,cohmat,nodecnc,edgecnc,subelem,subcnc)

        type(xquad_element),                      intent(in)  :: elem
        integer,                        optional, intent(out) :: curr_status
        integer,                        optional, intent(out) :: key
        integer,                        optional, intent(out) :: bulkmat
        integer,                        optional, intent(out) :: cohmat
        integer,            allocatable,optional, intent(out) :: nodecnc(:)
        integer,            allocatable,optional, intent(out) :: edgecnc(:)
        type(sub2d_element),allocatable,optional, intent(out) :: subelem(:)
        type(int_alloc_array),allocatable,optional,intent(out):: subcnc(:)

        if(present(curr_status)) curr_status=elem%curr_status
        if(present(key)) key=elem%key
        if(present(bulkmat)) bulkmat=elem%bulkmat
        if(present(cohmat)) cohmat=elem%cohmat

        if(present(nodecnc)) then
            allocate(nodecnc(nnode))
            nodecnc=elem%nodecnc
        end if

        if(present(edgecnc)) then
            allocate(edgecnc(nedge))
            edgecnc=elem%edgecnc
        end if

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

    end subroutine extract_xquad_element






    subroutine integrate_xquad_element(elem, K_matrix, F_vector)

        type(xquad_element),intent(inout)       :: elem
        real(kind=dp),allocatable,intent(out)   :: K_matrix(:,:), F_vector(:)


        ! local variables
        type(int_alloc_array), allocatable  :: subglbcnc(:)     ! glb cnc of sub element, used when elem is intact
        real(kind=dp),allocatable           :: Ki(:,:), Fi(:)   ! sub_elem K matrix and F vector

        integer :: i,j,l, elstat
        integer, allocatable :: dofcnc(:)


        ! initialize K & F
        allocate(K_matrix(ndof,ndof),F_vector(ndof))
        K_matrix=zero; F_vector=zero

        ! initialize local variables
        i=0; j=0; l=0
        elstat=0



        !---------------------------------------------------------------------!
        !               update sub element definitions
        !---------------------------------------------------------------------!

        ! if elem is not yet failed, check elem edge status variables and update elem status and sub elem cnc
        if(elem%curr_status<elfail) then
            call edge_status_partition(elem)
        end if

        !****** after edge status update, elem curr status can be any value from intact to failed, but intact-case subelem hasn't been created


        ! if elem is still intact after checking edge status (no broken edges), assign 1 quad subelem if not yet done
        if(elem%curr_status==intact) then
            if(.not.allocated(elem%subelem)) then
                allocate(elem%subelem(1))
                allocate(elem%subcnc(1))
                allocate(elem%subcnc(1)%array(4))   ! quad elem
                allocate(subglbcnc(1))
                allocate(subglbcnc(1)%array(4))
                ! sub elm 1 connec
                elem%subcnc(1)%array=[1,2,3,4]
                subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))
                ! create sub elements
                call prepare(elem%subelem(1),eltype='quad', matkey=elem%bulkmat, glbcnc=subglbcnc(1)%array)
            end if
        end if

        !******* reaching here, elem curr status can be any value from intact to failed, and in all cases, subelems have been created




        !---------------------------------------------------------------------!
        !       integrate and assemble sub element system arrays
        !---------------------------------------------------------------------!

        ! if elem is not yet failed, integrate sub elem and check the failure criterion, and repartition if necessary
        if(elem%curr_status<elfail) then
            ! store current status value
            elstat=elem%curr_status

            ! integrate sub elements and assemble into global matrix
            do i=1, size(elem%subelem)
                call integrate(elem%subelem(i),Ki,Fi)
                if(allocated(dofcnc)) deallocate(dofcnc)
                allocate(dofcnc(size(Fi))); dofcnc=0
                do j=1, size(elem%subcnc(i)%array) ! no. of nodes in sub elem i
                    do l=1, ndim
                        ! dof indices of the jth node of sub elem i 
                        dofcnc((j-1)*ndim+l)=(elem%subcnc(i)%array(j)-1)*ndim+l
                    end do
                end do
                call assembleKF(K_matrix,F_vector,Ki,Fi,dofcnc)
            end do

            !***** check failure criterion *****
            call failure_criterion_partition(elem)

            if(elem%curr_status/=elstat) then
            ! elem status changed, elem failed, partition changed, reintegrate subelems
                ! empty K and F for reuse
                K_matrix=zero; F_vector=zero
                ! integrate sub elements and assemble into global matrix
                do i=1, size(elem%subelem)
                    call integrate(elem%subelem(i),Ki,Fi)
                    if(allocated(dofcnc)) deallocate(dofcnc)
                    allocate(dofcnc(size(Fi))); dofcnc=0
                    do j=1, size(elem%subcnc(i)%array) ! no. of nodes in sub elem i
                        do l=1, ndim
                            ! dof indices of the jth node of sub elem i 
                            dofcnc((j-1)*ndim+l)=(elem%subcnc(i)%array(j)-1)*ndim+l
                        end do
                    end do
                    call assembleKF(K_matrix,F_vector,Ki,Fi,dofcnc)
                end do
            end if

        else if(elem%curr_status==elfail) then
        ! element is already failed, integrate and assemble subelem

            ! integrate sub elements and assemble into global matrix
            do i=1, size(elem%subelem)
                call integrate(elem%subelem(i),Ki,Fi)
                if(allocated(dofcnc)) deallocate(dofcnc)
                allocate(dofcnc(size(Fi))); dofcnc=0
                do j=1, size(elem%subcnc(i)%array) ! no. of nodes in sub elem i
                    do l=1, ndim
                        ! dof indices of the jth node of sub elem i 
                        dofcnc((j-1)*ndim+l)=(elem%subcnc(i)%array(j)-1)*ndim+l
                    end do
                end do
                call assembleKF(K_matrix,F_vector,Ki,Fi,dofcnc)
            end do

        else
            write(msg_file,*)'unsupported elem curr status value in xquad!'
            call exit_function
        end if




        !---------------------------------------------------------------------!
        !               deallocate local arrays
        !---------------------------------------------------------------------!
        if(allocated(Ki)) deallocate(Ki)
        if(allocated(Fi)) deallocate(Fi)
        if(allocated(subglbcnc)) deallocate(subglbcnc)
        if(allocated(dofcnc)) deallocate(dofcnc)


    end subroutine integrate_xquad_element







    subroutine edge_status_partition(elem)


    ! passed-in variables
    type(xquad_element), intent(inout) :: elem


    ! extracted variables, from glb libraries
    integer :: edgstat(nedge)               ! status variable array of element edges
    type(real_alloc_array) :: coord(nnode)  ! nodal coord arrays to store the coords of elem nodes extracted from glb node lib



    ! local variable

    integer :: nfailedge        ! no. of failed edges in the element
    integer :: ifedg(nedge)     ! index of failed edges in the element
    integer :: elstat           ! local copy of elem curr status
    integer :: jbe1,jbe2, jnode ! indices of 2 broken edges, and a variable to hold a glb node index
    integer :: iscross          ! indicator of line intersection; >0 if two lines intersect
    integer :: i, j, l, k       ! counters

    real(dp) :: xp1, yp1, xp2, yp2      ! (x,y) of point 1 and point 2 on the crack line
    real(dp) :: x1, y1, x2, y2          ! (x,y) of node 1 and node 2 of an element edge
    real(dp) :: xct, yct                ! (x,y) of a crack tip on an edge, i.e., intersection of crack line & edge
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
        jbe1=0; jbe2=0; jnode=0
        iscross=0
        i=0; j=0; l=0; k=0

        xp1=zero; yp1=zero
        xp2=zero; yp2=zero

        x1=zero; y1=zero
        x2=zero; y2=zero

        xct=zero; yct=zero
        theta=zero


!-----------------------------------------------------------------------!
!               EXTRACTION INTERFACE
!           extract variables from global libraries
!-----------------------------------------------------------------------!
        ! extract edge status variables from glb edge library
        edgstat(:)=lib_edge(elem%edgecnc(:))

        ! extract nodal coords from glb node library
        do i=1, nnode
            call extract(lib_node(elem%nodecnc(i)),x=coord(i)%array)
        end do
        
        ! extract material orientation (fibre angle)
        call extract(lib_mat(elem%bulkmat),theta=theta)



!-----------------------------------------------------------------------!
!       procedure calculations (pure)
!-----------------------------------------------------------------------!

!       find and store the broken edges' variables
        do i=1,nedge
            if(edgstat(i)/=intact) then
                nfailedge=nfailedge+1   ! update total no. of damaged edges
                ifedg(nfailedge)=i  ! update the indices of damaged edges
            end if
        end do


!       calculate elstat value from edge status variables

10      if(nfailedge==0) then
            ! no edge failed/damaged, do nothing
            continue

        else if(nfailedge==1) then
        ! could be 1st time wake elm, tip elm, ref elm and trans elm
        ! update the edge status, crack tip coords and elstat accordingly

            if(edgstat(ifedg(1))==egtrans) then
              ! edge marks the refinement end and the trans elem start
              ! elem is a trans elem, only this edge needs to be partitioned
                elstat=eltrans

            else
              ! another edge must be partitioned to form a ref/tip/wake elem
              ! find the other edge to be partitioned

                ! first, find the index of the broken edge
                jbe1=ifedg(1)

                ! find the first (or the second also can) fl. node on this edge
                jnode=topo(3,jbe1)

                ! store x,y values of this node in xp1, yp1 (legacy format)
                xp1=coord(jnode)%array(1)
                yp1=coord(jnode)%array(2)

                ! from theta (local fibre dir.), calculate another point on the crack line (could be any point along the line)
                xp2=xp1+cos(theta/halfcirc*pi)
                yp2=yp1+sin(theta/halfcirc*pi)

                ! next, find the other edge crossed by the crack line
                do i=1,nedge
                    if (i==jbe1) cycle ! the already broken edge, go to next edge

                    ! extract end node 1 coords of edge i
                    jnode=topo(1,i)
                    x1=coord(jnode)%array(1)
                    y1=coord(jnode)%array(2)

                    ! extract end node 2 coords of edge i
                    jnode=topo(2,i)
                    x2=coord(jnode)%array(1)
                    y2=coord(jnode)%array(2)

                    ! zero cross status and intersection coords for reuse
                    iscross=0
                    xct=zero
                    yct=zero

                    ! check intersection of the crack line and the edge
                    call klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct)
                    if (iscross>0) exit ! found the edge, no need to proceed
                end do

                ! update nfailedge & ifedg
                nfailedge=2     ! no. of broken edge is now 2
                ifedg(2)=i  ! index of the 2nd broken edge is i

                ! store it in jbe2 (safer)
                jbe2=i

                ! update the two fl. node coords on this edge
                jnode=topo(3,jbe2)
                coord(jnode)%array=[xct,yct]
                jnode=topo(4,jbe2)
                coord(jnode)%array=[xct,yct]

                ! update edge status variables
                if(edgstat(jbe1)==egref) then ! tip elem end, refinement elem. start (1st time)
                    elstat=elref
                    ! change 2nd broken edge status to 1 (trans elem start)
                    edgstat(jbe2)=egtrans

                else if(edgstat(jbe1)==3) then ! wake elem end, tip elem start (1st time)
                    elstat=eltip
                    ! change 2nd broken edge status to 2 (refinement start)
                    edgstat(jbe2)=egref

                else if(edgstat(jbe1)>=4) then ! wake elem start (1st time)
                    elstat=elwake !wake elem
                    ! change 2nd broken edge status to 3 (tip elem start)
                    edgstat(jbe2)=egtip

                else ! unknown edge status
                    write(msg_file,*)'unknown edge status!'
                    call exit_function
                end if
            endif

        else if(nfailedge==2) then
        ! could be cracked, wake, tip, refinement elem
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
                elstat=elfail
            else ! unknown combination
                write(msg_file,*)'unknown combination of 2 edge status!'
                call exit_function
            end if

        else if(nfailedge>2) then
        ! not yet supported

            write(msg_file,*)'more than 2 broken edges not yet supported!'

            ! partition only according to the first two broken edges
            nfailedge=2
            goto 10

        else

            write(msg_file,*)'unknown no. of nfailedge!'
            call exit_function

        end if



!-----------------------------------------------------------------------!
!                   UPDATE INTERFACE
!               update global libraries
!-----------------------------------------------------------------------!

!       update element curr_status and sub-element cnc matrices

        if(elstat>elem%curr_status) then
            elem%curr_status=elstat
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

    ! passed in variables
    type(xquad_element) :: elem

    ! extracted variables, from glb libraries
    integer :: edgstat(nedge)               ! status variable array of element edges
    type(real_alloc_array) :: coord(nnode)  ! nodal coord arrays to store the coords of elem nodes extracted from glb node lib

    integer     :: edg(4,nedge)             ! local copy of parameter topo (legacy format, no time to change...)
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


    real(dp) :: xp1, yp1, xp2, yp2      ! (x,y) of point 1 and point 2 on the crack line
    real(dp) :: x1, y1, x2, y2          ! (x,y) of node 1 and node 2 of an element edge
    real(dp) :: xct, yct                ! (x,y) of a crack tip on an edge, i.e., intersection of crack line & edge
    real(dp) :: xo, yo                  ! (x,y) of element centroid


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

        if(elem%curr_status==elfail) return ! elem already failed, no need to proceed

        ! initialize local variables

        edgstat=0; edg=0; xelm=zero; theta=zero

        nfailedge=0; ifedg=0; elstat=0
        jbe1=0; jbe2=0; jnode=0
        iscross=0
        i=0; j=0; l=0; k=0
        subelstat=0

        failed=.false.

        xp1=zero; yp1=zero
        xp2=zero; yp2=zero

        x1=zero; y1=zero
        x2=zero; y2=zero

        xct=zero; yct=zero

        xo=zero; yo=zero




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

        edg=topo

        ! extract material orientation (fibre angle)
        call extract(lib_mat(elem%bulkmat),theta=theta)
        
        elstat=elem%curr_status


!-----------------------------------------------------------------------!
!           check failure criterion on all sub elements
!   elem is judged failed a.l.a. 1 sub elem has reached failure onset
!-----------------------------------------------------------------------!
        do i=1, size(elem%subelem)
            call extract(elem%subelem(i),curr_status=subelstat)
            if(subelstat>intact) failed=.true.
            if(failed) exit
        end do


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
                do i=1,nedge
                    ! tip corods of edge i
                    x1=xelm(1,edg(1,i))
                    y1=xelm(2,edg(1,i))
                    x2=xelm(1,edg(2,i))
                    y2=xelm(2,edg(2,i))
                    iscross=0
                    xct=zero
                    yct=zero
                    call klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct)
                    if (iscross.gt.0) then
                        edgstat(i)=cohcrack  ! edge i cracked
                        nfailedge=nfailedge+1
                        ifedg(nfailedge)=i   ! store failed edge indices
                        xelm(1,edg(3,i))=xct ! store c tip coords on edge i fl. nd 1
                        xelm(2,edg(3,i))=yct
                        xelm(1,edg(4,i))=xct ! store c tip coords on edge i fl. nd 2
                        xelm(2,edg(4,i))=yct
                     endif
                     if(nfailedge.eq.2) exit ! found 2 broken edges already
                end do

            else ! element already partitioned
                !---------- find no. of broken edges -------------------------
                do i=1,nedge
                    if(edgstat(i)/=intact) then
                        nfailedge=nfailedge+1   ! update total no. of damaged edges
                        ifedg(nfailedge)=i  ! update the indices of damaged edges
                    end if
                end do
                !-------------------------------------------------------------
                ! find (xp1,yp1), (xp2,yp2)
                if (elstat .eq. eltrans) then ! transition elm, already has an edge partitioned
                    if(nfailedge .ne. 1) then
                         write(msg_file,*)'inconsistency in kplyfail,elstat=1!'
                         call exit_function
                    end if
                    j=ifedg(1)

                    ! update edgstat(j) & fnode(jfnd)
                    edgstat(j)=cohcrack ! cohesive crack

                    ! find xp1, yp1
                    xp1=xelm(1,edg(3,j))
                    yp1=xelm(2,edg(3,j))

                    ! find xp2, yp2
                    xp2=xp1+cos(theta/halfcirc*pi)
                    yp2=yp1+sin(theta/halfcirc*pi)

                     ! next, find the edge crossed by the crack line
                    do i=1,nedge
                        if (i.eq.j) cycle ! the already cracked edge,go to next edge
                        ! tip corods of edge i
                        x1=xelm(1,edg(1,i))
                        y1=xelm(2,edg(1,i))
                        x2=xelm(1,edg(2,i))
                        y2=xelm(2,edg(2,i))
                        iscross=0
                        xct=zero
                        yct=zero
                        call klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct)
                        if (iscross .gt. 0) then
                            edgstat(i)=cohcrack  ! edge i cracked
                            nfailedge=nfailedge+1
                            ifedg(nfailedge)=i   ! 2nd broken edge index is i
                            xelm(1,edg(3,i))=xct ! store c tip coords on edge i fl. nd 1
                            xelm(2,edg(3,i))=yct
                            xelm(1,edg(4,i))=xct ! store c tip coords on edge i fl. nd 2
                            xelm(2,edg(4,i))=yct
                        end if
                        if(nfailedge .eq. 2) exit ! found 2 broken edges already
                    end do
                else if (elstat.gt.eltrans .and. elstat.lt.elfail) then
                ! element already has two edges partitioned; just update edge status var to coh crack status
                    if(nfailedge.lt.2) then
                         write(msg_file,*)'inconsistency in kplyfail,elstat>1!'
                         call exit_function
                    end if
                    nfailedge=2 ! ignore more than 2 broken edges at the moment
                    ! update first edge status var & fl. node status var
                    j=ifedg(1)
                    edgstat(j)=cohcrack ! cohesive crack
                    ! update second edge status var & fl. node status var
                    j=ifedg(2)
                    edgstat(j)=cohcrack ! cohesive crack
                else
                    write(msg_file,*) 'unsupported elstat value for failure!'
                    call exit_function
                end if

            end if



            ! update the elstat to failed status value
            elstat=elfail

    !       update element curr_status and sub-element cnc matrices

            elem%curr_status=elstat
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

                l=edg(3,j)          ! elem lcl index of broken edge fl. node 1
                coord(l)%array(:)=xelm(:,l)
                k=elem%nodecnc(l)   ! global index of broken edge fl. node 1
                call update(lib_node(k),x=coord(l)%array)

                l=edg(4,j)          ! elem lcl index of broken edge fl. node 2
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
    type(xquad_element),    intent(inout)   :: elem
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


        case (1) !- one edge cracked, trans partition
            ! find the index of the broken edge
            ibe=ifedg(1)

            ! verify its status variable value
            if(edgstat(ibe)/=egtrans) then
                write(msg_file,*)'transition partition only accepts edgstat=1!'
                call exit_function
            end if

            ! allocate sub element arrays; in this case, 3 tri elements
            nsub=3
            if(allocated(elem%subelem)) deallocate(elem%subelem)
            if(allocated(elem%subcnc)) deallocate(elem%subcnc)
            if(allocated(subglbcnc)) deallocate(subglbcnc)
            allocate(elem%subelem(nsub))
            allocate(elem%subcnc(nsub))
            allocate(subglbcnc(nsub))
            do j=1, nsub
                allocate(elem%subcnc(j)%array(3))
                allocate(subglbcnc(j)%array(3))
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

            subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))

            ! sub elm 2 connec
            elem%subcnc(2)%array(1)=topo(1,e2)
            elem%subcnc(2)%array(2)=topo(2,e2)
            elem%subcnc(2)%array(3)=jnode

            subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))

            ! sub elm 3 connec
            elem%subcnc(3)%array(1)=topo(1,e3)
            elem%subcnc(3)%array(2)=topo(2,e3)
            elem%subcnc(3)%array(3)=jnode

            subglbcnc(3)%array(:)=elem%nodecnc(elem%subcnc(3)%array(:))

            ! create sub elements
            call prepare(elem%subelem(1),eltype='tri', matkey=elem%bulkmat, glbcnc=subglbcnc(1)%array)
            call prepare(elem%subelem(2),eltype='tri', matkey=elem%bulkmat, glbcnc=subglbcnc(2)%array)
            call prepare(elem%subelem(3),eltype='tri', matkey=elem%bulkmat, glbcnc=subglbcnc(3)%array)



        case (2) !- two edges cracked

            ibe1=min(ifedg(1),ifedg(2))   ! local edge index of 1st broken edge
            ibe2=max(ifedg(1),ifedg(2))   ! local edge index of 2nd broken edge
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
                    end select
                case(2)
                    select case(ibe2)
                        case(3)
                            nbulk=4
                            e1=2; e2=3; e3=4; e4=1
                        case(4)
                            nbulk=2
                            e1=2; e2=3; e3=4; e4=1
                    end select
                case(3)
                    if(ibe2==4) then
                        nbulk=4
                        e1=3; e2=4; e3=1; e4=2
                    end if
                case default
                        write(msg_file,*)'wrong broken edge in update subcnc xquad'
                        call exit_function
            end select


            select case(nbulk)
                case(2)
                ! two quad bulk subdomains
                    if(iscoh) then
                        nsub=3  ! two quad and one coh2d
                    else
                        nsub=2  ! only two quad
                    end if
                    if(allocated(elem%subelem)) deallocate(elem%subelem)
                    if(allocated(elem%subcnc)) deallocate(elem%subcnc)
                    if(allocated(subglbcnc)) deallocate(subglbcnc)
                    allocate(elem%subelem(nsub))
                    allocate(elem%subcnc(nsub))
                    allocate(subglbcnc(nsub))
                    do j=1, nsub
                        allocate(elem%subcnc(j)%array(4))
                        allocate(subglbcnc(j)%array(4))
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

                    subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))

                    ! sub elm 2 connec
                    elem%subcnc(2)%array(1)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(2)%array(1)=jnode1
                    elem%subcnc(2)%array(2)=topo(2,e1)
                    elem%subcnc(2)%array(3)=topo(1,e3)
                    elem%subcnc(2)%array(4)=topo(3,e3); if(edgstat(e3)<cohcrack) elem%subcnc(2)%array(4)=jnode3

                    subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))

                    ! create sub bulk elements
                    call prepare(elem%subelem(1),eltype='quad', matkey=elem%bulkmat, glbcnc=subglbcnc(1)%array)
                    call prepare(elem%subelem(2),eltype='quad', matkey=elem%bulkmat, glbcnc=subglbcnc(2)%array)

                    if(iscoh) then

                    ! sub elm 3 connec
                    elem%subcnc(3)%array(1)=topo(4,e3); if(edgstat(e3)<cohcrack) elem%subcnc(3)%array(1)=jnode3
                    elem%subcnc(3)%array(2)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(2)=jnode1
                    elem%subcnc(3)%array(3)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(3)=jnode1
                    elem%subcnc(3)%array(4)=topo(3,e3); if(edgstat(e3)<cohcrack) elem%subcnc(3)%array(4)=jnode3

                    subglbcnc(3)%array(:)=elem%nodecnc(elem%subcnc(3)%array(:))

                    call prepare(elem%subelem(3),eltype='coh2d', matkey=elem%cohmat, glbcnc=subglbcnc(3)%array)

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
                        allocate(elem%subcnc(j)%array(3))
                        allocate(subglbcnc(j)%array(3))
                        elem%subcnc(j)%array=0
                        subglbcnc(j)%array=0
                    end do
                    if(iscoh) then
                        allocate(elem%subcnc(5)%array(4))
                        allocate(subglbcnc(5)%array(4))
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

                    subglbcnc(1)%array(:)=elem%nodecnc(elem%subcnc(1)%array(:))

                    ! sub elm 2 connec
                    elem%subcnc(2)%array(1)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(2)%array(1)=jnode2
                    elem%subcnc(2)%array(2)=topo(1,e3)
                    elem%subcnc(2)%array(3)=topo(2,e3)

                    subglbcnc(2)%array(:)=elem%nodecnc(elem%subcnc(2)%array(:))

                    ! sub elm 3 connec
                    elem%subcnc(3)%array(1)=topo(1,e4)
                    elem%subcnc(3)%array(2)=topo(2,e4)
                    elem%subcnc(3)%array(3)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(3)%array(3)=jnode1

                    subglbcnc(3)%array(:)=elem%nodecnc(elem%subcnc(3)%array(:))

                    ! sub elm 4 connec
                    elem%subcnc(4)%array(1)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(4)%array(1)=jnode1
                    elem%subcnc(4)%array(2)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(4)%array(2)=jnode2
                    elem%subcnc(4)%array(3)=topo(2,e3)

                    subglbcnc(4)%array(:)=elem%nodecnc(elem%subcnc(4)%array(:))

                    ! create sub bulk elements
                    call prepare(elem%subelem(1),eltype='tri', matkey=elem%bulkmat, glbcnc=subglbcnc(1)%array)
                    call prepare(elem%subelem(2),eltype='tri', matkey=elem%bulkmat, glbcnc=subglbcnc(2)%array)
                    call prepare(elem%subelem(3),eltype='tri', matkey=elem%bulkmat, glbcnc=subglbcnc(3)%array)
                    call prepare(elem%subelem(4),eltype='tri', matkey=elem%bulkmat, glbcnc=subglbcnc(4)%array)

                    if(iscoh) then

                    ! sub elm 5 connec
                    elem%subcnc(5)%array(1)=topo(4,e1); if(edgstat(e1)<cohcrack) elem%subcnc(5)%array(1)=jnode1
                    elem%subcnc(5)%array(2)=topo(3,e2); if(edgstat(e2)<cohcrack) elem%subcnc(5)%array(2)=jnode2
                    elem%subcnc(5)%array(3)=topo(4,e2); if(edgstat(e2)<cohcrack) elem%subcnc(5)%array(3)=jnode2
                    elem%subcnc(5)%array(4)=topo(3,e1); if(edgstat(e1)<cohcrack) elem%subcnc(5)%array(4)=jnode1

                    subglbcnc(5)%array(:)=elem%nodecnc(elem%subcnc(5)%array(:))

                    call prepare(elem%subelem(5),eltype='coh2d', matkey=elem%cohmat, glbcnc=subglbcnc(5)%array)

                    end if

                case default
                    write(msg_file,*)'wrong nbulk in update subcnc xquad'
                    call exit_function
            end select


        !~case(3) !- three edges crack
        !~   write(msg_file,*) 'three-edge crack partition not yet supported!'
        !~   nfailedge=2
        !~   goto 10
        !~
        !~case(4) !- four edges crack
        !~   write(msg_file,*) 'four-edge crack partition not yet supported!'
        !~   nfailedge=2
        !~   goto 10

        case default
           write(msg_file,*) 'WARNING: xquad update subcnc case selection default!'

        end select


        ! deallocate local array

        if(allocated(subglbcnc)) deallocate(subglbcnc)


    end subroutine update_subcnc












end module xquad_element_module
