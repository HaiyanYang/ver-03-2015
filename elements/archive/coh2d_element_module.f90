    module coh2d_element_module
    use parameter_module
    use integration_point_module        ! integration point type
    
    implicit none
    private
    
    integer, parameter          :: ndim=2, nst=2, nnode=4, nig=2, ndof=ndim*nnode ! constants for type coh2d_element 
    real(kind=dp), parameter    :: dfail=one                                      ! max. degradation at final failure
    
    
    type, public :: coh2d_element 
        private
        
        integer :: curr_status=0
        integer :: key=0                            ! glb index of this element
        integer :: connec(nnode)=0                  ! node indices in their global arrays
        integer :: matkey=0                         ! material index in the global material arrays
        type(integration_point) :: ig_point(nig)    ! x, xi, weight, stress, strain, sdv; initialize in prepare procedure
        
        ! below are optional terms 
        
        type(sdv_array), allocatable :: sdv(:)
        
    end type
    
    
    
    
    interface empty
        module procedure empty_coh2d_element
    end interface
    
    interface prepare
        module procedure prepare_coh2d_element
    end interface
    
    interface integrate
        module procedure integrate_coh2d_element
    end interface
    
    interface extract
        module procedure extract_coh2d_element
    end interface
    
    


    public :: empty,prepare,integrate,extract
    
    
    
    
    contains
    
    
    ! this subroutine is used to format the element for use
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine empty_coh2d_element(elem)
    
        type(coh2d_element),intent(out) ::elem
        
        integer :: i
        i=0
        
        elem%curr_status=0
        elem%key=0
        elem%connec=0
        elem%matkey=0

        do i=1,nig
            call empty(elem%ig_point(i))
        end do
          
        if(allocated(elem%sdv)) deallocate(elem%sdv)
    
    end subroutine empty_coh2d_element
    
    
    
    
    ! this subroutine is used to prepare the connectivity and material lib index of the element
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine prepare_coh2d_element(elem,key,connec,matkey)
    
        type(coh2d_element),    intent(inout)   :: elem
        integer,                intent(in)      :: connec(nnode)
        integer,                intent(in)      :: key,matkey
        
        real(kind=dp)   :: x(ndim),u(ndim),stress(nst),strain(nst)
        integer         :: i
        x=zero; u=zero; stress=zero; strain=zero
        i=0
        
        elem%key=key
        elem%connec=connec
        elem%matkey=matkey
        
        do i=1,nig
            call update(elem%ig_point(i),x=x,u=u,stress=stress,strain=strain)
        end do
    
    end subroutine prepare_coh2d_element
    
    
    
    
    subroutine extract_coh2d_element(elem,curr_status,key,connec,matkey,ig_point,sdv)
    
        type(coh2d_element), intent(in) :: elem
        
        integer,                              optional, intent(out) :: curr_status,key, matkey
        integer,                 allocatable, optional, intent(out) :: connec(:)
        type(integration_point), allocatable, optional, intent(out) :: ig_point(:)
        type(sdv_array),         allocatable, optional, intent(out) :: sdv(:)
        
        if(present(curr_status)) curr_status=elem%curr_status
        
        if(present(key)) key=elem%key
        
        if(present(matkey)) matkey=elem%matkey
        
        if(present(connec)) then
            allocate(connec(nnode))
            connec=elem%connec
        end if
        
        if(present(ig_point)) then
            allocate(ig_point(nig))
            ig_point=elem%ig_point
        end if
        
        if(present(sdv)) then        
            if(allocated(elem%sdv)) then
                allocate(sdv(size(elem%sdv)))
                sdv=elem%sdv
            end if
        end if
         
    end subroutine extract_coh2d_element


    
    
    
    
    ! the integration subroutine, updates K matrix, F vector, integration point stress and strain
    ! as well as all the solution dependent variables (sdvs) at intg points and element
    subroutine integrate_coh2d_element(elem,K_matrix,F_vector,gauss,mnode)
        use toolkit_module                  ! global tools for element integration
        use lib_mat_module                  ! global material library
        use lib_node_module                 ! global node library
        use glb_clock_module                ! global analysis progress (curr. step, inc, time, dtime)
    
        type(coh2d_element),        intent(inout)   :: elem 
        real(kind=dp), allocatable, intent(out)     :: K_matrix(:,:), F_vector(:)
        logical,        optional,   intent(in)      :: gauss     ! true for gauss integration
        type(xnode),    optional,   intent(in)      :: mnode(:)  ! material (interpolated) coords and u
 
        !-----------------------------------------!
        ! local variables, extracted from glb libs
        !-----------------------------------------!
        
        ! - element nodes extracted from global node library
        type(xnode)     :: node(nnode)                  ! x, u, du, v, extra dof ddof etc
        
        ! - element material extracted from global material library
        type(material)  :: mat                          ! matname, mattype and matkey to glb mattype array   

        ! - glb clock step and increment no. extracted from glb clock module
        integer         :: curr_step, curr_inc
        
        
        !-----------------------------------------!
        ! - the rest are all pure local variables
        !-----------------------------------------!
        
        ! - variables extracted from element nodes
        real(kind=dp),allocatable :: xj(:),uj(:)        ! nodal vars extracted from glb lib_node array
        real(kind=dp)   :: coords(ndim,nnode)           ! coordinates of the element nodes, formed from xj of each node
        real(kind=dp)   :: u(ndof)                      ! element nodal disp. vector, formed from uj of each node
        
        ! - variables extracted from element material
        character(len=matnamelength)    :: matname      ! name of the material assigned to this element
        character(len=mattypelength)    :: mattype      ! type of the material
        integer                         :: typekey       ! index of the material in the material library of its type 
        
        ! - variables extracted from element isdv
        integer         :: nstep, ninc                  ! step and increment no. of the last iteration, stored in the element
        logical         :: last_converged               ! true if last iteration has converged: a new increment/step has started
        
        ! - variables extracted from intg point sdvs
        type(sdv_array),  allocatable   :: ig_sdv(:)
        
        
        ! - variables defined locally
        
        real(kind=dp)   :: igxi(ndim-1,nig),igwt(nig)   ! ig point natural coords and weights
        real(kind=dp)   :: tmpx(ndim), tmpu(ndim)       ! temporary arrays to hold x and u values for intg points
               
        real(kind=dp)   :: fn(nnode)                    ! shape functions
        real(kind=dp)   :: Nmatrix(ndim,ndof)           ! obtained from fn, to compute disp. jump acrss intfc: {u}_jump = [N]*{u}
        real(kind=dp)   :: normal(ndim),tangent(ndim)   ! normal and tangent vectors of the interface, obtained from coords
        real(kind=dp)   :: det                          ! determinant of jacobian (=length of element/2); 2 is length of ref. elem
        real(kind=dp)   :: Qmatrix(ndim,ndim)           ! rotation matrix from global to local coordinates (from normal & tangent)
        real(kind=dp)   :: ujump(ndim),delta(ndim)      ! {delta}=[Q]*{u}_jump, {delta} is the jump vector in lcl coords.
        real(kind=dp)   :: Dee(ndim,ndim)               ! material stiffness matrix [D]
        real(kind=dp)   :: Tau(ndim)                    ! {Tau}=[D]*{delta}, traction on the interface
        
        real(kind=dp)   :: QN(ndim,ndof),DQN(ndim,ndof)         ! [Q]*[N], [D]*[Q]*[N]
        real(kind=dp)   :: NtQt(ndof,ndim),NtQtDQN(ndof,ndof)   ! [N']*[Q'], [N']*[Q']*[D]*[Q]*[N] 
        real(kind=dp)   :: NtQtTau(ndof)                        ! [N']*[Q']*{Tau}
        integer         :: i,j,k,kig,igstat
      
        
        
        
        !------------------------------------------------!
        !           initialize all variables 
        !------------------------------------------------!
        
        ! intent(out) variables, automatically deallocated when passed in
        allocate(K_matrix(ndof,ndof),F_vector(ndof)); K_matrix=zero; F_vector=zero 
        
        ! integer counters
        i=0; j=0; k=0; kig=0
        
        ! local variables, extracted from glb libs
        do i=1,nnode
            call empty(node(i))
        end do 
        call empty(mat)
        curr_step=0; curr_inc=0
        
        ! pure local variables
        coords=zero; u=zero 
        matname=''; mattype=''; typekey=0 
        nstep=0; ninc=0; last_converged=.false.
        igxi=zero; igwt=zero
        tmpx=zero; tmpu=zero
        fn=zero; Nmatrix=zero
        tangent=zero; normal=zero 
        det=zero; Qmatrix=zero
        ujump=zero; delta=zero
        Dee=zero; Tau=zero
        QN=zero; DQN=zero; NtQt=zero; NtQtDQN=zero; NtQtTau=zero
        
        
        
        
        
        !------------------------------------------------!
        !   extract variables from global arrays 
        !------------------------------------------------!
        
        if(present(mnode)) then
            ! - extract nodes from passed-in node array
            node(:)=mnode(:)
        else
            ! - extract nodes from global node array 
            node(:)=lib_node(elem%connec(:))
        end if
        
        ! - extract material values from global material array
        mat=lib_mat(elem%matkey)
        
        ! - extract curr step and inc values from glb clock module
        call extract_glb_clock(kstep=curr_step,kinc=curr_inc)
        
        
        
        
        
        !------------------------------------------------!
        !   assign values to material, coords and u 
        !------------------------------------------------!      
        
        ! - extract x and u values from nodes and assign to local arrays 
        do j=1,nnode
            ! extract x and u values from nodes
            call extract(node(j),x=xj,u=uj)     
            ! assign x to coords matrix
            if(allocated(xj)) then
                coords(:,j)=xj(:)
            else
                write(msg_file,*)'WARNING: x not allocated for node:',elem%connec(j)
            end if
            ! and u to u vector
            if(allocated(uj)) then 
                u((j-1)*ndim+1:j*ndim)=uj(1:ndim)
            else
                write(msg_file,*)'WARNING: u not allocated for node:',elem%connec(j)
            end if    
        end do
        
        
        ! - extract values from mat (material type) and assign to local vars (matname, mattype & matkey)
        call extract(mat,matname,mattype,typekey) 
        
        
        ! - extract isdv values from element and assign to nstep and ninc
        if(.not.allocated(elem%sdv)) then   ! 1st iteration
            allocate(elem%sdv(1))
            allocate(elem%sdv(1)%i(2))      ! allocate integer sdv array
            elem%sdv(1)%i(1) = curr_step    ! store current step & increment in the integer sdv array
            elem%sdv(1)%i(2) = curr_inc
        end if
        nstep        = elem%sdv(1)%i(1)     ! extract the step & increment no. of the last iteration
        ninc         = elem%sdv(1)%i(2)
        
        
        ! check if last iteration has converged, and if so, update logical var. 
        ! and store new step & iteration values
        if(nstep.ne.curr_step .or. ninc.ne.curr_inc) then
            last_converged=.true.
            elem%sdv(1)%i(1) = curr_step    ! update the current step & increment no.
            elem%sdv(1)%i(2) = curr_inc
        end if
        
        
        
        
        !------------------------------------------------!
        !   compute Q matrix (rotation) and determinat 
        !------------------------------------------------!
        
        ! - compute tangent of the interface: node 2,3 coords - node 1,4 coords
        tangent(1)=half*(coords(1,2)+coords(1,3)) - half*(coords(1,1)+coords(1,4))
        tangent(2)=half*(coords(2,2)+coords(2,3)) - half*(coords(2,1)+coords(2,4))
        
        ! - normalize tangent vector
        call normalize(tangent,det) !-tangent vector normalized
        
        ! - compute determinant of the line element: actual length/reference length
        det=det/two !-ref. line starts from -1 to 1, length=2

        ! - compute normal of the interface
        normal(1)=-tangent(2)
        normal(2)=tangent(1)

        ! - compute Q matrix: rotate global coords to local coords; Q = transpose[normal,tangent]
        do j=1,2
            Qmatrix(1,j)=normal(j)
            Qmatrix(2,j)=tangent(j)
        end do
        
        
        
        
        
        
        !------------------------------------------------!
        !   perform integration at ig points
        !------------------------------------------------!
        
        ! - calculate ig point xi and weight
        if(present(gauss).and.gauss)  then
            call init_ig(igxi,igwt,gauss)
        else                          
            call init_ig(igxi,igwt)
        end if
         
        !-calculate strain,stress,stiffness,sdv etc. at each int point
      	do kig=1,nig 
        
            ! - empty relevant arrays for reuse
            fn=zero; Nmatrix=zero
            ujump=zero; delta=zero; dee=zero; Tau=zero
            QN=zero; NtQt=zero; DQN=zero; NtQtDQN=zero; NtQtTau=zero
            tmpx=zero; tmpu=zero
            
            !- get shape function matrix
            call init_shape(igxi(:,kig),fn) 
            
		    ! Nmatrix: ujump (at each int pnt) = Nmatrix*u
            do i = 1,ndim
                do j = 1,nnode/2
                Nmatrix(i,i+(j-1)*ndim)= -fn(j)
                Nmatrix(i,i+(nnode-j)*ndim)=fn(j)
                end do
            end do
            
	   	    ! calculate ujump: disp. jump of the two crack surface, in global coords
            ujump=matmul(Nmatrix,u)
            
            ! calculate separation delta in local coords: delta=Qmatrix*ujump
            delta=matmul(Qmatrix,ujump)
            
            ! - extract sdvs from integration points; ig_sdv automatically deallocated when passed in
            call extract(elem%ig_point(kig),sdv=ig_sdv)
            
            ! allocate ig_sdv arrays for 1st iteration of analysis
            if(.not.allocated(ig_sdv)) then
            ! allocate 2 sets of sdv arrays, 1 for converged sdvs and 1 for iterating sdvs
                allocate(ig_sdv(2))
            end if
            
            ! update converged sdvs (sdv1) with iterating sdvs (sdv2) when last iteration has converged
            ! and revalue iterating sdvs (sdv2) to the last converged sdvs (sdv1) if otherwise
            if(last_converged) then
                ig_sdv(1)=ig_sdv(2)
            else               
                ig_sdv(2)=ig_sdv(1)
            end if
            
            ! get D matrix dee accord. to material properties, and update intg point variables
            select case (mattype)
                case ('interface')
                    
                    ! calculate D matrix, update stress, and iterating sdv
                    call ddsdde(lib_interface(typekey), Dee, jump=delta, stress=Tau, sdv=ig_sdv(2), dfail=dfail) 
                    
                case default
                    write(msg_file,*) 'material type not supported for cohesive element!'
                    call exit_function
            end select
            
            

            !------------------------------------------------!
            !  add this ig point contributions to K and F
            !------------------------------------------------!
          
            QN      =   matmul(Qmatrix, Nmatrix)
            NtQt    =   transpose(QN)
            DQN     =   matmul(Dee, QN)
            NtQtDQN =   matmul(NtQt, DQN)
            NtQtTau =   matmul(NtQt, Tau)
		
            do j=1, ndof
                do k=1, ndof
                    K_matrix(j,k) = K_matrix(j,k) + NtQtDQN(j,k) * det * igwt(kig)
                end do
                F_vector(j) = F_vector(j) + NtQtTau(j) * igwt(kig) * det
            end do
               
            
            
            !------------------------------------------------!
            !   update ig point kig
            !------------------------------------------------!

            !- calculate integration point physical coordinates (initial)
            tmpx    =   matmul(coords,fn)
            
            !- calculate integration point displacement
            do j=1, ndim
                do i=1, nnode
                    tmpu(j) = tmpu(j) + fn(i) * u((i-1)*ndim+j)
                end do
            end do
            
            ! update element ig point arrays
            call update(elem%ig_point(kig),x=tmpx,u=tmpu,strain=delta,stress=Tau,sdv=ig_sdv)
            
            ! update elem curr status variable
            igstat=0
            if(allocated(ig_sdv(2)%i)) igstat=ig_sdv(2)%i(1)
            elem%curr_status=max(elem%curr_status,igstat)
            
            
       	end do !-looped over all int points. ig=nig

        if(allocated(xj)) deallocate(xj) 
        if(allocated(uj)) deallocate(uj) 
        if(allocated(ig_sdv)) deallocate(ig_sdv) 
          
    
    end subroutine integrate_coh2d_element
    
    
    
    
    
    
    
    
    
    
    
    
    ! the rest are private subroutines
    
    
    
    
    
    

    
 
    subroutine init_ig(xi,wt,gauss)

      real(kind=dp), intent(inout)  :: xi(ndim-1,nig), wt(nig)
      logical, optional, intent(in) :: gauss
      
      real(kind=dp) :: cn=zero
	
        cn=one !-newton cotes
    
        if(present(gauss).and.gauss) cn=0.5773502691896260_dp !-gauss ig point
    
        if (nig .eq. 2) then          
            xi(1,1)=-cn
            xi(1,2)=cn
            wt=one
        else
            write(msg_file,*) 'no. of integration points incorrect for coh2d_ig!'
            call exit_function
        end if

    end subroutine init_ig
    
    
    
    subroutine init_shape(igxi,f)
      
        real(kind=dp),intent(inout) :: f(nnode)
        real(kind=dp),intent(in) :: igxi(ndim-1)
        
        real(kind=dp) :: xi ! local variables
        xi=zero
        
        xi=igxi(1)
        
        f(1)=half*(one-xi)
        f(2)=half*(one+xi)
        f(3)=f(2)
        f(4)=f(1)

    end subroutine init_shape
    
    
    end module coh2d_element_module
