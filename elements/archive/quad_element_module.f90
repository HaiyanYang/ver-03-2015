    module quad_element_module
    use parameter_module
    use integration_point_module        ! integration point type
    
    implicit none
    private
    
    integer, parameter :: ndim=2, nst=3, nnode=4, nig=4, ndof=ndim*nnode ! constants for type quad_element 
    real(dp),parameter :: dfail=one                                      ! max. degradation at final failure 
    
    type, public :: quad_element 
        private
        
        integer :: curr_status=0
        integer :: key=0 ! glb index of this element
        integer :: connec(nnode)=0 ! node indices in their global arrays
        integer :: matkey=0 ! material index in the global material arrays
        real(kind=dp) :: theta=zero ! material (local) orientation for composite lamina
        type(integration_point) :: ig_point(nig) ! x, xi, weight, stress, strain, sdv; initialize in prepare procedure
        
        ! below are optional terms 
        
        type(sdv_array), allocatable :: sdv(:)
        
    end type
    
    
    
    
    interface empty
        module procedure empty_quad_element
    end interface
    
    interface prepare
        module procedure prepare_quad_element
    end interface
    
    interface integrate
        module procedure integrate_quad_element
    end interface
    
    interface extract
        module procedure extract_quad_element
    end interface
    
    


    public :: empty,prepare,integrate,extract
    
    
    
    
    contains
    
    
    ! this subroutine is used to format the element for use
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine empty_quad_element(elem)
    
        type(quad_element),intent(out) ::elem
        
        integer :: i
        i=0
        
        elem%curr_status=0
        elem%key=0
        elem%connec=0
        elem%matkey=0
        elem%theta=zero
        do i=1,nig
            call empty(elem%ig_point(i))
        end do       
   
        if(allocated(elem%sdv)) deallocate(elem%sdv) 
    
    end subroutine empty_quad_element
    
    
    
    
    ! this subroutine is used to prepare the connectivity and material lib index of the element
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine prepare_quad_element(elem,key,connec,matkey,theta)
    
        type(quad_element),     intent(inout)   :: elem
        integer,                intent(in)      :: connec(nnode)
        integer,                intent(in)      :: key,matkey
        real(kind=dp),optional, intent(in)      :: theta
        
        real(kind=dp)   :: x(ndim),u(ndim),stress(nst),strain(nst)
        integer         :: i
        x=zero; u=zero; stress=zero; strain=zero
        i=0
        
        elem%key=key
        elem%connec=connec
        elem%matkey=matkey
        if(present(theta)) elem%theta=theta
        
        do i=1,nig
            call update(elem%ig_point(i),x=x,u=u,stress=stress,strain=strain)
        end do
    
    end subroutine prepare_quad_element
    
    
    
    
    subroutine extract_quad_element(elem,curr_status,key,connec,matkey,theta,ig_point,sdv)
    
        type(quad_element), intent(in) :: elem
        
        integer,                              optional, intent(out) :: curr_status, key, matkey
        real(kind=dp),                        optional, intent(out) :: theta
        integer,                 allocatable, optional, intent(out) :: connec(:)
        type(integration_point), allocatable, optional, intent(out) :: ig_point(:)
        type(sdv_array),         allocatable, optional, intent(out) :: sdv(:)
        
        if(present(curr_status)) curr_status=elem%curr_status
        
        if(present(key)) key=elem%key
        
        if(present(matkey)) matkey=elem%matkey
        
        if(present(theta)) theta=elem%theta
        
        
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
    
    
    end subroutine extract_quad_element


    
    
    
    
    ! the integration subroutine, updates K matrix, F vector, integration point stress and strain
    ! as well as all the solution dependent variables (sdvs) at intg points and element
    subroutine integrate_quad_element(elem,K_matrix,F_vector,planestrain)
    use toolkit_module                  ! global tools for element integration
    use lib_mat_module                  ! global material library
    use lib_node_module                 ! global node library
    use glb_clock_module                ! global analysis progress (curr. step, inc, time, dtime)
    
        type(quad_element),intent(inout)        :: elem 
        real(kind=dp),allocatable,intent(out)   :: K_matrix(:,:), F_vector(:)
        logical, optional, intent(in)           :: planestrain
        
        ! the rest are all local variables
        
        ! variables to be extracted from global arrays
        type(xnode) :: node(nnode) ! x, u, du, v, extra dof ddof etc
        type(material) :: mat ! matname, mattype and matkey to glb mattype array
        character(len=matnamelength) :: matname
        character(len=mattypelength) :: mattype
        integer :: typekey
        ! - glb clock step and increment no. extracted from glb clock module
        integer         :: curr_step, curr_inc
        
        ! - variables extracted from element isdv
        integer         :: nstep, ninc                  ! step and increment no. of the last iteration, stored in the element
        logical         :: last_converged               ! true if last iteration has converged: a new increment/step has started
        
        ! - variables extracted from intg point sdvs
        type(sdv_array),  allocatable   :: ig_sdv(:)
        
        ! variables defined locally
        real(kind=dp)   :: igxi(ndim,nig),igwt(nig) ! ig point natural coords and weights
        real(kind=dp)   :: coords(ndim,nnode) ! coordinates of the element nodes
        real(kind=dp)   :: theta ! orientation of element local coordinates
        real(kind=dp)   :: u(ndof) ! nodal disp. vector
        logical         :: plstrain ! true for plane strain stiffness
        logical         :: failure ! true for failure analysis
        real(kind=dp)   :: fn(nnode),dn(nnode,ndim) ! shape functions & their deriv. physical space
        real(kind=dp)   :: jac(ndim,ndim),gn(nnode,ndim),detj ! jacobian & shape func. deriv. natural space
        real(kind=dp)   :: bee(nst,ndof),beet(ndof,nst) ! b matrix and its transpose
        real(kind=dp)   :: dee(nst,nst) ! d matrix
        real(kind=dp)   :: btd(ndof,nst),btdb(ndof,ndof) ! b'*d & b'*d*b
        real(kind=dp)   :: tmpx(ndim),tmpu(ndim),tmpstrain(nst),tmpstress(nst) ! temp. x, strain & stress arrays for intg pnts      
        real(kind=dp),allocatable :: xj(:),uj(:)! nodal vars extracted from glb lib_node array
        real(kind=dp),allocatable :: vec(:,:)
        real(kind=dp)   :: vecf(2), c, s, clength
        
        integer :: i,j,kig,igstat
        
        ! initialize variables
        allocate(K_matrix(ndof,ndof),F_vector(ndof))
        K_matrix=zero; F_vector=zero
        
        i=0; j=0; kig=0
        do i=1,nnode
            call empty(node(i))
        end do
        call empty(mat)
        curr_step=0; curr_inc=0
        
        igxi=zero; igwt=zero
        coords=zero; theta=zero; u=zero
        plstrain=.false.; failure=.false.
        
        ! copy nodes from global node array 
        node(:)=lib_node(elem%connec(:))
        
        ! assign values to local arrays from nodal values
        do j=1,nnode
        
            ! extract useful values from nodes
            call extract(node(j),x=xj,u=uj)
            
            ! assign values to coords matrix and u vector
            if(allocated(xj)) then
                coords(:,j)=xj(:)
            else
                write(msg_file,*)'WARNING: x not allocated for node:',elem%connec(j)
            end if
            
            if(allocated(uj)) then 
                u((j-1)*ndim+1:j*ndim)=uj(1:ndim)
            else
                write(msg_file,*)'WARNING: u not allocated for node:',elem%connec(j)
            end if
            
        end do
        
        ! extract material values from global material array
        mat=lib_mat(elem%matkey)
        call extract(mat,matname,mattype,typekey,theta)
        
        ! - extract curr step and inc values from glb clock module
        call extract_glb_clock(kstep=curr_step,kinc=curr_inc)
        
        ! extract plain strain variable from element
        if(present(planestrain)) plstrain=planestrain
        
        
        
        ! calculate approximate clength
        vecf=zero; c=zero; s=zero; clength=zero
        allocate(vec(2,2)); vec=zero
        vec(1:2,1)=coords(1:2,3)-coords(1:2,1) ! diagonal vector 1
        vec(1:2,2)=coords(1:2,4)-coords(1:2,2) ! diagonal vector 2 
        c=cos(pi*theta/halfcirc)
        s=sin(pi*theta/halfcirc)
        vecf=[c,s]
        clength=max(dot_product(vecf,vec(:,1)),dot_product(vecf,vec(:,2)))


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
        
        
        ! update ig point xi and weight
        call init_ig(igxi,igwt)
          
        
        
        !-calculate strain,stress,stiffness,sdv etc. at each int point
      	do kig=1,nig 
        
            ! empty relevant arrays for reuse
            fn=zero; dn=zero
            jac=zero; gn=zero; detj=zero
            bee=zero; beet=zero; dee=zero
            btd=zero; btdb=zero
            tmpx=zero; tmpu=zero; tmpstrain=zero; tmpstress=zero
            
            !- get shape matrix and derivatives
            call init_shape(igxi(:,kig),fn,dn) 
            
            !- calculate integration point physical coordinates (initial)
            tmpx=matmul(coords,fn)
            
            !- calculate integration point displacement
            do j=1,ndim
                do i=1,nnode
                    tmpu(j)=tmpu(j)+fn(i)*u((i-1)*ndim+j)
                end do
            end do
            
            ! get jacobian
            jac=matmul(coords,dn)
            
            !-get determinant of jacobian
            detj=determinant(jac)
            
            ! invert jac onto itself
            jac=inverse(jac,detj)
            
            ! calculate gradient of shape function matrix
            gn=matmul(dn,jac)
            
            !-obtain b matrix (nst*ndof) from rearranging terms of gn
            bee=beemat(gn)
            
            ! calculate global strains
            tmpstrain=matmul(bee,u)
            
            ! transfer strain to local coordinates
            if(theta/=zero) tmpstrain=lcl_strain(tmpstrain,theta)
            
            
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
            
            
            ! get D matrix dee accord. to material properties, and update intg point stresses
            select case (mattype)
                case ('isotropic')
                    
                    ! check if failure analysis is needed (check if strength parameters are present)
                    call extract(lib_iso(typekey),strength_active=failure)
                    
                    if(failure) write(msg_file,*) "WARNING: failure analysis is not yet supported for &
                    & quad_element type isotropic material; only linear elastic stiffness matrix is integrated."
                    
                    ! calculate D matrix, update tmpstress
                    call ddsdde(lib_iso(typekey),dee,strain=tmpstrain,stress=tmpstress,PlaneStrain=plstrain) 
                    
                case ('lamina')
                
                    ! check if failure analysis is needed (check if strength parameters are present)
                    call extract(lib_lamina(typekey),strength_active=failure)
                    
                    if(failure) then
                        ! calculate D matrix, update tmpstress and sdv
                        call ddsdde(lib_lamina(typekey),clength,dee,strain=tmpstrain,stress=tmpstress,PlaneStrain=plstrain &
                        & ,sdv=ig_sdv(2),dfail=dfail)
                    else
                        ! calculate D matrix, update tmpstress
                        call ddsdde(lib_lamina(typekey),dee=dee,strain=tmpstrain,stress=tmpstress,PlaneStrain=plstrain)
                    end if
                    
                case default
                    write(msg_file,*) 'material type not supported for tri element!'
                    call exit_function
            end select
            
            ! get D matrix in global coordinates deeg
            if(theta/=zero) dee=glb_dee(dee,theta)
            
            ! calculate B' D B
            beet=transpose(bee)
            btd=matmul(beet,dee)
            btdb=matmul(btd,bee)

            ! integrate and update K matrix
       		do i=1,ndof
          		do j=1,ndof
              		K_matrix(i,j) = K_matrix(i,j)+btdb(i,j)*detj*igwt(kig) !-gauss integration
          		end do
        	end do	
            
            
            ! update ig point arrays
            call update(elem%ig_point(kig),x=tmpx,u=tmpu,strain=tmpstrain,stress=tmpstress,sdv=ig_sdv)
            
            ! update elem curr status variable
            igstat=0
            if(allocated(ig_sdv(2)%i)) igstat=ig_sdv(2)%i(1)
            elem%curr_status=max(elem%curr_status,igstat)
            
       	end do !-looped over all int points. ig=nig
        
        F_vector=matmul(K_matrix,u) 
        
        
        ! deallocate local dynamic arrays
        if(allocated(xj)) deallocate(xj) 
        if(allocated(uj)) deallocate(uj) 
        if(allocated(vec)) deallocate(vec) 
        if(allocated(ig_sdv)) deallocate(ig_sdv)
        
        
    
    end subroutine integrate_quad_element
    
    
    
    
    
    
    
    
    
    
    
    
    ! the rest are private subroutines
    
    
    
    
    
    

    
 
    subroutine init_ig(xi,wt)

      real(kind=dp),intent(inout) :: xi(ndim,nig),wt(nig)
	
        if (nig .eq. 1) then
            xi(1,1)= zero
            xi(2,1)= zero
            wt = four
        else if(nig .eq. 4) then
            xi(1,1)= -root3
            xi(2,1)= -root3
            xi(1,2)= root3
            xi(2,2)= -root3
            xi(1,3)= -root3
            xi(2,3)= root3
            xi(1,4)= root3
            xi(2,4)= root3
            wt = one
        else
            write(msg_file,*) 'no. of integration points incorrect for quad_ig!'
            call exit_function
        end if

    end subroutine init_ig
    
    
    
    subroutine init_shape(igxi,f,df)
      
        real(kind=dp),intent(inout) :: f(nnode),df(nnode,ndim)
        real(kind=dp),intent(in) :: igxi(ndim)
        
        real(kind=dp) :: xi,eta ! local variables
        xi=zero
        eta=zero
        
        xi=igxi(1)
        eta=igxi(2)
        
        f(1)=quarter*(one-xi)*(one-eta)
        f(2)=quarter*(one+xi)*(one-eta)
        f(3)=quarter*(one+xi)*(one+eta)
        f(4)=quarter*(one-xi)*(one+eta)
        df(1,1) = -quarter*(one-eta)
        df(2,1) =  quarter*(one-eta)
        df(3,1) =  quarter*(one+eta)
        df(4,1) = -quarter*(one+eta)
        df(1,2) = -quarter*(one-xi)
        df(2,2) = -quarter*(one+xi)
        df(3,2) =  quarter*(one+xi)
        df(4,2) =  quarter*(one-xi)

    end subroutine init_shape
    
    
    end module quad_element_module
