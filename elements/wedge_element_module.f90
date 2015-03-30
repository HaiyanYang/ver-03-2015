    module wedge_element_module
    use parameter_module
    use integration_point_module        ! integration point type
    
    implicit none
    private
    
    integer, parameter :: ndim=3, nst=6, nnode=6, nig=6, ndof=ndim*nnode ! constants for type wedge_element 
    real(dp),parameter :: dfail=one                                      ! max. degradation at final failure 
   
    ! parameters used for calculating element characteristic length (clength)
    integer, parameter :: nedge=3
    integer, parameter :: edg(2,nedge)=reshape([1,2,2,3,3,1],[2,nedge])

 
    type, public :: wedge_element 
        private
        
        integer :: curr_status=0
        integer :: key=0 ! glb index of this element
        integer :: connec(nnode)=0 ! node indices in their global arrays
        integer :: matkey=0 ! material index in the global material arrays
        
        real(dp):: plyangle=zero        ! ply angle for composite lamina (rotation around z axis) 

        real(dp) :: stress(nst)=zero ! stresses for output
        real(dp) :: strain(nst)=zero ! strains for output

        type(integration_point) :: ig_point(nig) ! x, xi, weight, stress, strain, sdv; initialize in prepare procedure
        
        integer :: nstep=0, ninc=0      ! to store curr step and increment no.
        
        ! below are optional terms 
        
        type(sdv_array), allocatable :: sdv(:)
        
    end type
    
    
    
    
    interface empty
        module procedure empty_wedge_element
    end interface
    
    interface prepare
        module procedure prepare_wedge_element
    end interface
    
    interface integrate
        module procedure integrate_wedge_element
    end interface
    
    interface extract
        module procedure extract_wedge_element
    end interface
    
    


    public :: empty,prepare,integrate,extract
    
    
    
    
    contains
    
    
    ! this subroutine is used to format the element for use
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine empty_wedge_element(elem)
    
        type(wedge_element),intent(out) ::elem
        
        integer :: i
        i=0
        
        elem%curr_status=0
        elem%key=0
        elem%connec=0
        elem%matkey=0
        elem%plyangle=zero
        elem%stress=zero
        elem%strain=zero
        
        elem%nstep=0
        elem%ninc=0

        do i=1,nig
            call empty(elem%ig_point(i))
        end do
      
        if(allocated(elem%sdv)) deallocate(elem%sdv)
    
    end subroutine empty_wedge_element
    
    
    
    
    ! this subroutine is used to prepare the connectivity and material lib index of the element
    ! it is used in the initialize_lib_elem procedure in the lib_elem module
    subroutine prepare_wedge_element(elem,key,connec,matkey,plyangle)
    
        type(wedge_element),    intent(inout)   :: elem
        integer,                intent(in)      :: connec(nnode)
        integer,                intent(in)      :: key,matkey
        real(dp),               intent(in)      :: plyangle

        
        real(kind=dp)   :: x(ndim),u(ndim),stress(nst),strain(nst)
        integer         :: i
        x=zero; u=zero; stress=zero; strain=zero
        i=0
        
        elem%key=key
        elem%connec=connec
        elem%matkey=matkey
        elem%plyangle=plyangle
        
        do i=1,nig
            call update(elem%ig_point(i),x=x,u=u,stress=stress,strain=strain)
        end do
    
    end subroutine prepare_wedge_element
    
    
    
    
    subroutine extract_wedge_element(elem,curr_status,key,connec,matkey,plyangle,stress,strain,ig_point,sdv,nstep,ninc)
    
        type(wedge_element), intent(in) :: elem
        
        integer,                              optional, intent(out) :: curr_status, key, matkey
        real(dp),                             optional, intent(out) :: plyangle
        real(dp),                allocatable, optional, intent(out) :: stress(:), strain(:)
        integer,                 allocatable, optional, intent(out) :: connec(:)
        type(integration_point), allocatable, optional, intent(out) :: ig_point(:)
        type(sdv_array),         allocatable, optional, intent(out) :: sdv(:)
        integer,                        optional, intent(out) :: nstep, ninc
        
        if(present(curr_status)) curr_status=elem%curr_status
        
        if(present(key)) key=elem%key
        
        if(present(matkey)) matkey=elem%matkey
        
        if(present(plyangle)) plyangle=elem%plyangle
        
        if(present(stress)) then
            allocate(stress(nst))
            stress=elem%stress
        end if

        if(present(strain)) then
            allocate(strain(nst))
            strain=elem%strain
        end if
        
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
        
        if(present(nstep)) nstep=elem%nstep
        
        if(present(ninc)) ninc=elem%ninc
    
    
    end subroutine extract_wedge_element


    
    
    
    
    ! the integration subroutine, updates K matrix, F vector, integration point stress and strain
    ! as well as all the solution dependent variables (sdvs) at intg points and element
    subroutine integrate_wedge_element(elem,K_matrix,F_vector,nofailure)
    use toolkit_module                  ! global tools for element integration
    use lib_mat_module                  ! global material library
    use lib_node_module                 ! global node library
    use glb_clock_module                ! global analysis progress (curr. step, inc, time, dtime)
    
        type(wedge_element),intent(inout)       :: elem 
        real(kind=dp),allocatable,intent(out)   :: K_matrix(:,:), F_vector(:)
        logical,  optional,  intent(in)         :: nofailure
        
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
        logical         :: failure ! true for failure analysis
        real(kind=dp)   :: fn(nnode),dn(nnode,ndim) ! shape functions & their deriv. physical space
        real(kind=dp)   :: jac(ndim,ndim),gn(nnode,ndim),detj ! jacobian & shape func. deriv. natural space
        real(kind=dp)   :: bee(nst,ndof),beet(ndof,nst) ! b matrix and its transpose
        real(kind=dp)   :: dee(nst,nst) ! d matrix
        real(kind=dp)   :: btd(ndof,nst),btdb(ndof,ndof) ! b'*d & b'*d*b
        real(kind=dp)   :: tmpx(ndim),tmpu(ndim),tmpstrain(nst),tmpstress(nst) ! temp. x, strain & stress arrays for intg pnts      
        real(kind=dp),allocatable :: xj(:),uj(:)! nodal vars extracted from glb lib_node array
        
        integer :: i,j,kig,igstat
        
        ! variables used for calculating clength
        real(kind=dp)   :: clength,ctip(2,2)

        logical :: nofail
        
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
        failure=.false.
        
        nofail=.false.
        if(present(nofailure)) nofail=nofailure
        
        ! copy nodes from global node array 
        node(:)=lib_node(elem%connec(:))
        
        ! assign values to local arrays from nodal values
        do j=1,nnode
        
            ! extract useful values from nodes
            call extract(node(j),x=xj,u=uj)
            
            ! assign values to coords matrix and u vector
            if(allocated(xj)) then
                coords(:,j)=xj(:)
                deallocate(xj)
            else
                write(msg_file,*)'WARNING: x not allocated for node:',elem%connec(j)
            end if
            
            if(allocated(uj)) then 
                u((j-1)*ndim+1:j*ndim)=uj(1:ndim)
                deallocate(uj)
            else
                write(msg_file,*)'WARNING: u not allocated for node:',elem%connec(j)
            end if
            
        end do
        
        ! extract material values from global material array
        mat=lib_mat(elem%matkey)
        call extract(mat,matname,mattype,typekey)
        
        ! extract ply angle from element definition
        theta=elem%plyangle
        
        ! - extract curr step and inc values from glb clock module
        call extract_glb_clock(kstep=curr_step,kinc=curr_inc)
        
        ! - check if last iteration has converged, and update the current step & increment no.
        if(elem%nstep.ne.curr_step .or. elem%ninc.ne.curr_inc) then
            last_converged=.true.
            elem%nstep = curr_step
            elem%ninc = curr_inc
        end if
        
        
        
        !-----------------------------------------------------------!
        !           calculate approximate clength
        !-----------------------------------------------------------!
        ! initialize relevant variables
        clength=zero
        call elem_ctips_origin('wedge',theta,coords,edg,nedge,ctip)
        clength=sqrt((ctip(1,2)-ctip(1,1))**2+(ctip(2,2)-ctip(2,1))**2)
        
        !-----------------------------------------------------------!
        !-----------------------------------------------------------! 
        
        
        ! update ig point xi and weight
        call init_ig(igxi,igwt)
        
        ! zero elem curr status for update
        elem%curr_status=zero
          
        ! zero element stress and strain (used for output only) for update
        elem%stress=zero
        elem%strain=zero
        
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
                    
                    if(nofail) failure=.false.
                    
                    if(failure) write(msg_file,*) "WARNING: failure analysis is not yet supported for &
                    & wedge_element type isotropic material; only linear elastic stiffness matrix is integrated."
                    
                    ! calculate D matrix, update tmpstress
                    call ddsdde(lib_iso(typekey),dee,strain=tmpstrain,stress=tmpstress) 
                    
                case ('lamina')
                
                    ! check if failure analysis is needed (check if strength parameters are present)
                    call extract(lib_lamina(typekey),strength_active=failure)
                    
                    if(nofail) failure=.false.
                    
                    if(failure) then
                        ! calculate D matrix, update tmpstress and sdv
                        call ddsdde(lib_lamina(typekey),clength,dee,strain=tmpstrain,stress=tmpstress,sdv=ig_sdv(2),dfail=dfail)
                    else
                        ! calculate D matrix, update tmpstress
                        call ddsdde(lib_lamina(typekey),dee=dee,strain=tmpstrain,stress=tmpstress)
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
            !call update(elem%ig_point(kig),x=tmpx,u=tmpu,strain=tmpstrain,stress=tmpstress,sdv=ig_sdv)
            call update(elem%ig_point(kig),sdv=ig_sdv)
            
            ! update elem curr status variable
            igstat=0
            if(allocated(ig_sdv(2)%i)) igstat=ig_sdv(2)%i(1)
            elem%curr_status=max(elem%curr_status,igstat)

            ! update elem stress & strain (avg of ig point stress & strains)
            elem%stress=elem%stress+tmpstress/nig
            elem%strain=elem%strain+tmpstrain/nig
            
            deallocate(ig_sdv)
            
       	end do !-looped over all int points. ig=nig
        
        F_vector=matmul(K_matrix,u) 
        
        
        ! deallocate local dynamic arrays
        if(allocated(xj)) deallocate(xj) 
        if(allocated(uj)) deallocate(uj) 
        if(allocated(ig_sdv)) deallocate(ig_sdv)
        
        
    
    end subroutine integrate_wedge_element
    
    
    
    
    
    
    
    
    
    
    
    
    ! the rest are private subroutines
    
    
    
    
    
    

    
 
    subroutine init_ig(xi,wt)

      real(kind=dp),intent(inout) :: xi(ndim,nig),wt(nig)
	
        if (nig .eq. 2) then
            xi(1,1)= one_third
            xi(2,1)= one_third
            xi(3,1)= -root3
            xi(1,2)= one_third
            xi(2,2)= one_third
            xi(3,2)= root3
            wt = half
        else if(nig == 6) then
            xi(1,1)=one/six
            xi(2,1)=one/six
            xi(3,1)=-root3

            xi(1,2)=four/six
            xi(2,2)=one/six
            xi(3,2)=-root3

            xi(1,3)=one/six
            xi(2,3)=four/six
            xi(3,3)=-root3

            xi(1,4)=one/six
            xi(2,4)=one/six
            xi(3,4)=root3

            xi(1,5)=four/six
            xi(2,5)=one/six
            xi(3,5)=root3

            xi(1,6)=one/six
            xi(2,6)=four/six
            xi(3,6)=root3

            wt=one/six

!~        else if(nig == 9) then
!~            xi(1,1)=one/six
!~            xi(2,1)=one/six
!~            xi(3,1)=-sqrt(three/five)
!~
!~            xi(1,2)=four/six
!~            xi(2,2)=one/six
!~            xi(3,2)=-sqrt(three/five)
!~
!~            xi(1,3)=one/six
!~            xi(2,3)=four/six
!~            xi(3,3)=-sqrt(three/five)
!~
!~            xi(1,4)=one/six
!~            xi(2,4)=one/six
!~            xi(3,4)=sqrt(three/five)
!~
!~            xi(1,5)=four/six
!~            xi(2,5)=one/six
!~            xi(3,5)=sqrt(three/five)
!~
!~            xi(1,6)=one/six
!~            xi(2,6)=four/six
!~            xi(3,6)=sqrt(three/five)
!~
!~            xi(1,7)=one/six
!~            xi(2,7)=one/six
!~            xi(3,7)=zero
!~
!~            xi(1,8)=four/six
!~            xi(2,8)=one/six
!~            xi(3,8)=zero
!~
!~            xi(1,9)=one/six
!~            xi(2,9)=four/six
!~            xi(3,9)=zero
!~
!~            wt(1:6)=five/54._dp
!~            wt(7:9)=eight/54._dp

        else
            write(msg_file,*) 'no. of integration points incorrect for wedge_ig!'
            call exit_function
        end if

    end subroutine init_ig
    
    
    
    subroutine init_shape(igxi,f,df)
      
        real(kind=dp),intent(inout) :: f(nnode),df(nnode,ndim)
        real(kind=dp),intent(in) :: igxi(ndim)
        
        real(kind=dp) :: xi,eta,zeta ! local variables
        xi=zero
        eta=zero
        zeta=zero

        xi=igxi(1)
        eta=igxi(2)
        zeta=igxi(3)
        
        f(1)=half*(one-xi-eta)*(one-zeta)
        f(2)=half*xi*(one-zeta)
        f(3)=half*eta*(one-zeta)
        f(4)=half*(one-xi-eta)*(one+zeta)
        f(5)=half*xi*(one+zeta)
        f(6)=half*eta*(one+zeta)
        
        
        df(1,3) = -half*(one-xi-eta)
        df(2,3) = -half*xi
        df(3,3) = -half*eta
        df(4,3) =  half*(one-xi-eta)
        df(5,3) =  half*xi
        df(6,3) =  half*eta
        
        
        df(1,1) = -half*(one-zeta)
        df(2,1) =  half*(one-zeta)
        df(3,1) =  zero
        df(4,1) = -half*(one+zeta)
        df(5,1) =  half*(one+zeta)
        df(6,1) =  zero
        
        df(1,2) = -half*(one-zeta)
        df(2,2) =  zero
        df(3,2) =  half*(one-zeta)
        df(4,2) = -half*(one+zeta)
        df(5,2) =  zero
        df(6,2) =  half*(one+zeta)
        
        

    end subroutine init_shape
    
    
    end module wedge_element_module
