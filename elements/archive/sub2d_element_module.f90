    module sub2d_element_module
        use parameter_module
        use xnode_module
        use tri_element_module
        use quad_element_module
        use coh2d_element_module
    
        implicit none
        private
        
        integer, parameter :: ndim=2
    
        type,public :: sub2d_element
            private ! encapsulate components of this type
            
            integer                         :: curr_status=0
            
            character(len=eltypelength)     :: eltype=''    ! can be all types of 2D elements
            integer                         :: matkey=0     ! material index in glb material array
            real(kind=dp)                   :: theta=zero   ! fibre orientation (lamina)

            integer,allocatable             :: glbcnc(:)    ! sub_elem connec to global node library
            
            type(tri_element), allocatable  :: tri(:)       ! tri sub elements
            type(quad_element),allocatable  :: quad(:)      ! quad sub elements
            type(coh2d_element),allocatable :: coh2d(:)     ! 2D coh. sub elements
            
            real(kind=dp),allocatable       :: Tmatrix(:,:) ! interpolation matrix
            type(xnode),  allocatable       :: mnode(:)     ! interpolated (material domain) nodes
        end type
    
        interface empty
            module procedure empty_sub2d_element
        end interface
      
        interface prepare
            module procedure prepare_sub2d_element
        end interface
        
        interface integrate
            module procedure integrate_sub2d_element
        end interface
        
        interface extract
            module procedure extract_sub2d_element
        end interface
    
    


        public :: empty,prepare,integrate,extract
        
        
        
        
        
        
        
        
        contains
        
        
        
        
        
        
        
        
        subroutine empty_sub2d_element(elem)
        
            type(sub2d_element), intent(out) :: elem
            
            elem%curr_status=0
            elem%eltype=''    ! can be all types of 2D elements
            elem%matkey=0     ! material index in glb material array
            elem%theta=zero   ! fibre orientation (lamina)
            
            if(allocated(elem%glbcnc))  deallocate(elem%glbcnc)
            if(allocated(elem%tri))     deallocate(elem%tri)
            if(allocated(elem%quad))    deallocate(elem%quad)
            if(allocated(elem%coh2d))   deallocate(elem%coh2d)
            if(allocated(elem%Tmatrix)) deallocate(elem%Tmatrix)
            if(allocated(elem%mnode))   deallocate(elem%mnode)
         
            
        end subroutine empty_sub2d_element
        
        
        
        
        
        
        
        
        
        subroutine prepare_sub2d_element(elem,eltype,matkey,theta,glbcnc,Tmatrix,mnode)
        
            type(sub2d_element), intent(inout) :: elem
            
            character(len=*),intent(in)    :: eltype       ! can be all types of 2D elements
            integer,intent(in)                      :: matkey       ! material index in glb material array
            
            integer,intent(in)                      :: glbcnc(:)    ! sub_elem connec to global node library
            
            real(kind=dp),intent(in),optional       :: theta        ! fibre orientation (lamina)
            real(kind=dp),intent(in),optional       :: Tmatrix(:,:) ! interpolation matrix
            type(xnode),  intent(in),optional       :: mnode(:)     ! interpolated (material domain) nodes
            
            ! local variables
            integer :: lrow,lcol,nrow,ncol
            
            lrow=0; lcol=0; nrow=0; ncol=0
                       
            
            elem%eltype=eltype   
            elem%matkey=matkey     
             
            if(allocated(elem%glbcnc))  then
                if(size(elem%glbcnc)/=size(glbcnc)) then
                    deallocate(elem%glbcnc)
                    allocate(elem%glbcnc(size(glbcnc)))
                end if
            else
                allocate(elem%glbcnc(size(glbcnc)))
            end if
            elem%glbcnc=glbcnc
                
            if(present(theta)) elem%theta=theta
            
            if(present(Tmatrix)) then
                if(.not.present(mnode)) then
                    write(msg_file,*)'mnode need to be passed in tgt with Tmatrix in preparing sub2d elem'
                    call exit_function
                end if
            
                lrow=size(Tmatrix(:,1))
                lcol=size(Tmatrix(1,:))
            
                if(allocated(elem%Tmatrix))  then
                    nrow=size(elem%Tmatrix(:,1))
                    ncol=size(elem%Tmatrix(1,:))
                    if(nrow/=lrow .or. ncol/=lcol) then
                        deallocate(elem%Tmatrix)
                        !~allocate(elem%Tmatrix,source=Tmatrix) ! not yet supported on Centos 6.5 gfortran
                        allocate(elem%Tmatrix(lrow,lcol))
                    end if
                else
                    !~allocate(elem%Tmatrix,source=Tmatrix)
                    allocate(elem%Tmatrix(lrow,lcol))
                end if 
                elem%Tmatrix=Tmatrix
            end if
            
            if(present(mnode)) then
                if(.not.present(Tmatrix)) then
                    write(msg_file,*)'Tmatrix need to be passed in tgt with mnode in preparing sub2d elem'
                    call exit_function
                end if
            
                if(allocated(elem%mnode))  then
                    if(size(elem%mnode)/=size(mnode)) then
                        deallocate(elem%mnode)
                        allocate(elem%mnode(size(mnode)))
                    end if
                else
                    allocate(elem%mnode(size(mnode)))
                end if
                elem%mnode=mnode          
            end if 
            
        end subroutine prepare_sub2d_element
        
        
        
        
        
        
        
        
        
        
        subroutine extract_sub2d_element(elem,eltype,curr_status,matkey,theta,glbcnc,tri,quad,coh2d,Tmatrix,mnode)
        
            type(sub2d_element), intent(in) :: elem
            
            character(len=eltypelength),  intent(out),optional    :: eltype       ! can be all types of 2D elements
            integer,                    intent(out),optional    :: curr_status,matkey       ! material index in glb material array
            real(kind=dp),              intent(out),optional    :: theta        ! fibre orientation (lamina)
            
            integer,        allocatable,intent(out),optional    :: glbcnc(:)    ! sub_elem connec to global node library
            
            type(tri_element),  allocatable,intent(out),optional:: tri(:)       ! tri sub elements
            type(quad_element), allocatable,intent(out),optional:: quad(:)      ! quad sub elements
            type(coh2d_element),allocatable,intent(out),optional:: coh2d(:)     ! 2D coh. sub elements
            
            
            real(kind=dp),  allocatable,intent(out),optional    :: Tmatrix(:,:) ! interpolation matrix
            type(xnode),    allocatable,intent(out),optional    :: mnode(:)     ! interpolated (material domain) nodes

                       
            
            if(present(eltype)) eltype=elem%eltype
            if(present(curr_status)) curr_status=elem%curr_status
            if(present(matkey)) matkey=elem%matkey
            if(present(theta))  theta=elem%theta
            
            
            if(present(glbcnc)) then
                if(allocated(elem%glbcnc))  then
                    allocate(glbcnc(size(elem%glbcnc)))
                    glbcnc=elem%glbcnc
                end if
            end if    
            
            
            
            if(present(tri)) then
                if(allocated(elem%tri))  then
                    allocate(tri(size(elem%tri)))
                    tri=elem%tri
                end if
            end if
            
            
            if(present(quad)) then
                if(allocated(elem%quad))  then
                    allocate(quad(size(elem%quad)))
                    quad=elem%quad
                end if
            end if
            
            
            if(present(coh2d)) then
                if(allocated(elem%coh2d))  then
                    allocate(coh2d(size(elem%coh2d)))
                    coh2d=elem%coh2d
                end if
            end if
            
            
            if(present(Tmatrix)) then           
                if(allocated(elem%Tmatrix))  then
                    allocate(Tmatrix(size(elem%Tmatrix(:,1)),size(elem%Tmatrix(1,:))))
                    Tmatrix=elem%Tmatrix
                end if     
            end if
            
            
            if(present(mnode)) then           
                if(allocated(elem%mnode))  then
                    allocate(mnode(size(elem%mnode)))
                    mnode=elem%mnode
                end if         
            end if 
            
        end subroutine extract_sub2d_element
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine integrate_sub2d_element(elem,Kmatrix,Fvector,cohgauss)
        
            type(sub2d_element), intent(inout) :: elem
            real(dp), allocatable, intent(out) :: Kmatrix(:,:),Fvector(:)
            logical,  optional,  intent(in)    :: cohgauss
            
            ! local variables
            !real(dp), allocatable ::  xi(:), xe(:), ue(:), xm(:), um(:), coords(:,:)
            integer :: i,j,l, elstat
            logical :: gauss
            
            i=0; j=0; l=0; elstat=0
            gauss=.false.
            
            if(present(cohgauss)) gauss=cohgauss
            
            select case(elem%eltype)
                case('tri')
                    ! first call, prepare element
                    if(.not.allocated(elem%tri)) then
                        allocate(elem%tri(1))
                        call empty(elem%tri(1))
                        call prepare(elem%tri(1),key=0,connec=elem%glbcnc,matkey=elem%matkey)
                    end if
                    
                    call integrate(elem%tri(1),Kmatrix,Fvector)
                    
                    call extract(elem%tri(1),curr_status=elstat)
                    elem%curr_status=elstat
                    
                case('quad')
                    if(.not.allocated(elem%quad)) then
                        allocate(elem%quad(1))
                        call empty(elem%quad(1))
                        call prepare(elem%quad(1),key=0,connec=elem%glbcnc,matkey=elem%matkey)
                    end if
                    
                    call integrate(elem%quad(1),Kmatrix,Fvector)
                    
                    call extract(elem%quad(1),curr_status=elstat)
                    elem%curr_status=elstat
                    
                case('coh2d')
                    if(.not.allocated(elem%coh2d)) then 
                        allocate(elem%coh2d(1))
                        call empty(elem%coh2d(1))
                        ! check if the coh domain is an interpolated domain
                        if(allocated(elem%Tmatrix)) then
                            ! connec array is not needed any more; zero it instead
                            call prepare(elem%coh2d(1),key=0,connec=[0,0,0,0],matkey=elem%matkey)
                            ! material nodes need to be passed into the integration subroutine
                            if(.not.allocated(elem%mnode)) then
                                write(msg_file,*)'interoplated material nodes must be allocated in sub2d element module!'
                                call exit_function
                            end if
                        else
                            call prepare(elem%coh2d(1),key=0,connec=elem%glbcnc,matkey=elem%matkey)
                        end if
                    end if
                    
                    if(allocated(elem%Tmatrix)) then
                        call integrate(elem%coh2d(1),Kmatrix,Fvector,gauss,elem%mnode)
                    else
                        call integrate(elem%coh2d(1),Kmatrix,Fvector,gauss)
                    end if
                    
                    call extract(elem%coh2d(1),curr_status=elstat)
                    elem%curr_status=elstat
                    
                case default
                    write(msg_file,*)'unsupported elem type in sub2d element module!'
                    call exit_function
            end select
        
        end subroutine integrate_sub2d_element
        
        
        
    
    end module sub2d_element_module