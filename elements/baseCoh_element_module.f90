module baseCoh_element_module
    use parameter_module
    use xnode_module
    use coh3d6_element_module
    use coh3d8_element_module

    implicit none
    private
    
    integer, parameter :: ndim=3

    type,public :: baseCoh_element
        private ! encapsulate components of this type
        
        integer                         :: fstat=0
        
        character(len=eltypelength)     :: eltype=''    ! can be all types of 2D elements
        integer                         :: ID_matlist=0     ! material index in glb material array
        
        integer,allocatable             :: connec(:)    ! sub_elem connec to global node library
        
        type(coh3d6_element),allocatable :: coh3d6     ! 3D coh. sub elements
        type(coh3d8_element),allocatable :: coh3d8     ! 3D coh. sub elements
        
        real(kind=dp),allocatable       :: Tmatrix(:,:) ! interpolation matrix
        type(xnode),  allocatable       :: mnode(:)     ! interpolated (material domain) nodes
    end type

    interface empty
        module procedure empty_baseCoh_element
    end interface
  
    interface prepare
        module procedure prepare_baseCoh_element
    end interface
    
    interface update
        module procedure update_baseCoh_element
    end interface
    
    interface integrate
        module procedure integrate_baseCoh_element
    end interface
    
    interface extract
        module procedure extract_baseCoh_element
    end interface




    public :: empty,prepare,update,integrate,extract
    
    
    
    
    
    
    
    
    contains
    
    
    
    
    
    
    
    
    subroutine empty_baseCoh_element(elem)
    
        type(baseCoh_element), intent(out) :: elem
        
        elem%fstat=0
        elem%eltype=''    ! can be all types of 2D elements
        elem%ID_matlist=0     ! material index in glb material array
        elem%plyangle=zero

        
        if(allocated(elem%connec))      deallocate(elem%connec)
        if(allocated(elem%wedge))       deallocate(elem%wedge)
        if(allocated(elem%brick))       deallocate(elem%brick)
        if(allocated(elem%coh3d6))      deallocate(elem%coh3d6)
        if(allocated(elem%coh3d8))      deallocate(elem%coh3d8)
        if(allocated(elem%Tmatrix))     deallocate(elem%Tmatrix)
        if(allocated(elem%mnode))       deallocate(elem%mnode)
     
        
    end subroutine empty_baseCoh_element
    
    
    
    
    
    
    
    
    
    subroutine prepare_baseCoh_element(elem,eltype,ID_matlist,plyangle,connec,Tmatrix,mnode)
    
        type(baseCoh_element), intent(inout) :: elem
        
        character(len=*),intent(in)             :: eltype       ! can be all types of 2D elements
        integer,intent(in)                      :: ID_matlist       ! material index in glb material array
        real(dp),intent(in),optional            :: plyangle
        
        integer,intent(in)                      :: connec(:)    ! sub_elem connec to global node library
        

        real(kind=dp),intent(in),optional       :: Tmatrix(:,:) ! interpolation matrix
        type(xnode),  intent(in),optional       :: mnode(:)     ! interpolated (material domain) nodes
        
        ! local variables
        integer :: lrow,lcol,nrow,ncol
        
        lrow=0; lcol=0; nrow=0; ncol=0
                   
        
        elem%eltype=eltype   
        elem%ID_matlist=ID_matlist    
        
        if(present(plyangle)) elem%plyangle=plyangle
         
        if(allocated(elem%connec))  then
            if(size(elem%connec)/=size(connec)) then
                deallocate(elem%connec)
                allocate(elem%connec(size(connec)))
            end if
        else
            allocate(elem%connec(size(connec)))
        end if
        elem%connec=connec
            


        
        if(present(Tmatrix)) then
            if(.not.present(mnode)) then
                write(msg_file,*)'mnode need to be passed in tgt with Tmatrix in preparing baseCoh elem'
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
                write(msg_file,*)'Tmatrix need to be passed in tgt with mnode in preparing baseCoh elem'
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
        
    end subroutine prepare_baseCoh_element
    
    
    
    
    subroutine update_baseCoh_element(elem,mnode)
    
      type(baseCoh_element), intent(inout) :: elem
      type(xnode), intent(in)	:: mnode(:)
      
      if(allocated(elem%mnode)) then
        if(size(elem%mnode)==size(mnode)) then
          elem%mnode=mnode
        else
          write(msg_file,*)'mis-match in baseCoh mnode array size! No mnode update done'
        end if
      else
        allocate(elem%mnode(size(mnode)))
        elem%mnode=mnode
      end if
    
    end subroutine update_baseCoh_element
    
    
    
    
    
    subroutine extract_baseCoh_element(elem,eltype,fstat,ID_matlist,plyangle,connec,wedge,brick,coh3d6,coh3d8,Tmatrix,mnode)
    
        type(baseCoh_element), intent(in) :: elem
        
        character(len=eltypelength),  intent(out),optional    :: eltype       ! can be all types of 2D elements
        integer,                    intent(out),optional    :: fstat,ID_matlist       ! material index in glb material array
        real(dp),                   intent(out),optional    :: plyangle
        
        integer,        allocatable,intent(out),optional    :: connec(:)    ! sub_elem connec to global node library
                    
        type(wedge_element), allocatable,intent(out),optional:: wedge(:)       ! wedge sub elements
        type(brick_element), allocatable,intent(out),optional:: brick(:)      ! brick sub elements
        type(coh3d6_element),allocatable,intent(out),optional:: coh3d6(:)     ! 3D coh. sub elements
        type(coh3d8_element),allocatable,intent(out),optional:: coh3d8(:)     ! 3D coh. sub elements
        
        
        real(kind=dp),  allocatable,intent(out),optional    :: Tmatrix(:,:) ! interpolation matrix
        type(xnode),    allocatable,intent(out),optional    :: mnode(:)     ! interpolated (material domain) nodes

                   
        
        if(present(eltype)) eltype=elem%eltype
        if(present(fstat)) fstat=elem%fstat
        if(present(ID_matlist)) ID_matlist=elem%ID_matlist
        if(present(plyangle)) plyangle=elem%plyangle
        
        
        if(present(connec)) then
            if(allocated(elem%connec))  then
                allocate(connec(size(elem%connec)))
                connec=elem%connec
            end if
        end if    
            
            
        
        
        if(present(wedge)) then
            if(allocated(elem%wedge))  then
                allocate(wedge(size(elem%wedge)))
                wedge=elem%wedge
            end if
        end if
        
        
        if(present(brick)) then
            if(allocated(elem%brick))  then
                allocate(brick(size(elem%brick)))
                brick=elem%brick
            end if
        end if
        
        
        if(present(coh3d6)) then
            if(allocated(elem%coh3d6))  then
                allocate(coh3d6(size(elem%coh3d6)))
                coh3d6=elem%coh3d6
            end if
        end if
        
        
        if(present(coh3d8)) then
            if(allocated(elem%coh3d8))  then
                allocate(coh3d8(size(elem%coh3d8)))
                coh3d8=elem%coh3d8
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
        
    end subroutine extract_baseCoh_element
    
    
    
    
    
    
    
    
    
    
    
    
    subroutine integrate_baseCoh_element(elem,Kmatrix,Fvector,nofailure,cohgauss)
    
        type(baseCoh_element), intent(inout) :: elem
        real(dp), allocatable, intent(out) :: Kmatrix(:,:),Fvector(:)
        logical,  optional,  intent(in)    :: nofailure
        logical,  optional,  intent(in)    :: cohgauss
        
        ! local variables
        !real(dp), allocatable ::  xi(:), xe(:), ue(:), xm(:), um(:), coords(:,:)
        integer :: i,j,l, elstat
        logical :: nofail, gauss
        
        ! temp K and F arrays for coh sub elem integration; Tmatrix relating mat
        ! nodes' dof array and num nodes dof array
        real(dp), allocatable :: Kmatrix2(:,:), Fvector2(:), Tmatrixfull(:,:)
        
        i=0; j=0; l=0; elstat=0
        nofail=.false. ; gauss=.false.
        
        if(present(nofailure)) nofail=nofailure
        if(present(cohgauss)) gauss=cohgauss
        
        select case(elem%eltype)
            case('wedge')
                ! first call, prepare element
                if(.not.allocated(elem%wedge)) then
                    allocate(elem%wedge(1))
                    call empty(elem%wedge(1))
                    call prepare(elem%wedge(1),key=0,connec=elem%connec,ID_matlist=elem%ID_matlist,plyangle=elem%plyangle)
                end if
                
                call integrate(elem%wedge(1),Kmatrix,Fvector,nofail)
                
                call extract(elem%wedge(1),fstat=elstat)
                elem%fstat=elstat
                
            case('brick')
                if(.not.allocated(elem%brick)) then
                    allocate(elem%brick(1))
                    call empty(elem%brick(1))
                    call prepare(elem%brick(1),key=0,connec=elem%connec,ID_matlist=elem%ID_matlist,plyangle=elem%plyangle)
                end if
                
                call integrate(elem%brick(1),Kmatrix,Fvector,nofail)
                
                call extract(elem%brick(1),fstat=elstat)
                elem%fstat=elstat
                
            case('coh3d6')
                if(.not.allocated(elem%coh3d6)) then 
                    allocate(elem%coh3d6(1))
                    call empty(elem%coh3d6(1))
                    ! check if the coh domain is an interpolated domain
                    if(allocated(elem%Tmatrix)) then
                        ! connec array is not needed any more; zero it instead
                        call prepare(elem%coh3d6(1),key=0,connec=[0,0,0,0,0,0],ID_matlist=elem%ID_matlist)
                        ! material nodes need to be passed into the integration subroutine
                        if(.not.allocated(elem%mnode)) then
                            write(msg_file,*)'interoplated material nodes must be allocated in baseCoh element module!'
                            call exit_function
                        end if
                    else
                        call prepare(elem%coh3d6(1),key=0,connec=elem%connec,ID_matlist=elem%ID_matlist)
                    end if
                end if
                
                if(allocated(elem%Tmatrix)) then
                
                    call integrate(elem%coh3d6(1),Kmatrix2,Fvector2,nofail,gauss,elem%mnode)
                    
                    if(allocated(Tmatrixfull)) deallocate(Tmatrixfull)
                    allocate(Tmatrixfull(6*ndim,7*ndim))  ! 6 mat nodal dofs and 7 num nodal dofs
                    Tmatrixfull=zero
                    ! Tmatrix corresponding to bottom nodes (3 mat nodes and 4 num nodes)
                    do i=1, 3    
                        do j=1, 4
                            do l=1, ndim
                                Tmatrixfull((i-1)*ndim+l,(j-1)*ndim+l)=elem%Tmatrix(i,j)
                            end do
                        end do   
                    end do
                    ! Tmatrix corresponding to top 3 nodes
                    do i=4, 6    
                        do l=1, ndim
                            Tmatrixfull((i-1)*ndim+l,i*ndim+l)=one
                        end do  
                    end do
                    
                    if(allocated(Kmatrix)) deallocate(Kmatrix)
                    allocate(Kmatrix(7*ndim,7*ndim)); Kmatrix=zero
                    
                    if(allocated(Fvector)) deallocate(Fvector)
                    allocate(Fvector(7*ndim)); Fvector=zero
                    
                    Kmatrix=matmul(matmul(transpose(Tmatrixfull),Kmatrix2),Tmatrixfull)
                    Fvector=matmul(transpose(Tmatrixfull),Fvector2)
                    
                else
                    call integrate(elem%coh3d6(1),Kmatrix,Fvector,nofail,gauss)
                end if
                
                call extract(elem%coh3d6(1),fstat=elstat)
                elem%fstat=elstat
                
            case('coh3d8')
                if(.not.allocated(elem%coh3d8)) then 
                    allocate(elem%coh3d8(1))
                    call empty(elem%coh3d8(1))
                    ! check if the coh domain is an interpolated domain
                    if(allocated(elem%Tmatrix)) then
                        ! connec array is not needed any more; zero it instead
                        call prepare(elem%coh3d8(1),key=0,connec=[0,0,0,0,0,0,0,0],ID_matlist=elem%ID_matlist)
                        ! material nodes need to be passed into the integration subroutine
                        if(.not.allocated(elem%mnode)) then
                            write(msg_file,*)'interoplated material nodes must be allocated in baseCoh element module!'
                            call exit_function
                        end if
                    else
                        call prepare(elem%coh3d8(1),key=0,connec=elem%connec,ID_matlist=elem%ID_matlist)
                    end if
                end if
                
                if(allocated(elem%Tmatrix)) then
                
                    call integrate(elem%coh3d8(1),Kmatrix2,Fvector2,nofail,gauss,elem%mnode)
                    
                    if(allocated(Tmatrixfull)) deallocate(Tmatrixfull)
                    allocate(Tmatrixfull(8*ndim,8*ndim))  ! 8 mat nodal dofs and 8 num nodal dofs
                    Tmatrixfull=zero
                    ! Tmatrix corresponding to bottom nodes (4 mat nodes and 4 num nodes)
                    do i=1, 4
                        do j=1, 4
                            do l=1, ndim
                                Tmatrixfull((i-1)*ndim+l,(j-1)*ndim+l)=elem%Tmatrix(i,j)
                            end do
                        end do   
                    end do
                    ! Tmatrix corresponding to top 4 nodes
                    do i=5, 8    
                        do l=1, ndim
                            Tmatrixfull((i-1)*ndim+l,(i-1)*ndim+l)=one
                        end do  
                    end do
                    
                    if(allocated(Kmatrix)) deallocate(Kmatrix)
                    allocate(Kmatrix(8*ndim,8*ndim)); Kmatrix=zero
                    
                    if(allocated(Fvector)) deallocate(Fvector)
                    allocate(Fvector(8*ndim)); Fvector=zero
                    
                    Kmatrix=matmul(matmul(transpose(Tmatrixfull),Kmatrix2),Tmatrixfull)
                    Fvector=matmul(transpose(Tmatrixfull),Fvector2)
                    
                else
                    call integrate(elem%coh3d8(1),Kmatrix,Fvector,nofail,gauss)
                end if
                
                call extract(elem%coh3d8(1),fstat=elstat)
                elem%fstat=elstat
                
                
            case default
                write(msg_file,*)'unsupported elem type in baseCoh element module!'
                call exit_function
        end select
        
        if(allocated(Kmatrix2)) deallocate(Kmatrix2)
        if(allocated(Fvector2)) deallocate(Fvector2)
        if(allocated(Tmatrixfull)) deallocate(Tmatrixfull)
    
    end subroutine integrate_baseCoh_element
    
    
    

end module baseCoh_element_module
