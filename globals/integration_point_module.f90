    module integration_point_module
    use parameter_module
    
    implicit none
    private
    
    type, public :: integration_point ! everything needed for storage/output at ig point
        
        private
    
        !~real(kind=dp) :: x(ndim)=zero ! physical coordinates
        !~real(kind=dp) :: stress(nst)=zero ! stresses for output
        !~real(kind=dp) :: strain(nst)=zero ! strains for output
        
        real(kind=dp),allocatable :: x(:) ! physical coordinates
        real(kind=dp),allocatable :: u(:) ! displacement
        real(kind=dp),allocatable :: stress(:) ! stresses for output
        real(kind=dp),allocatable :: strain(:) ! strains for output
        
        ! below are optional terms 
        
        type(sdv_array), allocatable :: sdv(:) ! allocatable sets of sdv real, integer & logical arrays
        
        !~real(kind=dp),allocatable :: rsdv(:) ! sdvs for calculation and output
        !~integer,allocatable :: isdv(:)
        !~logical,allocatable :: lsdv(:)
        
    end type
    
    
    
    interface empty
        module procedure empty_ig
    end interface
      
    interface update
        module procedure update_ig
    end interface
      
    interface extract
        module procedure extract_ig
    end interface
      

    public :: empty,update,extract
    
    
    
    contains
    
    subroutine empty_ig(ig_point)
        type(integration_point),intent(out) :: ig_point
        
        if(allocated(ig_point%x)) deallocate(ig_point%x)
        if(allocated(ig_point%u)) deallocate(ig_point%u)
        if(allocated(ig_point%stress)) deallocate(ig_point%stress)
        if(allocated(ig_point%strain)) deallocate(ig_point%strain)
        if(allocated(ig_point%sdv)) deallocate(ig_point%sdv)
        !~if(allocated(ig_point%rsdv)) deallocate(ig_point%rsdv)
        !~if(allocated(ig_point%isdv)) deallocate(ig_point%isdv)
        !~if(allocated(ig_point%lsdv)) deallocate(ig_point%lsdv)
    
    end subroutine empty_ig
    
    
    
    !~subroutine empty_ig(ig_point)
    !~    type(integration_point),intent(out) :: ig_point
    !~    
    !~    ig_point%x=zero
    !~    ig_point%stress=zero
    !~    ig_point%strain=zero
    !~    
    !~    if(allocated(ig_point%rsdv)) then
    !~        deallocate(ig_point%rsdv)
    !~    end if
    !~    
    !~    if(allocated(ig_point%isdv)) then
    !~        deallocate(ig_point%isdv)
    !~    end if
    !~    
    !~    if(allocated(ig_point%lsdv)) then
    !~        deallocate(ig_point%lsdv)
    !~    end if
    !~
    !~end subroutine empty_ig
  
    
    subroutine update_ig(ig_point,x,u,stress,strain,sdv)
        type(integration_point),intent(inout) :: ig_point
        real(kind=dp),optional,intent(in)   :: x(:),u(:)
        real(kind=dp),optional,intent(in)   :: stress(:),strain(:)
        type(sdv_array),optional,intent(in) :: sdv(:)
        !~real(kind=dp),optional,intent(in) :: rsdv(:)
        !~integer,optional,intent(in) :: isdv(:)
        !~logical,optional,intent(in) :: lsdv(:)
        
        !~if(present(x)) ig_point%x=x
        !~
        !~if(present(stress)) ig_point%stress=stress
        !~
        !~if(present(strain)) ig_point%strain=strain
        
        
        
        if(present(x)) then        
            if(allocated(ig_point%x)) then
                if(size(x)==size(ig_point%x)) then
                    ig_point%x=x
                else
                    deallocate(ig_point%x)
                    allocate(ig_point%x(size(x)))
                    ig_point%x=x
                end if
            else
                allocate(ig_point%x(size(x)))
                ig_point%x=x
            end if
        end if 

        if(present(u)) then        
            if(allocated(ig_point%u)) then
                if(size(u)==size(ig_point%u)) then
                    ig_point%u=u
                else
                    deallocate(ig_point%u)
                    allocate(ig_point%u(size(u)))
                    ig_point%u=u
                end if
            else
                allocate(ig_point%u(size(u)))
                ig_point%u=u
            end if
        end if
        
        if(present(stress)) then        
            if(allocated(ig_point%stress)) then
                if(size(stress)==size(ig_point%stress)) then
                    ig_point%stress=stress
                else
                    deallocate(ig_point%stress)
                    allocate(ig_point%stress(size(stress)))
                    ig_point%stress=stress
                end if
            else
                allocate(ig_point%stress(size(stress)))
                ig_point%stress=stress
            end if
        end if 
        
        if(present(strain)) then        
            if(allocated(ig_point%strain)) then
                if(size(strain)==size(ig_point%strain)) then
                    ig_point%strain=strain
                else
                    deallocate(ig_point%strain)
                    allocate(ig_point%strain(size(strain)))
                    ig_point%strain=strain
                end if
            else
                allocate(ig_point%strain(size(strain)))
                ig_point%strain=strain
            end if
        end if    

        if(present(sdv)) then        
            if(allocated(ig_point%sdv)) then
                if(size(sdv)==size(ig_point%sdv)) then
                    ig_point%sdv=sdv
                else
                    deallocate(ig_point%sdv)
                    allocate(ig_point%sdv(size(sdv)))
                    ig_point%sdv=sdv
                end if
            else
                allocate(ig_point%sdv(size(sdv)))
                ig_point%sdv=sdv
            end if
        end if

    
        
        !~if(present(rsdv)) then        
        !~    if(allocated(ig_point%rsdv)) then
        !~        if(size(rsdv)==size(ig_point%rsdv)) then
        !~            ig_point%rsdv=rsdv
        !~        else
        !~            deallocate(ig_point%rsdv)
        !~            allocate(ig_point%rsdv(size(rsdv)))
        !~            ig_point%rsdv=rsdv
        !~        end if
        !~    else
        !~        allocate(ig_point%rsdv(size(rsdv)))
        !~        ig_point%rsdv=rsdv
        !~    end if
        !~end if    
        !~
        !~if(present(isdv)) then        
        !~    if(allocated(ig_point%isdv)) then
        !~        if(size(isdv)==size(ig_point%isdv)) then
        !~            ig_point%isdv=isdv
        !~        else
        !~            deallocate(ig_point%isdv)
        !~            allocate(ig_point%isdv(size(isdv)))
        !~            ig_point%isdv=isdv
        !~        end if
        !~    else
        !~        allocate(ig_point%isdv(size(isdv)))
        !~        ig_point%isdv=isdv
        !~    end if
        !~end if 
        !~
        !~if(present(lsdv)) then        
        !~    if(allocated(ig_point%lsdv)) then
        !~        if(size(lsdv)==size(ig_point%lsdv)) then
        !~            ig_point%lsdv=lsdv
        !~        else
        !~            deallocate(ig_point%lsdv)
        !~            allocate(ig_point%lsdv(size(lsdv)))
        !~            ig_point%lsdv=lsdv
        !~        end if
        !~    else
        !~        allocate(ig_point%lsdv(size(lsdv)))
        !~        ig_point%lsdv=lsdv
        !~    end if
        !~end if 
    
    end subroutine update_ig
    
    
    
    
    subroutine extract_ig(ig_point,x,u,stress,strain,sdv)
        type(integration_point),intent(in) :: ig_point
        real(kind=dp),allocatable,optional,intent(out)   :: x(:),u(:)
        real(kind=dp),allocatable,optional,intent(out)   :: stress(:),strain(:)
        type(sdv_array),allocatable,optional,intent(out) :: sdv(:)
        !~real(kind=dp),allocatable,optional,intent(out) :: rsdv(:)
        !~integer,allocatable,optional,intent(out) :: isdv(:)
        !~logical,allocatable,optional,intent(out) :: lsdv(:)
         
        
        if(present(x)) then        
            if(allocated(ig_point%x)) then
                allocate(x(size(ig_point%x)))
                x=ig_point%x
            end if
        end if
        
        if(present(u)) then        
            if(allocated(ig_point%u)) then
                allocate(u(size(ig_point%u)))
                u=ig_point%u
            end if
        end if
        
        if(present(stress)) then        
            if(allocated(ig_point%stress)) then
                allocate(stress(size(ig_point%stress)))
                stress=ig_point%stress
            end if
        end if 
        
        if(present(strain)) then        
            if(allocated(ig_point%strain)) then
                allocate(strain(size(ig_point%strain)))
                strain=ig_point%strain
            end if
        end if
                
        if(present(sdv)) then        
            if(allocated(ig_point%sdv)) then
                allocate(sdv(size(ig_point%sdv)))
                sdv=ig_point%sdv
            end if
        end if
                
                
                
        !~if(present(rsdv)) then        
        !~    if(allocated(ig_point%rsdv)) then
        !~        allocate(rsdv(size(ig_point%rsdv)))
        !~        rsdv=ig_point%rsdv
        !~    end if
        !~end if    
        !~
        !~if(present(isdv)) then        
        !~    if(allocated(ig_point%isdv)) then
        !~        allocate(isdv(size(ig_point%isdv)))
        !~        isdv=ig_point%isdv
        !~    end if
        !~end if 
        !~
        !~if(present(lsdv)) then        
        !~    if(allocated(ig_point%lsdv)) then
        !~        allocate(lsdv(size(ig_point%lsdv)))
        !~        lsdv=ig_point%lsdv
        !~    end if
        !~end if 
    
    end subroutine extract_ig
    
    
    end module integration_point_module