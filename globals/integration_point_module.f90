module integration_point_module
!
!  Purpose:
!  to define a generic integration point (NOT used currently)
!
!    
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    11/04/15  B. Y. Chen            Original code
!

use parameter_module, only : DP

implicit none

private


type, public :: integration_point
    
    private
    
    real(DP),allocatable :: x(:)       ! physical coordinates
    real(DP),allocatable :: u(:)       ! displacement
    real(DP),allocatable :: stress(:)  ! stresses for output
    real(DP),allocatable :: strain(:)  ! strains for output
    type(sdv_array), allocatable :: sdv(:)
    
end type integration_point



interface empty
    module procedure empty_ig
end interface empty
  
interface update
    module procedure update_ig
end interface update
  
interface extract
    module procedure extract_ig
end interface extract
  

public :: empty, update, extract



contains




pure subroutine empty_ig(ig_point)
    type(integration_point),intent(out) :: ig_point
    
    if(allocated(ig_point%x)) deallocate(ig_point%x)
    if(allocated(ig_point%u)) deallocate(ig_point%u)
    if(allocated(ig_point%stress)) deallocate(ig_point%stress)
    if(allocated(ig_point%strain)) deallocate(ig_point%strain)
    if(allocated(ig_point%sdv)) deallocate(ig_point%sdv)

end subroutine empty_ig



pure subroutine update_ig(ig_point,x,u,stress,strain,sdv)
    type(integration_point),intent(inout) :: ig_point
    real(DP),optional,intent(in)   :: x(:),u(:)
    real(DP),optional,intent(in)   :: stress(:),strain(:)
    type(sdv_array),optional,intent(in) :: sdv(:) 
    
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

end subroutine update_ig



pure subroutine extract_ig(ig_point,x,u,stress,strain,sdv)
    type(integration_point),intent(in) :: ig_point
    real(DP),allocatable,optional,intent(out)   :: x(:),u(:)
    real(DP),allocatable,optional,intent(out)   :: stress(:),strain(:)
    type(sdv_array),allocatable,optional,intent(out) :: sdv(:)
     
    
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

end subroutine extract_ig




end module integration_point_module