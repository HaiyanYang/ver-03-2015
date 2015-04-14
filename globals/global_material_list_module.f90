!***************************************!                  
!   the global lists of materials       !                  
!***************************************!                  
                                                                      
include "materials/lamina_material_module.f90"                 
include "materials/cohesive_material_module.f90"                               
                                                           
module global_material_list_module                                      
use lamina_material_module                                     
use cohesive_material_module
                                                           
implicit none                                              
save                                                       
                                                           
type(lamina_material),   allocatable :: global_lamina_list(:)
type(cohesive_material), allocatable :: global_cohesive_list(:)   

contains

  subroutine empty_global_material_list()  
                            
    if(allocated(global_lamina_list))   deallocate(global_lamina_list)
    if(allocated(global_cohesive_list)) deallocate(global_cohesive_list)

  end subroutine empty_global_material_list

end module global_material_list_module                  
