!***************************************!                  
!   the global library of materials     !                  
!***************************************!                  
                                                                      
include "materials/lamina_material_module.f90"                 
include "materials/cohesive_material_module.f90"              
include "materials/global_matkey_module.f90"                    
                                                           
module lib_mat_module                                      
use lamina_material_module                                     
use cohesive_material_module                                  
use global_matkey_module                                        
                                                           
implicit none                                              
save                                                       
                                                           
type(global_matkey),     allocatable :: list_matkey(:)
type(lamina_material),   allocatable :: list_lamina_mat(:)
type(cohesive_material), allocatable :: list_cohesive_mat(:)   

contains

  subroutine empty_lib_mat()  
                            
    if(allocated(list_matkey))       deallocate(list_matkey)
    if(allocated(list_lamina_mat))   deallocate(list_lamina_mat)
    if(allocated(list_cohesive_mat)) deallocate(list_cohesive_mat)

  end subroutine empty_lib_mat

end module lib_mat_module                  
