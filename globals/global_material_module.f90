!***************************************!                  
!   the global lists of materials       !                  
!***************************************!                           
                                                           
module global_material_module                                      
use lamina_material_module,   only: lamina_material                                     
use cohesive_material_module, only: cohesive_material
                                                           
implicit none                                              
save                                                       
                                                           
type(lamina_material),   allocatable :: UDSinglePly_material 
type(cohesive_material), allocatable :: matrixCrack_material, interface_material

contains

  subroutine empty_global_material()  
                            
    if(allocated(UDSinglePly_material)) deallocate(UDSinglePly_material)
    if(allocated(matrixCrack_material)) deallocate(matrixCrack_material)
    if(allocated(interface_material))   deallocate(interface_material)

  end subroutine empty_global_material

end module global_material_module
