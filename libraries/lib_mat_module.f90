    !***************************************!                  
    !   the global library of materials     !                  
    !***************************************!                  
                                                               
    include "materials/isotropic_type_module.f90"              
    include "materials/lamina_type_module.f90"                 
    include "materials/interface_type_module.f90"              
    include "materials/material_module.f90"                    
                                                               
    module lib_mat_module                                      
    use parameter_module                                       
    use isotropic_type_module                                  
    use lamina_type_module                                     
    use interface_type_module                                  
    use material_module                                        
                                                               
    implicit none                                              
    save                                                       
                                                               
    type(material),       allocatable    :: lib_mat(:)         
    type(isotropic_type), allocatable    :: lib_iso(:)         
    type(lamina_type),    allocatable    :: lib_lamina(:)      
    type(interface_type), allocatable    :: lib_interface(:)   
    
    contains
    
    subroutine empty_lib_mat()  
                              
    if(allocated(lib_mat)) deallocate(lib_mat)
    if(allocated(lib_iso)) deallocate(lib_iso)
    if(allocated(lib_lamina)) deallocate(lib_lamina)
    if(allocated(lib_interface)) deallocate(lib_interface)

    end subroutine empty_lib_mat
    
    end module lib_mat_module                  
