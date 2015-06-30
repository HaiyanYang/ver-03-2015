!***************************************!                             
!   the global library of elements      !                             
!***************************************!                                                        
                                                                  
module global_elem_list_module                                                
use fBrickLam_element_module, only: fBrickLam_element                                            
                                                                      
implicit none                                                         
save                                                                  
                                                                      
type(fBrickLam_element), allocatable :: global_elem_list(:)
integer,                 allocatable :: elem_node_connec(:,:)
integer,                 allocatable :: elem_edge_connec(:,:)

contains

  subroutine empty_global_elem_list()

    if(allocated(global_elem_list)) deallocate(global_elem_list)
    if(allocated(elem_node_connec)) deallocate(elem_node_connec)
    if(allocated(elem_edge_connec)) deallocate(elem_edge_connec)

  end subroutine empty_global_elem_list
                                                    
end module global_elem_list_module
