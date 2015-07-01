!***************************************!                             
!   the list of elements                                 
!***************************************!                                                        
                                                                  
module elem_list_module                                                
use fBrickLam_elem_module, only: fBrickLam_elem                                            
                                                                      
implicit none                                                         
save                                                                  
                                                                      
type(fBrickLam_elem), allocatable :: elem_list(:)
integer,              allocatable :: elem_node_connec(:,:)
integer,              allocatable :: elem_edge_connec(:,:)

contains

  subroutine empty_elem_list()

    if(allocated(elem_list))        deallocate(elem_list)
    if(allocated(elem_node_connec)) deallocate(elem_node_connec)
    if(allocated(elem_edge_connec)) deallocate(elem_edge_connec)

  end subroutine empty_elem_list
                                                    
end module elem_list_module
