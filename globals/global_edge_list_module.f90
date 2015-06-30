!***************************************! 
!   the global library of edges         ! 
!***************************************! 
                                          
module global_edge_list_module                     
                                          
implicit none                             
save                                      
                                          
integer, allocatable :: global_edge_list(:)    

contains

  subroutine empty_global_edge_list()

    if(allocated(global_edge_list)) deallocate(global_edge_list)

  end subroutine empty_global_edge_list                                          

end module global_edge_list_module                
