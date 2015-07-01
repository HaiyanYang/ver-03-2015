!***************************************! 
!   the list of edges         
!***************************************! 
                                          
module edge_list_module                     
                                          
implicit none                             
save                                      
                                          
integer, allocatable :: edge_list(:)    

contains

  subroutine empty_edge_list()

    if(allocated(edge_list)) deallocate(edge_list)

  end subroutine empty_edge_list                                          

end module edge_list_module                
