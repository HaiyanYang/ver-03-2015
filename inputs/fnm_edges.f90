subroutine fnm_edges()            
use edge_list_module, only: edge_list 
                                  
  integer :: nedge=0              
                                  
  nedge=42          
  allocate(edge_list(nedge))      
  edge_list = 0                   
                                  
end subroutine fnm_edges
