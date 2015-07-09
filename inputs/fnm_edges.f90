subroutine set_fnm_edges()            
use edge_list_module, only: edge_list 
                                  
  integer :: nedge=0              
                                  
  nedge=14          
  allocate(edge_list(nedge))      
  edge_list = 0                   
                                  
end subroutine set_fnm_edges      
