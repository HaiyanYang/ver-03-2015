subroutine set_fnm_nodes()            
use parameter_module, only: DP, ZERO  
use node_list_module, only: node_list 
use fnode_module,     only: update    
                                      
  integer :: nnode=0                  
  integer :: i=0                      
                                      
  nnode=24              
  allocate(node_list(nnode))          
  call update(node_list(1),        x=[-1.0_DP,-1.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(2),        x=[1.0_DP,-1.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(3),        x=[-1.0_DP,1.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(4),        x=[1.0_DP,1.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(5),        x=[-1.0_DP,-1.0_DP,0.1_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(6),        x=[1.0_DP,-1.0_DP,0.1_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(7),        x=[-1.0_DP,1.0_DP,0.1_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(8),        x=[1.0_DP,1.0_DP,0.1_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(9),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(10),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(11),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(12),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(13),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(14),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(15),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(16),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(17),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(18),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(19),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(20),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(21),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(22),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(23),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])
  call update(node_list(24),        x=[0.0_DP,0.0_DP,0.0_DP],        u=[ZERO,ZERO,ZERO])

end subroutine set_fnm_nodes
