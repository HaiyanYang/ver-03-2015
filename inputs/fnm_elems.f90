subroutine set_fnm_elems()                                 
use parameter_module,      only: DP                        
use elem_list_module,      only: layup, elem_list,&        
                      & elem_node_connec, elem_edge_connec 
use fBrickLam_elem_module, only: plyblock_layup, set       
                                                           
  integer :: nelem   = 0                                   
  integer :: elnnode = 0                                   
  integer :: elnedge = 0                                   
  integer :: nplyblk = 0                                   
  integer, allocatable :: nodecnc(:), edgecnc(:)           
                                                           
  nelem   =2                                
  elnnode =56                              
  elnedge =16                              
  nplyblk =2                                
  allocate(elem_list(nelem))                               
  allocate(elem_node_connec(elnnode,nelem))                
  allocate(elem_edge_connec(elnedge,nelem))                
  allocate(nodecnc(elnnode))                               
  allocate(edgecnc(elnedge))                               
  allocate(layup(nplyblk))                                 
  nodecnc = 0                                              
  edgecnc = 0                                              
                                                           
 layup(1)=plyblock_layup(angle=45._DP,nplies=1) 
 layup(2)=plyblock_layup(angle=-45._DP,nplies=1) 
                                                           

  nodecnc=[ &
& 1,2,4,3,7,8,10,9,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,41,42,44,43,47,48,50,49,53,54,55,56,57,58,59,60,61,62, &
& 63,64,65,66,67,68,81,82,83,84,85,86,87,88 &
& ]
  edgecnc=[ &
& 1,2,3,4,5,6,7,8,15,16,17,18,19,20,21,22 &
& ]
  call set(elem_list(1), NPLYBLKS=2)
  elem_node_connec(:,1)=nodecnc(:)
  elem_edge_connec(:,1)=edgecnc(:)


  nodecnc=[ &
& 3,4,6,5,9,10,12,11,18,17,29,30,31,32,33,34,26,25,35,36,37,38,39,40,43,44,46,45,49,50,52,51,58,57,69,70,71,72,73,74,66, &
& 65,75,76,77,78,79,80,83,89,90,91,87,92,93,94 &
& ]
  edgecnc=[ &
& 3,9,10,11,7,12,13,14,17,23,24,25,21,26,27,28 &
& ]
  call set(elem_list(2), NPLYBLKS=2)
  elem_node_connec(:,2)=nodecnc(:)
  elem_edge_connec(:,2)=edgecnc(:)

end subroutine set_fnm_elems
