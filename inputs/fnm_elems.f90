subroutine fnm_elems()                                     
use elem_list_module,      only: elem_list,&               
                      & elem_node_connec, elem_edge_connec 
use fBrickLam_elem_module, only: set                       
                                                           
  integer :: nelem   = 0                                   
  integer :: elnnode = 0                                   
  integer :: elnedge = 0                                   
                                                           
  nelem   =2                                
  elnnode =56                              
  elnedge =16                              
  allocate(elem_list(nelem))                               
  allocate(elem_node_connec(elnnode,nelem))                
  allocate(elem_edge_connec(elnedge,nelem))                
  nodecnc=[ &
& 1,2,4,3,7,8,10,9,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,41,42,44,43,47,48,50,49,53,54,55,56,57,58,59,60,61,62, &
& 63,64,65,66,67,68,81,82,83,84,85,86,87,88 &
& ]
  call set(elem_list(1), NPLYBLKS=2,& 
& node_connec=nodecnc, layup=plyblock_layup(angle=,nplies=),    & 

  nodecnc=[ &
& 3,4,6,5,9,10,12,11,20,19,29,30,31,32,33,34,28,27,35,36,37,38,39,40,43,44,46,45,49,50,52,51,60,59,69,70,71,72,73,74,68, &
& 67,75,76,77,78,79,80,77,89,90,91,73,92,93,94 &
& ]
  call set(elem_list(2), NPLYBLKS=2,& 
& node_connec=nodecnc, layup=plyblock_layup(angle=,nplies=),    & 

end subroutine fnm_elems
