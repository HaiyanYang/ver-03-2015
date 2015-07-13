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
                                                           
  nelem   =1                                
  elnnode =88                              
  elnedge =24                              
  nplyblk =3                                
  allocate(elem_list(nelem))                               
  allocate(elem_node_connec(elnnode,nelem))                
  allocate(elem_edge_connec(elnedge,nelem))                
  allocate(nodecnc(elnnode))                               
  allocate(edgecnc(elnedge))                               
  allocate(layup(nplyblk))                                 
  nodecnc = 0                                              
  edgecnc = 0                                              
                                                           
 layup(1)=plyblock_layup(angle=30._DP,nplies=1) 
 layup(2)=plyblock_layup(angle=0._DP,nplies=1) 
 layup(3)=plyblock_layup(angle=-30._DP,nplies=1) 
                                                           

  nodecnc=[ &
& 1,2,4,3,5,6,8,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,27,29,30,32,31,33,34,35,36,37,38,39,40,41,42,43, &
& 44,45,46,47,48,49,50,52,51,53,54,56,55,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83, &
& 84,85,86,87,88 &
& ]
  edgecnc=[ &
& 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 &
& ]
  call set(elem_list(1), NPLYBLKS=3)
  elem_node_connec(:,1)=nodecnc(:)
  elem_edge_connec(:,1)=edgecnc(:)

end subroutine set_fnm_elems
