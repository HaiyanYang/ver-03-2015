subroutine fnm_elems()                                     
use parameter_module,      only: DP                        
use elem_list_module,      only: elem_list,&               
                      & elem_node_connec, elem_edge_connec 
use fBrickLam_elem_module, only: plyblock_layup, set       
                                                           
  integer :: nelem   = 0                                   
  integer :: elnnode = 0                                   
  integer :: elnedge = 0                                   
  integer :: nplyblk = 0                                   
  integer, allocatable :: nodecnc(:), edgecnc(:)           
  type(plyblock_layup), allocatable :: layup(:)            
                                                           
  nelem   =2                                
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
                                                           
 layup(1)=plyblock_layup(angle=0._DP,nplies=1) 
 layup(2)=plyblock_layup(angle=90._DP,nplies=1) 
 layup(3)=plyblock_layup(angle=0._DP,nplies=1) 
                                                           

  nodecnc=[ &
& 1,2,4,3,7,8,10,9,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,41,42,44,43,47,48,50,49,53,54,55,56,57,58,59,60,61,62, &
& 63,64,65,66,67,68,81,82,84,83,87,88,90,89,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,121,122,123,124,125, &
& 126,127,128,135,136,137,138,139,140,141,142 &
& ]
  edgecnc=[ &
& 1,2,3,4,5,6,7,8,15,16,17,18,19,20,21,22,29,30,31,32,33,34,35,36 &
& ]
  call set(elem_list(1), NPLYBLKS=3,& 
& node_connec=nodecnc, layup=layup)
  elem_node_connec(:,1)=nodecnc(:)
  elem_edge_connec(:,1)=edgecnc(:)


  nodecnc=[ &
& 3,4,6,5,9,10,12,11,20,19,29,30,31,32,33,34,28,27,35,36,37,38,39,40,43,44,46,45,49,50,52,51,60,59,69,70,71,72,73,74,68, &
& 67,75,76,77,78,79,80,83,84,86,85,89,90,92,91,100,99,109,110,111,112,113,114,108,107,115,116,117,118,119,120,117,129,130, &
& 131,113,132,133,134,131,143,144,145,127,146,147,148 &
& ]
  edgecnc=[ &
& 3,9,10,11,7,12,13,14,17,23,24,25,21,26,27,28,31,37,38,39,35,40,41,42 &
& ]
  call set(elem_list(2), NPLYBLKS=3,& 
& node_connec=nodecnc, layup=layup)
  elem_node_connec(:,2)=nodecnc(:)
  elem_edge_connec(:,2)=edgecnc(:)

end subroutine fnm_elems
