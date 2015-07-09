subroutine set_fnm_elems()                                 
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
  elnnode =24                              
  elnedge =8                              
  nplyblk =1                                
  allocate(elem_list(nelem))                               
  allocate(elem_node_connec(elnnode,nelem))                
  allocate(elem_edge_connec(elnedge,nelem))                
  allocate(nodecnc(elnnode))                               
  allocate(edgecnc(elnedge))                               
  allocate(layup(nplyblk))                                 
  nodecnc = 0                                              
  edgecnc = 0                                              
                                                           
 layup(1)=plyblock_layup(angle=90._DP,nplies=1) 
                                                           

  nodecnc=[ &
& 1,2,4,3,7,8,10,9,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 &
& ]
  edgecnc=[ &
& 1,2,3,4,5,6,7,8 &
& ]
  call set(elem_list(1), NPLYBLKS=1,& 
& node_connec=nodecnc, layup=layup)
  elem_node_connec(:,1)=nodecnc(:)
  elem_edge_connec(:,1)=edgecnc(:)


  nodecnc=[ &
& 3,4,6,5,9,10,12,11,20,19,29,30,31,32,33,34,28,27,35,36,37,38,39,40 &
& ]
  edgecnc=[ &
& 3,9,10,11,7,12,13,14 &
& ]
  call set(elem_list(2), NPLYBLKS=1,& 
& node_connec=nodecnc, layup=layup)
  elem_node_connec(:,2)=nodecnc(:)
  elem_edge_connec(:,2)=edgecnc(:)

end subroutine set_fnm_elems
