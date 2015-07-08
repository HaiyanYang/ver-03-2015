include 'globals/parameter_module.f90'
include 'globals/global_clock_module.f90'
include 'globals/global_toolkit_module.f90'
include 'object_materials/lamina_material_module.f90'
include 'object_materials/cohesive_material_module.f90'
include 'object_node/fnode_module.f90'
include 'object_elements/base_elements/brickPly_elem_module.f90'
include 'object_elements/base_elements/wedgePly_elem_module.f90'
include 'object_elements/base_elements/coh8Crack_elem_module.f90'
include 'object_elements/base_elements/coh6Delam_elem_module.f90'
include 'object_elements/base_elements/coh8Delam_elem_module.f90'
include 'object_elements/base_elements/abstPly_elem_module.f90'
include 'object_elements/base_elements/abstDelam_elem_module.f90'
include 'object_elements/fBrickPly_elem_module.f90'
include 'object_elements/fCoh8Delam_subelem_module.f90'
include 'object_elements/fCoh8Delam_elem_module.f90'
include 'object_elements/fBrickLam_elem_module.f90'
include 'datalists/material_list_module.f90'
include 'datalists/node_list_module.f90'
include 'datalists/edge_list_module.f90'
include 'datalists/elem_list_module.f90'
include 'inputs/input_module.f90'
include 'outputs/output_module.f90'

program compile_base_elements
! Purpose:
! to perform unit testing on global_clock_module
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    09/04/15  B. Y. Chen            Original code
!
use parameter_module
use fnode_module
use lamina_material_module
use cohesive_material_module
use brickPly_elem_module
use wedgePly_elem_module
use coh6Delam_elem_module
use coh8Delam_elem_module
use coh8Crack_elem_module
use abstPly_elem_module
use abstDelam_elem_module
use fCoh8Delam_subelem_module
use fCoh8Delam_elem_module
use fBrickPly_elem_module
use fBrickLam_elem_module
use material_list_module
use node_list_module
use edge_list_module
use elem_list_module
use input_module
use output_module

implicit none

  real(DP),    allocatable :: Kmat(:,:), Fvec(:)
  type(fnode), allocatable :: nodes(:)
  integer,     allocatable :: edge_status(:)
  integer,     allocatable :: node_cnc(:), edge_cnc(:)
  integer                  :: nnode, nedge
  integer                  :: istat, elstat
  character(len=MSGLENGTH) :: emsg
  character(len=DIRLENGTH) :: dir
  type(fBrickPly_elem), allocatable :: plyblks(:)
  type(fCoh8Delam_elem),allocatable :: interfs(:)

  nnode     = 0
  nedge     = 0
  istat     = STAT_SUCCESS
  emsg      = ''
  elstat    = 0

  call set_fnm_materials()
  call set_fnm_nodes()
  call set_fnm_edges()
  call set_fnm_elems()

  call update(node_list(3),u=[ZERO,0.5_DP,ZERO])
  call update(node_list(4),u=[ZERO,0.5_DP,ZERO])
  call update(node_list(7),u=[ZERO,0.5_DP,ZERO])
  call update(node_list(8),u=[ZERO,0.5_DP,ZERO])
  call update(node_list(13),u=[ZERO,0.5_DP,ZERO])
  call update(node_list(14),u=[ZERO,0.5_DP,ZERO])
  call update(node_list(21),u=[ZERO,0.5_DP,ZERO])
  call update(node_list(22),u=[ZERO,0.5_DP,ZERO])

  call display(node_list(22))

  nnode = size(elem_node_connec(:,1))
  nedge = size(elem_edge_connec(:,1))
  allocate(nodes(nnode),node_cnc(nnode))
  allocate(edge_status(nedge),edge_cnc(nedge))

  node_cnc(:) = elem_node_connec(:,1)
  edge_cnc(:) = elem_edge_connec(:,1)
  nodes       = node_list(node_cnc)
  edge_status = edge_list(edge_cnc)



! integrate the element
  call integrate(elem_list(1),nodes,edge_status,UDSinglePly_material, &
  &  matrixCrack_material, interface_material, Kmat, Fvec, istat, emsg)

  call extract(elem_list(1),curr_status=elstat)

!  dir='C:\Users\mpecb\Documents\GitHub\ver-03-2015\outputs\'
!  call output(kstep=1,kinc=1,outdir=dir)

  print*, istat
  print*, emsg
  print*, elstat

end program compile_base_elements
