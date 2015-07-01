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

implicit none



end program compile_base_elements
