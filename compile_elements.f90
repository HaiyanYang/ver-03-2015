include 'globals/parameter_module.f90'
include 'globals/global_clock_module.f90'
include 'globals/global_node_list_module.f90'
include 'globals/global_material_list_module.f90'
include 'globals/global_toolkit_module.f90'
include 'elements/base_elements/brick_element_module.f90'
include 'elements/base_elements/wedge_element_module.f90'
include 'elements/base_elements/coh3d6_element_module.f90'
include 'elements/base_elements/coh3d8_element_module.f90'
include 'elements/basePly_element_module.f90'

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
!
use brick_element_module ! use everything
use wedge_element_module ! use everything
use coh3d6_element_module
use coh3d8_element_module
use basePly_element_module

implicit none



end program compile_base_elements
