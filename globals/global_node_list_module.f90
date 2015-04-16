!***************************************!
!   the global list of nodes            !
!***************************************!

include "globals/xnode_module.f90"

module global_node_list_module
use xnode_module

implicit none
save

type(xnode), allocatable :: global_node_list(:)


contains

  subroutine empty_global_node_list()

    if(allocated(global_node_list)) deallocate(global_node_list)

  end subroutine empty_global_node_list

end module global_node_list_module
