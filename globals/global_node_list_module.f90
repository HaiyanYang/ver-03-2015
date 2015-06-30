!***************************************!
!   the global list of nodes            !
!***************************************!

module global_node_list_module
use xnode_module, only: xnode

implicit none
save

type(xnode), allocatable :: global_node_list(:)


contains

  subroutine empty_global_node_list()

    if(allocated(global_node_list)) deallocate(global_node_list)

  end subroutine empty_global_node_list

end module global_node_list_module
