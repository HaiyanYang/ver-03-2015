!***************************************!
!   the global list of nodes            !
!***************************************!

module global_node_list_module
use fnode_module, only: fnode

implicit none
save

type(fnode), allocatable :: global_node_list(:)


contains

  subroutine empty_global_node_list()

    if(allocated(global_node_list)) deallocate(global_node_list)

  end subroutine empty_global_node_list

end module global_node_list_module
