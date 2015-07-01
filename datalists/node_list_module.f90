!***************************************!
!   the global list of nodes            !
!***************************************!

module node_list_module
use fnode_module, only: fnode

implicit none
save

type(fnode), allocatable :: node_list(:)


contains

  subroutine empty_node_list()

    if(allocated(node_list)) deallocate(node_list)

  end subroutine empty_node_list

end module node_list_module
