    !***************************************! 
    !   the global library of edges         ! 
    !***************************************! 
                                              
    module lib_edge_module                    
    use parameter_module                      
                                              
    implicit none                             
    save                                      
                                              
    integer,allocatable :: lib_edge(:)    

    contains
    
    subroutine empty_lib_edge()
    
    if(allocated(lib_edge)) deallocate(lib_edge)
    
    end subroutine empty_lib_edge                                          

    end module lib_edge_module                
