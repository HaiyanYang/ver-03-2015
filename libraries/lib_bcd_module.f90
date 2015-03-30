    !***************************************! 
    !   the global library of bcds         ! 
    !***************************************! 
                                              
    module lib_bcd_module                    
    use parameter_module                      
                                              
    implicit none                             
    save                                      
                                              
    integer,allocatable :: lib_bcdnodes(:)    

    contains
    
    subroutine empty_lib_bcd()
    
    if(allocated(lib_bcdnodes)) deallocate(lib_bcdnodes)
    
    end subroutine empty_lib_bcd                                          

    end module lib_bcd_module                
