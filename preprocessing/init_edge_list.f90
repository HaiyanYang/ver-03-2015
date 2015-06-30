    subroutine initialize_lib_edge()          
                                              
        integer :: nedge=0                    
        integer :: i=0                        
                                              
        nedge=15252   
        allocate(lib_edge(nedge))   
        lib_edge=0                  
    end subroutine initialize_lib_edge        
