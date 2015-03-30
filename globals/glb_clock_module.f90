    !***************************************!
    !   the global progress of analysis     !
    !                                       !
    !   this module is continuously updated !
    !   at the start of each increment      !
    !***************************************!
    
    module glb_clock_module
    use parameter_module
    
    implicit none
    private
    save

    integer       :: curr_step, curr_inc
    real(kind=dp) :: curr_time, curr_dtime 
    real(kind=dp) :: curr_pnewdt
    
    
    public :: extract_glb_clock, update_glb_clock, initialize_glb_clock
    
    
    
    
    contains
    
    
    
    
    subroutine initialize_glb_clock()
    
        curr_step=0
        curr_inc=0
        curr_time=zero
        curr_dtime=zero
        curr_pnewdt=one
    
    end subroutine
    
    
    
    
    subroutine extract_glb_clock(kstep,kinc,time,dtime,pnewdt)
        integer, optional, intent(out) :: kstep, kinc
        real(dp),optional, intent(out) :: time, dtime, pnewdt
        
        if(present(kstep)) kstep=curr_step
        
        if(present(kinc)) kinc=curr_inc
        
        if(present(time)) time=curr_time
        
        if(present(dtime)) dtime=curr_dtime
        
        if(present(pnewdt)) pnewdt=curr_pnewdt
    
    end subroutine
    
    
    
    subroutine update_glb_clock(kstep,kinc,time,dtime,pnewdt)
        integer, optional, intent(in) :: kstep, kinc
        real(dp),optional, intent(in) :: time, dtime ,pnewdt
        
        if(present(kstep)) curr_step=kstep
        
        if(present(kinc)) curr_inc=kinc
        
        if(present(time)) curr_time=time
        
        if(present(dtime)) curr_dtime=dtime
        
        if(present(pnewdt)) curr_pnewdt=pnewdt
    
    end subroutine
    
    
    
    
    !~private
    !~
    !~type, public :: glb_clock_type
    !~private
    !~
    !~integer       :: kstep, kinc
    !~real(kind=dp) :: time, dtime
    !~
    !~end type
    !~
    !~interface empty
    !~    module procedure empty_glb_clock
    !~end interface
    !~
    !~interface update
    !~    module procedure update_glb_clock
    !~end interface
    !~
    !~interface extract
    !~    module procedure extract_glb_clock
    !~end interface
    !~
    !~public :: empty, update, extract
    !~
    !~
    !~
    !~
    !~contains
    !~
    !~
    !~
    !~
    !~subroutine empty_glb_clock(glb_clock)
    !~
    !~type(glb_clock_type), intent(out) :: glb_clock
    !~
    !~glb_clock%kstep=0
    !~glb_clock%kinc=0
    !~glb_clock%time=zero
    !~glb_clock%dtime=zero
    !~
    !~end subroutine
    !~
    !~
    !~subroutine update_glb_clock(glb_clock,kstep,kinc,time,dtime)
    !~
    !~type(glb_clock_type), intent(inout) :: glb_clock
    !~integer, optional, intent(in)       :: kstep, kinc
    !~real(kind=dp), optional, intent(in) :: time, dtime
    !~
    !~if(present(kstep))  glb_clock%kstep=kstep
    !~
    !~if(present(kinc))   glb_clock%kinc=kinc
    !~
    !~if(present(time))   glb_clock%time=time
    !~
    !~if(present(dtime))  glb_clock%dtime=dtime
    !~
    !~end subroutine
    !~
    !~
    !~subroutine extract_glb_clock(glb_clock,kstep,kinc,time,dtime)
    !~
    !~type(glb_clock_type), intent(in) :: glb_clock
    !~integer, optional, intent(out)       :: kstep, kinc
    !~real(kind=dp), optional, intent(out) :: time, dtime
    !~
    !~if(present(kstep))  kstep=glb_clock%kstep
    !~
    !~if(present(kinc))   kinc=glb_clock%kinc
    !~
    !~if(present(time))   time=glb_clock%time
    !~
    !~if(present(dtime))  dtime=glb_clock%dtime
    !~
    !~end subroutine
    
    
    end module glb_clock_module