module precrack_module ! precrack info for all subroutines & all elements
    use parameter_module
        
    implicit none
    private


        !~integer                     :: nprc         ! no. of precracks
        !~integer,        allocatable :: prctype(:)   ! types of precracks
        !~real(kind=dp),  allocatable :: prctip(:,:)  ! tip coords of precracks

    type, public :: precrack2d
        private
        integer :: ctype=0      ! 3: weak discon. 4: cohesive 5: strong discon. 
        real(dp):: ctip(4)=zero ! (x1,y1)-(x2,y2) of precrack2d's two tips
    end type precrack2d

    
    interface create
        module procedure create_precrack2d
    end interface
      
    interface extract
        module procedure extract_precrack2d
    end interface
      

    public :: create,extract






    contains 





    subroutine create_precrack2d(prc,ctype,ctip)
    
        type(precrack2d), intent(out) :: prc
        integer,        intent(in)  :: ctype
        real(dp),       intent(in)  :: ctip(4)
    
        prc%ctype=ctype
        prc%ctip=ctip

    end subroutine create_precrack2d
    
    
    
    
    
    subroutine extract_precrack2d(prc,ctype,ctip)
    
        type(precrack2d),                   intent(in)     :: prc
        integer,                optional, intent(out)    :: ctype
        real(dp), allocatable,  optional, intent(out)    :: ctip(:)
    
        if(present(ctype)) ctype=prc%ctype
        
        if(present(ctip)) then
            allocate(ctip(4))
            ctip=prc%ctip
        end if

    end subroutine extract_precrack2d




end module precrack_module
