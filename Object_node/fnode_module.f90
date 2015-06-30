module fnode_module
!
!  Purpose:
!   this module defines an extended node object, which contains coordinates,
!   displacement and its increment, velocity, acceleration, and extra dof and
!   its increment, ddof; this module can be used for analysis of any dimension
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    10/04/15  B. Y. Chen            Original code
!
use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

implicit none

private

type, public :: fnode ! a node with enriched d.o.f
  private ! hide components from external operation
  ! list of type components:
  ! - x         : coordinates
  ! - u, du     : displacement vector and its increment
  ! - v, a      : velocity and acceleration
  ! - dof, ddof : additional dof and incremental dof
  real(DP), allocatable :: x(:)
  real(DP), allocatable :: u(:),   du(:)
  real(DP), allocatable :: v(:),   a(:)
  real(DP), allocatable :: dof(:), ddof(:)
end type fnode

type, public :: fnode_alloc_array
  type(fnode), allocatable :: array(:)
end type fnode_alloc_array

interface empty
  module procedure empty_fnode
end interface

interface update
  module procedure update_fnode
end interface

interface extract
  module procedure extract_fnode
end interface

interface display
  module procedure display_fnode
end interface

interface operator(+)
  module procedure plus_fnode
end interface

interface operator(-)
  module procedure minus_fnode
end interface

interface operator(*)
  module procedure ratio_fnode
end interface


public :: empty, update, extract, display, &
        & operator(+), operator(-), operator(*)




contains



  
  elemental subroutine empty_fnode(this_fnode)
  ! Purpose:
  ! to deallocate all the components of this object
  
    type(fnode), intent(inout) :: this_fnode
    
    if(allocated(this_fnode%x))     deallocate(this_fnode%x)
    if(allocated(this_fnode%u))     deallocate(this_fnode%u)
    if(allocated(this_fnode%du))    deallocate(this_fnode%du)
    if(allocated(this_fnode%v))     deallocate(this_fnode%v)
    if(allocated(this_fnode%a))     deallocate(this_fnode%a)
    if(allocated(this_fnode%dof))   deallocate(this_fnode%dof)
    if(allocated(this_fnode%ddof))  deallocate(this_fnode%ddof)

  end subroutine empty_fnode


  
  elemental function plus_fnode(fnode1, fnode2) result(fnode3)
  ! Purpose:
  ! to add two fnode objects, component by component
  
    type(fnode), intent(in) :: fnode1, fnode2
    type(fnode)             :: fnode3
    
    
    if (allocated(fnode1%x) .and. allocated(fnode2%x)) then
      if (size(fnode1%x) == size(fnode2%x)) then
        allocate(fnode3%x(size(fnode1%x)))
        fnode3%x = fnode1%x + fnode2%x
      end if
    end if


    if (allocated(fnode1%u) .and. allocated(fnode2%u)) then
      if (size(fnode1%u) == size(fnode2%u)) then
        allocate(fnode3%u(size(fnode1%u)))
        fnode3%u = fnode1%u + fnode2%u
      end if
    end if
    
    if (allocated(fnode1%du) .and. allocated(fnode2%du)) then
      if (size(fnode1%du) == size(fnode2%du)) then
        allocate(fnode3%du(size(fnode1%du)))
        fnode3%du = fnode1%du + fnode2%du
      end if
    end if
    
    if (allocated(fnode1%v) .and. allocated(fnode2%v)) then
      if (size(fnode1%v) == size(fnode2%v)) then
        allocate(fnode3%v(size(fnode1%v)))
        fnode3%v = fnode1%v + fnode2%v
      end if
    end if
    
    if (allocated(fnode1%a) .and. allocated(fnode2%a)) then
      if (size(fnode1%a) == size(fnode2%a)) then
        allocate(fnode3%a(size(fnode1%a)))
        fnode3%a = fnode1%a + fnode2%a
      end if
    end if
    
    if (allocated(fnode1%dof) .and. allocated(fnode2%dof)) then
      if (size(fnode1%dof) == size(fnode2%dof)) then
        allocate(fnode3%dof(size(fnode1%dof)))
        fnode3%dof = fnode1%dof + fnode2%dof
      end if
    end if
    
    if (allocated(fnode1%ddof) .and. allocated(fnode2%ddof)) then
      if (size(fnode1%ddof) == size(fnode2%ddof)) then
        allocate(fnode3%ddof(size(fnode1%ddof)))
        fnode3%ddof = fnode1%ddof + fnode2%ddof
      end if
    end if
    
    
  end function plus_fnode


  
  elemental function minus_fnode(fnode1, fnode2) result(fnode3)
  ! Purpose:
  ! to subtract two fnode objects (second from first), component by component
  
    type(fnode), intent(in) :: fnode1, fnode2
    type(fnode)             :: fnode3
    
    
    if (allocated(fnode1%x) .and. allocated(fnode2%x)) then
      if (size(fnode1%x) == size(fnode2%x)) then
        allocate(fnode3%x(size(fnode1%x)))
        fnode3%x = fnode1%x - fnode2%x
      end if
    end if


    if (allocated(fnode1%u) .and. allocated(fnode2%u)) then
      if (size(fnode1%u) == size(fnode2%u)) then
        allocate(fnode3%u(size(fnode1%u)))
        fnode3%u = fnode1%u - fnode2%u
      end if
    end if
    

    if (allocated(fnode1%du) .and. allocated(fnode2%du)) then
      if (size(fnode1%du) == size(fnode2%du)) then
        allocate(fnode3%du(size(fnode1%du)))
        fnode3%du = fnode1%du - fnode2%du
      end if
    end if
    
    if (allocated(fnode1%v) .and. allocated(fnode2%v)) then
      if (size(fnode1%v) == size(fnode2%v)) then
        allocate(fnode3%v(size(fnode1%v)))
        fnode3%v = fnode1%v - fnode2%v
      end if
    end if
    
    if (allocated(fnode1%a) .and. allocated(fnode2%a)) then
      if (size(fnode1%a) == size(fnode2%a)) then
        allocate(fnode3%a(size(fnode1%a)))
        fnode3%a = fnode1%a - fnode2%a
      end if
    end if
    
    if (allocated(fnode1%dof) .and. allocated(fnode2%dof)) then
      if (size(fnode1%dof) == size(fnode2%dof)) then
        allocate(fnode3%dof(size(fnode1%dof)))
        fnode3%dof = fnode1%dof - fnode2%dof
      end if
    end if
    
    if (allocated(fnode1%ddof) .and. allocated(fnode2%ddof)) then
      if (size(fnode1%ddof) == size(fnode2%ddof)) then
        allocate(fnode3%ddof(size(fnode1%ddof)))
        fnode3%ddof = fnode1%ddof - fnode2%ddof
      end if
    end if
    
    
  end function minus_fnode


  
  elemental function ratio_fnode(r, fnode1) result(fnode3)
  ! Purpose:
  ! to multiply all the components of an fnode object by a constant
  
    type(fnode), intent(in) :: fnode1
    real(DP),    intent(in) :: r
    type(fnode)             :: fnode3
    
    if (allocated(fnode1%x)) then
      allocate(fnode3%x(size(fnode1%x)))
      fnode3%x = r * fnode1%x
    end if


    if (allocated(fnode1%u)) then
      allocate(fnode3%u(size(fnode1%u)))
      fnode3%u = r * fnode1%u
    end if
    
    if (allocated(fnode1%du)) then
      allocate(fnode3%du(size(fnode1%du)))
      fnode3%du = r * fnode1%du
    end if
    
    if (allocated(fnode1%v)) then
      allocate(fnode3%v(size(fnode1%v)))
      fnode3%v = r * fnode1%v
    end if
    
    if (allocated(fnode1%a)) then
      allocate(fnode3%a(size(fnode1%a)))
      fnode3%a = r * fnode1%a
    end if
    
    if (allocated(fnode1%dof)) then
      allocate(fnode3%dof(size(fnode1%dof)))
      fnode3%dof = r * fnode1%dof
    end if
    
    if (allocated(fnode1%ddof)) then
      allocate(fnode3%ddof(size(fnode1%ddof)))
      fnode3%ddof = r * fnode1%ddof
    end if
    
    
  end function ratio_fnode



  pure subroutine update_fnode (this_fnode, istat, emsg, x, u, du, v, a, &
  & dof, ddof)
  ! Purpose:
  ! to update the components of this fnode; it is used both before and during
  ! analysis to set the initial component values and update the runtime 
  ! component values, respectively. status and error messages are needed to 
  ! flag an error when inputs are not valid
  
    type(fnode),              intent(inout) :: this_fnode
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    real(DP),       optional, intent(in)    :: x(:), u(:), du(:), v(:), a(:)
    real(DP),       optional, intent(in)    :: dof(:), ddof(:)
    
    ! initialize intent(out) & local variables
    istat = STAT_SUCCESS  ! default
    emsg  = ''
    
    if (present(x)) then
      if (allocated(this_fnode%x)) then
        if (size(x) == size(this_fnode%x)) then
          this_fnode%x = x
        else
          deallocate(this_fnode%x)
          allocate(this_fnode%x(size(x)))
          this_fnode%x = x
        end if         
      else
        allocate(this_fnode%x(size(x)))
        this_fnode%x = x           
      end if
    end if
    
    ! u can exist independent of x: floating node
    if (present(u)) then
      if (allocated(this_fnode%u)) then
        if (size(u) == size(this_fnode%u)) then
          this_fnode%u = u
        else
          deallocate(this_fnode%u)
          allocate(this_fnode%u(size(u)))
          this_fnode%u = u
        end if         
      else
        allocate(this_fnode%u(size(u)))
        this_fnode%u = u           
      end if
    end if
    
    ! du cannot exist without u
    if(present(du)) then
      ! check if du dimension matches u dimension
      if(.not.allocated(this_fnode%u)) then
        emsg = "u must be allocated before du in update_fnode, fnode module"
        istat = STAT_FAILURE
        return
      else if(size(du)/=size(this_fnode%u)) then
        emsg = "du size must be consistent with u size in update_fnode,&
        &fnode module"
        istat = STAT_FAILURE
        return
      end if           
      ! update du values
      if(allocated(this_fnode%du)) then
        if(size(du)==size(this_fnode%du)) then
          this_fnode%du=du
        else
          deallocate(this_fnode%du)
          allocate(this_fnode%du(size(du)))
          this_fnode%du=du
        end if         
      else
        allocate(this_fnode%du(size(du)))
        this_fnode%du=du
      end if
    end if
    
    ! v cannot exist without u
    if(present(v)) then
      ! check if v dimension matches u dimension
      if(.not.allocated(this_fnode%u)) then
        emsg = "u must be allocated before v in update_fnode, fnode module"
        istat = STAT_FAILURE
        return
      else if(size(v)/=size(this_fnode%u)) then
        emsg = "v size must be consistent with u size in update_fnode, &
        &fnode module"
        istat = STAT_FAILURE
        return
      end if           
      ! update v values
      if(allocated(this_fnode%v)) then
        if(size(v)==size(this_fnode%v)) then
          this_fnode%v=v
        else
          deallocate(this_fnode%v)
          allocate(this_fnode%v(size(v)))
          this_fnode%v=v
        end if         
      else
        allocate(this_fnode%v(size(v)))
        this_fnode%v=v
      end if
    end if
    
    ! a cannot exist without v
    if(present(a)) then
      ! check if a dimension matches v dimension
      if(.not.allocated(this_fnode%v)) then
        emsg = "v must be allocated before a in update_fnode, fnode module"
        istat = STAT_FAILURE
        return
      else if(size(a)/=size(this_fnode%v)) then
        emsg = "a size must be consistent with v size in update_fnode, &
        &fnode module"
        istat = STAT_FAILURE
        return
      end if           
      ! update a values
      if(allocated(this_fnode%a)) then
        if(size(a)==size(this_fnode%a)) then
          this_fnode%a=a
        else
          deallocate(this_fnode%a)
          allocate(this_fnode%a(size(a)))
          this_fnode%a=a
        end if         
      else
        allocate(this_fnode%a(size(a)))
        this_fnode%a=a
      end if
    end if

    ! extra dof can exist without x: floating node extra dof
    if(present(dof)) then
      if(allocated(this_fnode%dof)) then
        if(size(dof)==size(this_fnode%dof)) then
          this_fnode%dof=dof
        else
          deallocate(this_fnode%dof)
          allocate(this_fnode%dof(size(dof)))
          this_fnode%dof=dof
        end if         
      else
        allocate(this_fnode%dof(size(dof)))
        this_fnode%dof=dof           
      end if
    end if
    
    ! ddof cannot exist without dof
    if(present(ddof)) then
      ! check if ddof dimension matches dof dimension
      if(.not.allocated(this_fnode%dof)) then
        emsg = "dof must be allocated before ddof in update_fnode, fnode module"
        istat = STAT_FAILURE
        return
      else if(size(ddof)/=size(this_fnode%dof)) then
        emsg = "ddof size must be consistent with dof size in update_fnode, &
        &fnode module"
        istat = STAT_FAILURE
        return
      end if           
      ! update ddof values
      if(allocated(this_fnode%ddof)) then
        if(size(ddof)==size(this_fnode%ddof)) then
          this_fnode%ddof=ddof
        else
          deallocate(this_fnode%ddof)
          allocate(this_fnode%ddof(size(ddof)))
          this_fnode%ddof=ddof
        end if         
      else
        allocate(this_fnode%ddof(size(ddof)))
        this_fnode%ddof=ddof           
      end if
    end if

  end subroutine update_fnode 
  
  

  pure subroutine extract_fnode (this_fnode, x, u, du, v, a, dof, ddof)
  ! Purpose:
  ! to extract all the components of this fnode
  
    type(fnode),                     intent(in)  :: this_fnode
    real(DP), allocatable, optional, intent(out) :: x(:),  u(:)
    real(DP), allocatable, optional, intent(out) :: du(:), v(:), a(:)
    real(DP), allocatable, optional, intent(out) :: dof(:), ddof(:)
    
    
    if(present(x)) then
      if(allocated(this_fnode%x)) then
        allocate(x(size(this_fnode%x)))
        x=this_fnode%x
      end if
    end if
    
    if(present(u)) then
      if(allocated(this_fnode%u)) then
        allocate(u(size(this_fnode%u)))
        u=this_fnode%u
      end if
    end if
    
    if(present(du)) then
      if(allocated(this_fnode%du)) then
        allocate(du(size(this_fnode%du)))
        du=this_fnode%du
      end if
    end if
    
    if(present(v)) then
      if(allocated(this_fnode%v)) then
        allocate(v(size(this_fnode%v)))
        v=this_fnode%v
      end if
    end if
    
    if(present(a)) then
      if(allocated(this_fnode%a)) then
        allocate(a(size(this_fnode%a)))
        a=this_fnode%a
      end if
    end if
    
    if(present(dof)) then
      if(allocated(this_fnode%dof)) then
        allocate(dof(size(this_fnode%dof)))
        dof=this_fnode%dof
      end if
    end if
    
    if(present(ddof)) then
      if(allocated(this_fnode%ddof)) then
        allocate(ddof(size(this_fnode%ddof)))
        ddof=this_fnode%ddof
      end if
    end if
    

  end subroutine extract_fnode 
  
  

  
  subroutine display_fnode (this)
  ! Purpose:
  ! to display this fnode object's components on cmd window
  ! this is useful for debugging
 
    type(fnode), intent(in) :: this
    
    ! local variable to set the output format
    character(len=20) :: display_fmt
    integer :: i
    
    ! initialize local variable
    display_fmt = ''
    i = 0

    ! set display format for real
    ! ES for real (scientific notation)
    ! 10 is width, 3 is no. of digits aft decimal point
    ! note that for scientific real, ESw.d, w>=d+7
    display_fmt = '(ES10.3)' 
    
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the components of the inquired fnode object :'
    write(*,'(1X, A)') ''
    
    write(*,'(1X, A)') '- x of this node is: '
    if (allocated(this%x)) then
      do i = 1, size(this%x)
        write(*,display_fmt,advance="no") this%x(i) 
      end do
      write(*,'(1X, A)') ''
    else
      write(*,'(1X, A)') 'x of this node is not allocated'
    end if
    write(*,'(1X, A)') ''
    
    write(*,'(1X, A)') '- u of this node is: '
    if (allocated(this%u)) then
      do i = 1, size(this%u)
        write(*,display_fmt,advance="no") this%u(i)
      end do
      write(*,'(1X, A)') ''
    else
      write(*,'(1X, A)') 'u of this node is not allocated'
    end if
    write(*,'(1X, A)') ''

    write(*,'(1X, A)') '- du of this node is: '
    if (allocated(this%du)) then
      do i = 1, size(this%du)
        write(*,display_fmt,advance="no") this%du(i)
      end do
      write(*,'(1X, A)') ''
    else
      write(*,'(1X, A)') 'du of this node is not allocated'
    end if
    write(*,'(1X, A)') ''
    
    write(*,'(1X, A)') '- v of this node is: '
    if (allocated(this%v)) then
      do i = 1, size(this%v)
        write(*,display_fmt,advance="no") this%v(i)
      end do
      write(*,'(1X, A)') ''
    else
      write(*,'(1X, A)') 'v of this node is not allocated'
    end if
    write(*,'(1X, A)') ''
    
    write(*,'(1X, A)') '- a of this node is: '
    if (allocated(this%a)) then
      do i = 1, size(this%a)
        write(*,display_fmt,advance="no") this%a(i)
      end do
      write(*,'(1X, A)') ''
    else
      write(*,'(1X, A)') 'a of this node is not allocated'
    end if
    write(*,'(1X, A)') ''
    
    write(*,'(1X, A)') '- dof of this node is: '
    if (allocated(this%dof)) then
      do i = 1, size(this%dof)
        write(*,display_fmt,advance="no") this%dof(i)
      end do
      write(*,'(1X, A)') ''
    else
      write(*,'(1X, A)') 'dof of this node is not allocated'
    end if
    write(*,'(1X, A)') ''
    
    write(*,'(1X, A)') '- ddof of this node is: '
    if (allocated(this%ddof)) then
      do i = 1, size(this%ddof)
        write(*,display_fmt,advance="no") this%ddof(i)
      end do
      write(*,'(1X, A)') ''
    else
      write(*,'(1X, A)') 'ddof of this node is not allocated'
    end if
    write(*,'(1X, A)') ''

  end subroutine display_fnode
  
  




end module fnode_module