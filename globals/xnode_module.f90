module xnode_module
!
!  Purpose:
!   this module defines an extended node object, which contains coordinates,
!   displacement and its increment, velocity, acceleration, and extra dof and
!   its increment, ddof
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    10/04/15  B. Y. Chen            Original code
!
use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

implicit none

private

type, public :: xnode ! a node with enriched d.o.f
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
end type xnode

interface empty
  module procedure empty_xnode
end interface

interface update
  module procedure update_xnode
end interface

interface extract
  module procedure extract_xnode
end interface

interface operator(+)
  module procedure plus_xnode
end interface

interface operator(-)
  module procedure minus_xnode
end interface

interface operator(*)
  module procedure ratio_xnode
end interface


public :: empty, update, extract, operator(+), operator(-), operator(*)




contains



  
  pure subroutine empty_xnode(this_xnode)
  
    type(xnode), intent(out) :: this_xnode
    
    if(allocated(this_xnode%x))     deallocate(this_xnode%x)
    if(allocated(this_xnode%u))     deallocate(this_xnode%u)
    if(allocated(this_xnode%du))    deallocate(this_xnode%du)
    if(allocated(this_xnode%v))     deallocate(this_xnode%v)
    if(allocated(this_xnode%a))     deallocate(this_xnode%a)
    if(allocated(this_xnode%dof))   deallocate(this_xnode%dof)
    if(allocated(this_xnode%ddof))  deallocate(this_xnode%ddof)

  end subroutine empty_xnode


  
  pure function plus_xnode(xnode1, xnode2) result(xnode3)
  
    type(xnode), intent(in) :: xnode1, xnode2
    type(xnode)             :: xnode3
    
    
    if (allocated(xnode1%x) .and. allocated(xnode2%x)) then
      if (size(xnode1%x) == size(xnode2%x)) then
        allocate(xnode3%x(size(xnode1%x)))
        xnode3%x = xnode1%x + xnode2%x
      end if
    end if


    if (allocated(xnode1%u) .and. allocated(xnode2%u)) then
      if (size(xnode1%u) == size(xnode2%u)) then
        allocate(xnode3%u(size(xnode1%u)))
        xnode3%u = xnode1%u + xnode2%u
      end if
    end if
    
    if (allocated(xnode1%du) .and. allocated(xnode2%du)) then
      if (size(xnode1%du) == size(xnode2%du)) then
        allocate(xnode3%du(size(xnode1%du)))
        xnode3%du = xnode1%du + xnode2%du
      end if
    end if
    
    if (allocated(xnode1%v) .and. allocated(xnode2%v)) then
      if (size(xnode1%v) == size(xnode2%v)) then
        allocate(xnode3%v(size(xnode1%v)))
        xnode3%v = xnode1%v + xnode2%v
      end if
    end if
    
    if (allocated(xnode1%a) .and. allocated(xnode2%a)) then
      if (size(xnode1%a) == size(xnode2%a)) then
        allocate(xnode3%a(size(xnode1%a)))
        xnode3%a = xnode1%a + xnode2%a
      end if
    end if
    
    if (allocated(xnode1%dof) .and. allocated(xnode2%dof)) then
      if (size(xnode1%dof) == size(xnode2%dof)) then
        allocate(xnode3%dof(size(xnode1%dof)))
        xnode3%dof = xnode1%dof + xnode2%dof
      end if
    end if
    
    if (allocated(xnode1%ddof) .and. allocated(xnode2%ddof)) then
      if (size(xnode1%ddof) == size(xnode2%ddof)) then
        allocate(xnode3%ddof(size(xnode1%ddof)))
        xnode3%ddof = xnode1%ddof + xnode2%ddof
      end if
    end if
    
    
  end function plus_xnode


  
  pure function minus_xnode(xnode1, xnode2) result(xnode3)
  
    type(xnode), intent(in) :: xnode1, xnode2
    type(xnode)             :: xnode3
    
    
    if (allocated(xnode1%x) .and. allocated(xnode2%x)) then
      if (size(xnode1%x) == size(xnode2%x)) then
        allocate(xnode3%x(size(xnode1%x)))
        xnode3%x = xnode1%x - xnode2%x
      end if
    end if


    if (allocated(xnode1%u) .and. allocated(xnode2%u)) then
      if (size(xnode1%u) == size(xnode2%u)) then
        allocate(xnode3%u(size(xnode1%u)))
        xnode3%u = xnode1%u - xnode2%u
      end if
    end if
    

    if (allocated(xnode1%du) .and. allocated(xnode2%du)) then
      if (size(xnode1%du) == size(xnode2%du)) then
        allocate(xnode3%du(size(xnode1%du)))
        xnode3%du = xnode1%du - xnode2%du
      end if
    end if
    
    if (allocated(xnode1%v) .and. allocated(xnode2%v)) then
      if (size(xnode1%v) == size(xnode2%v)) then
        allocate(xnode3%v(size(xnode1%v)))
        xnode3%v = xnode1%v - xnode2%v
      end if
    end if
    
    if (allocated(xnode1%a) .and. allocated(xnode2%a)) then
      if (size(xnode1%a) == size(xnode2%a)) then
        allocate(xnode3%a(size(xnode1%a)))
        xnode3%a = xnode1%a - xnode2%a
      end if
    end if
    
    if (allocated(xnode1%dof) .and. allocated(xnode2%dof)) then
      if (size(xnode1%dof) == size(xnode2%dof)) then
        allocate(xnode3%dof(size(xnode1%dof)))
        xnode3%dof = xnode1%dof - xnode2%dof
      end if
    end if
    
    if (allocated(xnode1%ddof) .and. allocated(xnode2%ddof)) then
      if (size(xnode1%ddof) == size(xnode2%ddof)) then
        allocate(xnode3%ddof(size(xnode1%ddof)))
        xnode3%ddof = xnode1%ddof - xnode2%ddof
      end if
    end if
    
    
  end function minus_xnode


  
  pure function ratio_xnode(r, xnode1) result(xnode3)
  
    type(xnode), intent(in) :: xnode1
    real(DP),    intent(in) :: r
    type(xnode)             :: xnode3
    
    if (allocated(xnode1%x)) then
      allocate(xnode3%x(size(xnode1%x)))
      xnode3%x = r * xnode1%x
    end if


    if (allocated(xnode1%u)) then
      allocate(xnode3%u(size(xnode1%u)))
      xnode3%u = r * xnode1%u
    end if
    
    if (allocated(xnode1%du)) then
      allocate(xnode3%du(size(xnode1%du)))
      xnode3%du = r * xnode1%du
    end if
    
    if (allocated(xnode1%v)) then
      allocate(xnode3%v(size(xnode1%v)))
      xnode3%v = r * xnode1%v
    end if
    
    if (allocated(xnode1%a)) then
      allocate(xnode3%a(size(xnode1%a)))
      xnode3%a = r * xnode1%a
    end if
    
    if (allocated(xnode1%dof)) then
      allocate(xnode3%dof(size(xnode1%dof)))
      xnode3%dof = r * xnode1%dof
    end if
    
    if (allocated(xnode1%ddof)) then
      allocate(xnode3%ddof(size(xnode1%ddof)))
      xnode3%ddof = r * xnode1%ddof
    end if
    
    
  end function ratio_xnode



  pure subroutine update_xnode (this_xnode, istat, emsg, x, u, du, v, a, &
  & dof, ddof)
  ! Purpose:
  ! to update the components of this xnode; it is used both before and during
  ! analysis to set the initial component values and update the runtime 
  ! component values, respectively.
  
    type(xnode),              intent(inout) :: this_xnode
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    real(DP),       optional, intent(in)    :: x(:), u(:), du(:), v(:), a(:)
    real(DP),       optional, intent(in)    :: dof(:), ddof(:)
    
    ! initialize intent(out) & local variables
    istat = STAT_SUCCESS  ! default
    emsg  = ''
    
    if (present(x)) then
      if (allocated(this_xnode%x)) then
        if (size(x) == size(this_xnode%x)) then
          this_xnode%x = x
        else
          deallocate(this_xnode%x)
          allocate(this_xnode%x(size(x)))
          this_xnode%x = x
        end if         
      else
        allocate(this_xnode%x(size(x)))
        this_xnode%x = x           
      end if
    end if
    
    ! u can exist independent of x: floating node
    if (present(u)) then
      if (allocated(this_xnode%u)) then
        if (size(u) == size(this_xnode%u)) then
          this_xnode%u = u
        else
          deallocate(this_xnode%u)
          allocate(this_xnode%u(size(u)))
          this_xnode%u = u
        end if         
      else
        allocate(this_xnode%u(size(u)))
        this_xnode%u = u           
      end if
    end if
    
    ! du cannot exist without u
    if(present(du)) then
      ! check if du dimension matches u dimension
      if(.not.allocated(this_xnode%u)) then
        emsg = "**u must be allocated before du in update_xnode**"
        istat = STAT_FAILURE
        return
      else if(size(du)/=size(this_xnode%u)) then
        emsg = "**du size must be consistent with u size in update_xnode**"
        istat = STAT_FAILURE
        return
      end if           
      ! update du values
      if(allocated(this_xnode%du)) then
        if(size(du)==size(this_xnode%du)) then
          this_xnode%du=du
        else
          deallocate(this_xnode%du)
          allocate(this_xnode%du(size(du)))
          this_xnode%du=du
        end if         
      else
        allocate(this_xnode%du(size(du)))
        this_xnode%du=du
      end if
    end if
    
    ! v cannot exist without u
    if(present(v)) then
      ! check if v dimension matches u dimension
      if(.not.allocated(this_xnode%u)) then
        emsg = "**u must be allocated before v in update_xnode**"
        istat = STAT_FAILURE
        return
      else if(size(v)/=size(this_xnode%u)) then
        emsg = "**v size must be consistent with u size in update_xnode**"
        istat = STAT_FAILURE
        return
      end if           
      ! update v values
      if(allocated(this_xnode%v)) then
        if(size(v)==size(this_xnode%v)) then
          this_xnode%v=v
        else
          deallocate(this_xnode%v)
          allocate(this_xnode%v(size(v)))
          this_xnode%v=v
        end if         
      else
        allocate(this_xnode%v(size(v)))
        this_xnode%v=v
      end if
    end if
    
    ! a cannot exist without v
    if(present(a)) then
      ! check if a dimension matches v dimension
      if(.not.allocated(this_xnode%v)) then
        emsg = "**v must be allocated before a in update_xnode**"
        istat = STAT_FAILURE
        return
      else if(size(a)/=size(this_xnode%v)) then
        emsg = "**a size must be consistent with v size in update_xnode**"
        istat = STAT_FAILURE
        return
      end if           
      ! update a values
      if(allocated(this_xnode%a)) then
        if(size(a)==size(this_xnode%a)) then
          this_xnode%a=a
        else
          deallocate(this_xnode%a)
          allocate(this_xnode%a(size(a)))
          this_xnode%a=a
        end if         
      else
        allocate(this_xnode%a(size(a)))
        this_xnode%a=a
      end if
    end if

    ! extra dof can exist without x: floating node extra dof
    if(present(dof)) then
      if(allocated(this_xnode%dof)) then
        if(size(dof)==size(this_xnode%dof)) then
          this_xnode%dof=dof
        else
          deallocate(this_xnode%dof)
          allocate(this_xnode%dof(size(dof)))
          this_xnode%dof=dof
        end if         
      else
        allocate(this_xnode%dof(size(dof)))
        this_xnode%dof=dof           
      end if
    end if
    
    ! ddof cannot exist without dof
    if(present(ddof)) then
      ! check if ddof dimension matches dof dimension
      if(.not.allocated(this_xnode%dof)) then
        emsg = "**dof must be allocated before ddof in update_xnode**"
        istat = STAT_FAILURE
        return
      else if(size(ddof)/=size(this_xnode%dof)) then
        emsg = "**ddof size must be consistent with dof size in update_xnode**"
        istat = STAT_FAILURE
        return
      end if           
      ! update ddof values
      if(allocated(this_xnode%ddof)) then
        if(size(ddof)==size(this_xnode%ddof)) then
          this_xnode%ddof=ddof
        else
          deallocate(this_xnode%ddof)
          allocate(this_xnode%ddof(size(ddof)))
          this_xnode%ddof=ddof
        end if         
      else
        allocate(this_xnode%ddof(size(ddof)))
        this_xnode%ddof=ddof           
      end if
    end if

  end subroutine update_xnode 
  
  

  pure subroutine extract_xnode (this_xnode, x, u, du, v, a, dof, ddof)
  
    type(xnode),                     intent(in)  :: this_xnode
    real(DP), allocatable, optional, intent(out) :: x(:),  u(:)
    real(DP), allocatable, optional, intent(out) :: du(:), v(:), a(:)
    real(DP), allocatable, optional, intent(out) :: dof(:), ddof(:)
    
    
    if(present(x)) then
      if(allocated(this_xnode%x)) then
        allocate(x(size(this_xnode%x)))
        x=this_xnode%x
      end if
    end if
    
    if(present(u)) then
      if(allocated(this_xnode%u)) then
        allocate(u(size(this_xnode%u)))
        u=this_xnode%u
      end if
    end if
    
    if(present(du)) then
      if(allocated(this_xnode%du)) then
        allocate(du(size(this_xnode%du)))
        du=this_xnode%du
      end if
    end if
    
    if(present(v)) then
      if(allocated(this_xnode%v)) then
        allocate(v(size(this_xnode%v)))
        v=this_xnode%v
      end if
    end if
    
    if(present(a)) then
      if(allocated(this_xnode%a)) then
        allocate(a(size(this_xnode%a)))
        a=this_xnode%a
      end if
    end if
    
    if(present(dof)) then
      if(allocated(this_xnode%dof)) then
        allocate(dof(size(this_xnode%dof)))
        dof=this_xnode%dof
      end if
    end if
    
    if(present(ddof)) then
      if(allocated(this_xnode%ddof)) then
        allocate(ddof(size(this_xnode%ddof)))
        ddof=this_xnode%ddof
      end if
    end if
    

end subroutine extract_xnode 
  
  


end module xnode_module