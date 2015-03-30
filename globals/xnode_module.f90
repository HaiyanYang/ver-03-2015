    module xnode_module
    use parameter_module
    
      implicit none
      
      private
    
      type, public :: xnode ! a node with enriched d.o.f
      
        private ! hide components from external operation
        
        integer :: key=0
        real(kind=dp),allocatable :: x(:) ! coordinates, density
        real(kind=dp),allocatable :: u(:), du(:), v(:), a(:) ! displacements, incremental disp. velocity and acceleration
        real(kind=dp),allocatable :: dof(:), ddof(:) ! additional dof and incremental dof
        
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
      

      public :: empty,update,extract,operator(+),operator(-),operator(*)



      
      contains





      
      ! empty a node
      subroutine empty_xnode(this_xnode)
      
      	type(xnode),intent(out) :: this_xnode
      	
        this_xnode%key=0
        if(allocated(this_xnode%x)) deallocate(this_xnode%x)
        if(allocated(this_xnode%u)) deallocate(this_xnode%u)
        if(allocated(this_xnode%du)) deallocate(this_xnode%du)
        if(allocated(this_xnode%v)) deallocate(this_xnode%v)
        if(allocated(this_xnode%a)) deallocate(this_xnode%a)
        if(allocated(this_xnode%dof)) deallocate(this_xnode%dof)
        if(allocated(this_xnode%ddof)) deallocate(this_xnode%ddof)

      end subroutine empty_xnode


 
      
      function plus_xnode(xnode1,xnode2) result(xnode3)
      
      	type(xnode),intent(in) :: xnode1,xnode2
        type(xnode)            :: xnode3
      	
        !xnode3%key=xnode1%key+xnode2%key
        
        if(allocated(xnode1%x) .and. allocated(xnode2%x)) then
            if(size(xnode1%x)==size(xnode2%x)) then
                allocate(xnode3%x(size(xnode1%x)))
                xnode3%x=xnode1%x+xnode2%x
            else
                write(msg_file,*)'nodes x dimension not consistent in plus_xnode'
                call exit_function
            end if
        else
            write(msg_file,*)'one node x not allocated in plus_xnode'
            call exit_function
        end if
 
 
        if(allocated(xnode1%u) .and. allocated(xnode2%u)) then
            if(size(xnode1%u)==size(xnode2%u)) then
                allocate(xnode3%u(size(xnode1%u)))
                xnode3%u=xnode1%u+xnode2%u
            else
                write(msg_file,*)'nodes u dimension not consistent in plus_xnode'
                call exit_function
            end if
        else
            write(msg_file,*)'one node u not allocated in plus_xnode'
            call exit_function
        end if
        
  
        !~if(allocated(xnode%du)) 
        !~if(allocated(xnode%v)) 
        !~if(allocated(xnode%a)) 
        !~if(allocated(xnode%dof)) 
        !~if(allocated(xnode%ddof)) 
        
        
      end function plus_xnode
    




      
      function minus_xnode(xnode1,xnode2) result(xnode3)
      
      	type(xnode),intent(in) :: xnode1,xnode2
        type(xnode)            :: xnode3
      	
        !xnode3%key=xnode1%key+xnode2%key
        
        if(allocated(xnode1%x) .and. allocated(xnode2%x)) then
            if(size(xnode1%x)==size(xnode2%x)) then
                allocate(xnode3%x(size(xnode1%x)))
                xnode3%x=xnode1%x-xnode2%x
            else
                write(msg_file,*)'nodes x dimension not consistent in minus_xnode'
                call exit_function
            end if
        else
            write(msg_file,*)'one node x not allocated in minus_xnode'
            call exit_function
        end if
 
 
        if(allocated(xnode1%u) .and. allocated(xnode2%u)) then
            if(size(xnode1%u)==size(xnode2%u)) then
                allocate(xnode3%u(size(xnode1%u)))
                xnode3%u=xnode1%u-xnode2%u
            else
                write(msg_file,*)'nodes u dimension not consistent in minus_xnode'
                call exit_function
            end if
        else
            write(msg_file,*)'one node u not allocated in minus_xnode'
            call exit_function
        end if
        
  
        !~if(allocated(xnode%du)) 
        !~if(allocated(xnode%v)) 
        !~if(allocated(xnode%a)) 
        !~if(allocated(xnode%dof)) 
        !~if(allocated(xnode%ddof)) 
        
        
      end function minus_xnode
    





      
      function ratio_xnode(r,xnode1) result(xnode3)
      
      	type(xnode),intent(in) :: xnode1
        real(dp),   intent(in) :: r
        type(xnode)            :: xnode3
      	
        !xnode3%key=xnode1%key+xnode2%key
        
        if(allocated(xnode1%x)) then
            allocate(xnode3%x(size(xnode1%x)))
            xnode3%x=r*xnode1%x
        else
            write(msg_file,*)'node x not allocated in ratio_xnode'
            call exit_function
        end if
 
 
        if(allocated(xnode1%u)) then
            allocate(xnode3%u(size(xnode1%u)))
            xnode3%u=r*xnode1%u
        else
            write(msg_file,*)'one node u not allocated in ratio_xnode'
            call exit_function
        end if
        
  
        !~if(allocated(xnode%du)) 
        !~if(allocated(xnode%v)) 
        !~if(allocated(xnode%a)) 
        !~if(allocated(xnode%dof)) 
        !~if(allocated(xnode%ddof)) 
        
        
      end function ratio_xnode
    





      ! update a node
      subroutine update_xnode(this_xnode,key,x,u,du,v,a,dof,ddof)
      
      	type(xnode),intent(inout) :: this_xnode
        integer,optional,intent(in) :: key
        real(kind=dp),optional,intent(in) :: x(:),u(:),du(:),v(:),a(:)
        real(kind=dp),optional,intent(in) :: dof(:),ddof(:)
        
        if(present(key)) this_xnode%key=key
        
        if(present(x)) then
            if(allocated(this_xnode%x)) then
                if(size(x)==size(this_xnode%x)) then
                    this_xnode%x=x
                else
                    deallocate(this_xnode%x)
                    allocate(this_xnode%x(size(x)))
                    this_xnode%x=x
                end if         
            else
                allocate(this_xnode%x(size(x)))
                this_xnode%x=x           
            end if
        end if
        
        ! u can exist independent of x: floating node
        if(present(u)) then
            if(allocated(this_xnode%u)) then
                if(size(u)==size(this_xnode%u)) then
                    this_xnode%u=u
                else
                    deallocate(this_xnode%u)
                    allocate(this_xnode%u(size(u)))
                    this_xnode%u=u
                end if         
            else
                allocate(this_xnode%u(size(u)))
                this_xnode%u=u           
            end if
        end if
        
        ! du cannot exist without u
        if(present(du)) then
            ! check if du dimension matches u dimension
            if(.not.allocated(this_xnode%u)) then
                write(msg_file,*)"**u must be allocated before du in update_xnode**"
                call exit_function
            else if(size(du)/=size(this_xnode%u)) then
                write(msg_file,*)"**du size must be consistent with u size in update_xnode**"
                call exit_function
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
                write(msg_file,*)"**u must be allocated before v in update_xnode**"
                call exit_function
            else if(size(v)/=size(this_xnode%u)) then
                write(msg_file,*)"**v size must be consistent with u size in update_xnode**"
                call exit_function
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
                write(msg_file,*)"**v must be allocated before a in update_xnode**"
                call exit_function
            else if(size(a)/=size(this_xnode%v)) then
                write(msg_file,*)"**a size must be consistent with v size in update_xnode**"
                call exit_function
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
                write(msg_file,*)"**dof must be allocated before ddof in update_xnode**"
                call exit_function
            else if(size(ddof)/=size(this_xnode%dof)) then
                write(msg_file,*)"**ddof size must be consistent with dof size in update_xnode**"
                call exit_function
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
      
      
      
      
      
      
           ! update a node
    subroutine extract_xnode(this_xnode,key,x,u,du,v,a,dof,ddof)
      
      	type(xnode),intent(in) :: this_xnode
        integer,optional,intent(out) :: key
        real(kind=dp),allocatable,optional,intent(out) :: x(:),u(:),du(:),v(:),a(:)
        real(kind=dp),allocatable,optional,intent(out) :: dof(:),ddof(:)
        
        if(present(key)) key=this_xnode%key
        
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