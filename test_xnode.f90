include 'globals/parameter_module.f90'
include 'globals/xnode_module.f90'

program test_xnode
! Purpose:
! to perform unit testing on xnode_module
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    11/04/15  B. Y. Chen            Original code
!
!
use parameter_module, only : DP, ZERO, ONE, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE
use xnode_module ! use everything

implicit none

integer, parameter :: ndim = 3, ndof = 2

type(xnode) :: node1, node2, node3
real(DP), allocatable :: x(:)
real(DP), allocatable :: u(:),   du(:)
real(DP), allocatable :: v(:),   a(:)
real(DP), allocatable :: dof(:), ddof(:)

integer               :: istat
character(MSGLENGTH)  :: emsg

character(len=20) :: display_fmt
character(len=10) :: cndim, cndof

! initialize local variables
! all derived types have been initialized in definition

allocate(x(ndim), u(ndim), du(ndim), v(ndim), a(ndim))
allocate(dof(ndof), ddof(ndof))

x  = ZERO
u  = ZERO
du = ZERO
v  = ZERO
a  = ZERO
dof  = ZERO
ddof = ZERO

x  = ONE
u  = ONE
du = ONE
v  = ONE
a  = ONE
dof  = ONE
ddof = ONE

istat       = STAT_SUCCESS
emsg        = ''

display_fmt = ''
cndim = ''
cndof = ''

! store ndim and ndof as strings in cndim and cndof
write(cndim,'(i5)') ndim
write(cndof,'(i5)') ndof

write(*,*) 'display node 1 before any update:'
call display (node1)

call update (node1, istat=istat, emsg=emsg, x=x, u=u)
if (istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if
write(*,*) 'display node 1 after updates on x and u:'
call display (node1)

call update (node1, istat=istat, emsg=emsg, du=du, v=v, a=a, dof=dof, ddof=ddof)
if (istat == STAT_FAILURE) then
  write(*,*) emsg
  return
end if
write(*,*) 'display node 1 after updates on all components:'
call display (node1)

node2 = 0.6_DP * node1
write(*,*) 'display node 2 as 0.6 * node1:'
call display (node2)

node3 = node1 + node2
write(*,*) 'display node 3 as a sum of node1 and node2:'
call display (node3)

node3 = node1 - 0.5_DP * node2
write(*,*) 'display node 3 as node1 - 0.5 * node2:'
call display (node3)

write(*,*) 'check extracted values of node3:'
call extract(node3, x=x, u=u, du=du, v=v, a=a, dof=dof, ddof=ddof)
display_fmt = '(1X, A,'//trim(adjustl(cndim))//'ES10.3)'
write(*,display_fmt) '- x    :', x
write(*,display_fmt) '- u    :', u
write(*,display_fmt) '- du   :', du
write(*,display_fmt) '- v    :', v
write(*,display_fmt) '- a    :', a
display_fmt = '(1X, A,'//trim(adjustl(cndof))//'ES10.3)'
write(*,display_fmt) '- dof  :', dof
write(*,display_fmt) '- ddof :', ddof

call empty (node3)
write(*,*) 'display node 3 after being emptied:'
call display (node3)



end program test_xnode
