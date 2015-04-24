module tst_mod

implicit none

type, public :: tst
  integer :: i
end type

end module tst_mod

program test_array_arg_indexing
use tst_mod

implicit none

integer, parameter :: e(4) = [4,8,2,3]

integer :: a(4) = [1,2,3,4]
integer :: b(8) = [1,2,3,4,5,6,7,8]
integer :: c(4), d(4)
type(tst) :: f(8)
integer :: j

call print_func(a([4,1,2,3])) ! ok
call print_func(b([4,8,2,3])) ! ok

c = [4,8,2,3]
call print_func(b(c)) ! ok

d = b(c)
call print_func(d) ! ok

call print_func(b(e)) ! ok

do j=1, 8
  f(j)%i = j
end do

write(*,*) '----------------------------------------------'
call print_tst(f([1,2,3,4])) ! ok
call print_tst(f([3,2,1,4])) ! ok
call print_tst(f(c)) ! ok
call print_tst(f(e)) ! ok

end program

subroutine print_func(a)

!integer, intent(in) :: a(4) ! ok
integer, intent(inout) :: a(4) !NOT OK!!

write(*,*) a

end subroutine

subroutine print_tst(a)
use tst_mod

type(tst), intent(in) :: a(4)    ! ok
!type(tst), intent(inout) :: a(4) !NOT OK!!

write(*,*) a

end subroutine

! test conclusion:
! fortran allows argument passing with vector subscripts, but ONLY for 
! intent(in), NOT for intent out or inout or volatile or asynchronous