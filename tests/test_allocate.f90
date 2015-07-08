program test_allocate
type tA
  integer, allocatable :: a(:)
end type
type tB
  type(tA), allocatable :: b(:)
end type

type(tB), allocatable :: x(:), y(:)

allocate(x(5))

allocate(y(10))

y = x ! works

!allocate(y, source = x) ! doesnt work

print*, size(y)

end program