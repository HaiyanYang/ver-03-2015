    include 'parameter_module.f90'
    include 'xnode_module.f90'
    include 'material_module.f90'
    include 'integration_point_module.f90'

    program testdatahiding
    use parameter_module
    use xnode_module
    use material_module
    use integration_point_module

    implicit none

    type(xnode) :: node1,node2
    real(kind=dp) :: x(3),u(3),du(3),v(3),a(3),dof(6),ddof(6)
    real(kind=dp),allocatable :: xa(:),ua(:),va(:),aa(:),dofa(:),ddofa(:)


    call empty(node1)
    call empty(node2)

    u=zero
    du=zero
    v=zero
    a=one
    dof=ten
    ddof=three

    call update(node1,x=[one,one,one],u=u,du=du,v=v,a=a,dof=dof,ddof=ddof)
    node2=node1
    call export(node2,x=xa,u=ua,v=va,a=aa,dof=dofa,ddof=ddofa)

    print*,dofa,ddofa


!    type(isotropic_type) :: lib_iso(3)
!    integer :: i
!    logical :: mstat,sstat,tstat
!    i=0
!
!    do i=1,size(lib_iso)
!        call empty(lib_iso(i))
!    end do
!
!    call update(lib_iso(1),isotropic_modulus(E=one,nu=zero))
!
!    call export(lib_iso(1),modulus_active=mstat,strength_active=sstat,toughness_active=tstat)
!    !call export(lib_iso(1),mstat,sstat,tstat) ! this doesn't work
!
!    print*,mstat,sstat,tstat


!    integer :: ndim=2, nst=3
!    type(integration_point) :: igpnt
!    real(kind=dp),allocatable :: x(:),rsdv(:)
!
!    call update(igpnt,x=[zero,one],rsdv=[one,two])
!    call export(igpnt,x=x,rsdv=rsdv)
!
!    print*,x
!    print*,rsdv
!
!    call empty(igpnt)
!
!    x=zero;rsdv=ten
!
!    call export(igpnt,x=x,rsdv=rsdv)
!
!    print*,x
!    print*,rsdv









    end program testdatahiding
