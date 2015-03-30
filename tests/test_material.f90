! fortran test program

! include useful modules

        include 'parameter_module.f90'
        include 'material_module.f90'

        program test
        use parameter_module
        use material_module

        implicit none

        type(material) :: mat
        type(isotropic_type) :: iso(2)
        type(lamina_type) :: lamina(2)
        type(interlayer_type) :: interlayer(2)


        integer :: i

        i=0

        mat=material(matname='user1',mattype='isotropic',matkey=1)
        print*,'material info:','name=',mat%matname,'type=',mat%mattype,'key=',mat%matkey

        call update(iso(1),modulus=isotropic_modulus(E=70._dp,nu=0.2_dp))
        call update(lamina(1),modulus=lamina_modulus(E1=one,E2=three,nu12=zero,nu23=zero,G12=five,G23=four))
        call update(interlayer(1),modulus=interlayer_modulus(Dnn=one,Dtt=two,Dll=two))

        print*,'material modulus:', iso(1)%modulus, iso(1)%modulus_active
        print*,'material modulus:', lamina(1)%modulus, lamina(1)%modulus_active
        print*,'material modulus:', interlayer(1)%modulus, interlayer(1)%modulus_active

        print*,'material modulus:',  iso(1)%strength_active
        print*,'material modulus:',  lamina(1)%strength_active
        print*,'material modulus:',  interlayer(1)%strength_active


        end program test
