    module interface_type_module
    use parameter_module
    
    implicit none
    private


    type,public :: interface_modulus
        real(kind=dp) :: Dnn, Dtt, Dll ! penalty stiffness
    end type
    
    type,public :: interface_strength
        real(kind=dp) :: tau_nc, tau_tc, tau_lc ! strengths, t and l dir. values may be different
    end type
    
    type,public :: interface_toughness
        real(kind=dp) :: Gnc, Gtc, Glc ! toughness values, t and l dir. values may be different
        real(kind=dp) :: eta ! mixed-mode law component; eta_nt/nl/tl are the same; may be different in future
    end type
    

    type,public :: interface_type
        type(interface_modulus) :: modulus
        type(interface_strength) :: strength
        type(interface_toughness) :: toughness
        logical :: modulus_active=.false.
        logical :: strength_active=.false.
        logical :: toughness_active=.false.
    end type
    

    
    interface empty
        module procedure empty_interface
    end interface
      
    interface update
        module procedure update_interface
    end interface
    
    interface ddsdde
        module procedure ddsdde_interface
    end interface
    
   
    public :: empty,update,ddsdde
  


  
    contains
    
    
    
      subroutine empty_interface(this_interface)
      
      	type(interface_type),intent(inout) :: this_interface
        
        this_interface%modulus_active=.false.
        this_interface%strength_active=.false.
        this_interface%toughness_active=.false.
      	
      end subroutine empty_interface
 


      subroutine update_interface(this_interface,modulus,strength,toughness)
      
      	type(interface_type),intent(inout) :: this_interface
        type(interface_modulus),optional,intent(in) :: modulus
        type(interface_strength),optional,intent(in) :: strength
        type(interface_toughness),optional,intent(in) :: toughness
        
        if(present(modulus)) then
                this_interface%modulus=modulus  
                this_interface%modulus_active=.true.
        end if
        
        if(present(strength)) then
                this_interface%strength=strength
                this_interface%strength_active=.true.
        end if
        
        if(present(toughness)) then
                this_interface%toughness=toughness
                this_interface%toughness_active=.true.
        end if
        

      end subroutine update_interface   
      
      
!********************************************************************************************
!************************ subroutine ddsdde ************************************************
!***** returns the tangent stiffness matrix of the material ******************************
!********************************************************************************************
      subroutine ddsdde_interface(this_mat,dee,jump,stress,sdv,dfail)

        type(interface_type),       intent(in)  :: this_mat
        real(kind=dp),              intent(out) :: dee(:,:)
        real(kind=dp),              intent(in)  :: jump(:)
        real(kind=dp),              intent(out) :: stress(:)
        type(sdv_array), optional,  intent(inout) :: sdv
        real(kind=dp),  optional,   intent(in)  :: dfail
        
        ! local variables
        real(kind=dp) :: Dnn0, Dtt0, Dll0
        real(kind=dp) :: sigma(size(dee(:,1)))  ! local replica of stress
        real(kind=dp) :: dmax                   ! local replica of dfail
        integer :: nst, fstat
        
        
        ! initialize variables
        dee=zero; stress=zero ! intent(out) vars
        Dnn0=zero; Dtt0=zero; Dll0=zero
        sigma=zero; dmax=zero
        nst=0; fstat=0
        
        
        
        ! find no. of jumps
        nst=size(dee(:,1))
        
        ! if modulus is not defined, give an error msg and exit
        if(.not.this_mat%modulus_active) then
            write(msg_file,*) 'interface material modulus undefined!'
            call exit_function
        end if
        
        
        ! extract the original linear elasticity stiffness
        Dnn0=this_mat%modulus%Dnn
        Dtt0=this_mat%modulus%Dtt
        Dll0=this_mat%modulus%Dll

        ! form the intact D matrix of the interface
        if (nst .eq. 2) then ! 2D problem
            dee(1,1)=Dnn0
            dee(2,2)=Dtt0
        else if (nst .eq. 3) then ! 3D problem
            dee(1,1)=Dnn0
            dee(2,2)=Dtt0
            dee(3,3)=Dll0
        else
            write(msg_file,*) 'no. of jumps not supported for interface ddsdde!'
            call exit_function
        end if
        
        ! if strength is not present, give a warning and do linear elasticity with stress
        if(.not.this_mat%strength_active) then
            !write(msg_file,*) 'WARNING: strength parameters are not present in interface!'
            !write(msg_file,*) 'Only linear elastic stiffness matrix and stress can be calculated.' 
            ! calculate stress based on intact stiffness
            stress=matmul(dee,jump)
            return ! exit the program
        end if
        
        ! if jump or sdv are not passed in, only linear elasticity can be done
        if(.not.present(sdv)) then
            !write(msg_file,*) 'WARNING: sdv are not present in interface!'
            !write(msg_file,*) 'Only linear elastic stiffness matrix and stress can be calculated.' 
            ! calculate stress based on intact stiffness
            stress=matmul(dee,jump)
            return ! exit the program
        end if    
        
        
        ! --------------------------------------------------------- !
        ! -- reaching here, sdv and strength are present -- !
        ! -- failure criterion can be done                          !
        ! --------------------------------------------------------- !

             
        ! extract max. degradation at total failure
        if(present(dfail))  then
            dmax=dfail ! input parameter, valued at the coh. elem.
        else                
            dmax=one
        end if
      
        
        ! toughness parameters are not active, do only strength failure criterion
        if(.not.this_mat%toughness_active) then
            ! extract current values of failure status variable
            if(.not.allocated(sdv%i)) then 
                allocate(sdv%i(1)); sdv%i=0 ! 1st iteration
            end if
            fstat=sdv%i(1)
            ! check and update fstat
            if(fstat<cohmat_onset) then
                ! calculate stress based on jump
                ! dee is defined based on intact stiffness
                sigma=matmul(dee,jump)               
                ! go through failure criterion and update fstat
                call FailureCriterion(sigma,this_mat%strength,fstat)
            end if          
            ! update D matrix accord. to fstat and crack open/close status
            if(fstat>=cohmat_onset) dee=dee*(one-dmax)
            if(jump(1)<zero) dee(1,1)=Dnn0 ! crack closes, no damage
            ! update stress
            stress=matmul(dee,jump)
            ! update sdv
            sdv%i(1)=fstat          
            ! exit program 
            return
        end if



        ! --------------------------------------------------------- !
        ! reaching here, all variables are present
        ! do cohesive damage evolution to update sdv, Dee and sigma
        ! all inputs are compulsary here, no optional args
        ! --------------------------------------------------------- !


        Call CohesiveLaw(sdv,Dee,sigma,this_mat%strength, &
        & this_mat%toughness,jump,dmax)
        
        ! update stress
        stress=sigma 
        
        ! exit program
        return
        
      end subroutine ddsdde_interface 
      
      
      
      
      
      
      
      subroutine CohesiveLaw(sdv,Dee,sigma,strength,toughness,jump,dmax)
      
        type(sdv_array),        intent(inout)   :: sdv
        real(dp),               intent(inout)   :: Dee(:,:), sigma(:)
        type(interface_strength), intent(in)    :: strength
        type(interface_toughness),intent(in)    :: toughness
        real(dp),                 intent(in)    :: jump(:), dmax
      
        ! local variables
        real(dp) :: Dnn0, dm, u0, uf, Gnc, Gtc, Glc, Gsc, eta, findex, &
        & Gn, Gt, Gl, Gs, bk, Gmc, u_eff, T_eff, T0, dm2
        integer  :: nst, fstat, i, j, l

        real(dp) :: Dres(size(dee(:,1)),size(dee(1,:)))

        ! initialize local variables
        Dnn0=zero; dm=zero; u0=zero; uf=zero; Gnc=zero; Gtc=zero
        Glc=zero;  Gsc=zero;  eta=zero; findex=zero
        Gn=zero;  Gt=zero; Gl=zero; Gs=zero; bk=zero
        Gmc=zero; u_eff=zero; T_eff=zero; T0=zero; dm2=zero
        nst=0; fstat=0; i=0; j=0; l=0

        Dres=zero


        ! extract Dnn0 and nst
        Dnn0=dee(1,1)
        nst=size(dee(:,1))
        
        ! extract current values of failure status variable
        if(.not.allocated(sdv%i)) then
            allocate(sdv%i(1)); sdv%i=0 ! 1st iteration
        end if
        fstat=sdv%i(1)
        
        do i=1, nst
            Dres(i,i)= Kres
        end do


        ! --------------------------------------------------------- !
        ! if already failed, calculate directly Dee and Sigma and exit
        ! --------------------------------------------------------- !
        if(fstat==cohmat_failed) then    ! already failed
            ! calculate penalty stiffness
            dee=dee*(one-dmax)
            dee=max(dee,Dres)
            if(jump(1)<zero) dee(1,1)=Dnn0 ! crack closes, no damage
            ! calculate stress
            sigma=matmul(dee,jump)
            ! exit program
            return
        end if
        
        
        
        ! --------------------------------------------------------- !
        ! the following assumes linear cohesive law
        ! other types of cohesive law can also be used in the future
        ! just need to put in a selection criterion and algorithm
        ! --------------------------------------------------------- ! 
        ! e.g.: select case (toughness%CLtype)
        !           case(0)
        !               linear
        !           case(1)
        !               exponential
        !           ...
        
        
        ! extract damage parameters of last converged iteration
        if(.not.allocated(sdv%r)) then ! 1st iteration
            allocate(sdv%r(3))
            sdv%r=zero
        end if
        dm=sdv%r(1)
        u0=sdv%r(2)
        uf=sdv%r(3)
        
        ! extract toughness parameters
        Gnc=toughness%Gnc
        Gtc=toughness%Gtc
        Glc=toughness%Glc
        eta=toughness%eta 
     
     
        ! check and update fstat and damage variables
        
        ! still intact
        if(fstat==intact) then
            ! calculate stress based on jump
            ! dee is defined based on intact stiffness
            sigma=matmul(dee,jump)
            
            ! go through failure criterion and update fstat
            call FailureCriterion(sigma,strength,fstat,findex)
            
            ! failure onset
            if(fstat==cohmat_onset) then
            ! calculate u0, uf and go to next control
            
                ! normal and shear strain energy density
                Gn=half*max(zero,sigma(1))*max(zero,jump(1))
                if(nst==2)  then
                    Gt=half*sigma(2)*jump(2)
                    Gs=Gt
                else        
                    !Gs=half*(sigma(2)*jump(2)+sigma(3)*jump(3))
                    Gt=half*sigma(2)*jump(2)
                    Gl=half*sigma(3)*jump(3)
                    Gs=Gt+Gl
                end if
                
                ! BK ratio
                bk=Gs/(Gn+Gs)
                
                
                ! effective jump and traction
                if(nst==2)  then
                    u_eff=sqrt(max(zero,jump(1))**2+jump(2)**2)
                else        
                    u_eff=sqrt(max(zero,jump(1))**2+jump(2)**2+jump(3)**2)
                end if
                T_eff=two*(Gn+Gs)/u_eff
                
                ! effective jump and traction at failure onset, taking into acc. of overshoot
                u0=u_eff/sqrt(findex)
                T0=T_eff/sqrt(findex)
                
                ! mixed mode fracture toughness (BK formula)
                !Gsc=sqrt(Gtc**2+Glc**2) ! a quadratic avg of the two shear toughness
                if(Gs>tiny(one)) then
                    Gsc=Gtc*(Gt/Gs)+Glc*(Gl/Gs)
                else
                    Gsc=half*Gtc+half*Glc
                end if
                Gmc=Gnc+(Gsc-Gnc)*(bk**eta)
                
                ! effective jump at final failure
                uf=two*Gmc/T0
                
                
            else if(fstat==cohmat_failed) then
            ! this is the case where strength is close to zero
                dm=dmax
                
            else ! fstat==intact
            ! fstat remains intact; dm, u0 and uf remain zero as initialized
                continue
            end if
        end if
        
        ! failure started
        if(fstat==cohmat_onset) then

            ! calculate dm
            if(uf <= u0 + tiny(one)) then ! brittle failure
                dm=dmax
            else
                ! effective jump and traction
                if(nst==2)  then
                    u_eff=sqrt(max(zero,jump(1))**2+jump(2)**2)
                else        
                    u_eff=sqrt(max(zero,jump(1))**2+jump(2)**2+jump(3)**2)
                end if

                ! linear cohesive softening law 
                if(u_eff > tiny(one)) then
                    dm2=uf/u_eff*(u_eff-u0)/(uf-u0)
                    dm=max(dm,dm2) ! must be larger than last equilibrium dm
                end if
            end if
        
            ! check dm and update fstat
            if (dm>=dmax) then
                dm=dmax
                fstat=cohmat_failed
            end if
        end if
     

        
        ! update D matrix
        dee=dee*(one-dm)
        if(jump(1)<zero) dee(1,1)=Dnn0 ! crack closes, no damage

        dee=max(dee,Dres)
        
        ! update stress
        sigma=matmul(dee,jump)
        
        ! update sdv
        sdv%i(1)=fstat
        sdv%r(1)=dm
        sdv%r(2)=u0
        sdv%r(3)=uf
        
      end subroutine CohesiveLaw
      
      
      subroutine FailureCriterion(sigma,strength,fstat,findex)
      
        real(dp),                   intent(in) :: sigma(:)
        type(interface_strength),   intent(in) :: strength
        integer,                   intent(out) :: fstat
        real(dp), optional,        intent(out) :: findex
      
        ! local variables
        integer     :: nst
        real(dp)    :: tau_nc, tau_tc, tau_lc, findex2
        
        ! initialize variables
        if(present(findex)) findex=zero                     ! intent (out)
        tau_nc=zero; tau_tc=zero; tau_lc=zero; findex2=zero ! local
        nst=0; fstat=0
        
        ! extract nst
        nst=size(sigma)
      
        ! extract strength parameters

        tau_nc=strength%tau_nc
        tau_tc=strength%tau_tc
        tau_lc=strength%tau_lc
        
        ! --------------------------------------------------------- !
        ! the following assumes quadratic stress criterion
        ! other criteria can also be used in the future
        ! just need to put in a selection criterion and algorithm
        ! --------------------------------------------------------- ! 
        ! e.g.: select case (strength%FC)
        !           case(0)
        !               quad. stress
        !           case(1)
        !               max. stress
        !           ...
        
        if(min(tau_nc,tau_tc,tau_lc) > tiny(one)) then
            if(nst==3) then 
            findex2=sqrt((max(sigma(1),zero)/tau_nc)**2 + (sigma(2)/tau_tc)**2 + (sigma(3)/tau_lc)**2)
            else
            findex2=sqrt((max(sigma(1),zero)/tau_nc)**2 + (sigma(2)/tau_tc)**2)
            end if
            
            if(findex2>=one) fstat=cohmat_onset
        else
            ! strength close to zero == already failed
            fstat=cohmat_failed
        end if
        
        if(present(findex)) findex=findex2
      
      end subroutine FailureCriterion
      
      
    end module interface_type_module
