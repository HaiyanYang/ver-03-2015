module cohesive_material_module
!
!  Purpose:
!    define an object to represent a cohesive interface material
!    with the associated procedures to empty, update, extract and display
!    its components, to integrate local stiffness matrix, and
!    to update local solution-dependent variables (damage variables)
!    
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    05/04/15  B. Y. Chen            Original code
!

! list out ALL the parameter used, no matter how long:
! DP                : precision of real number, Double Precision
! ZERO, ONE, TWO    : self explanatory
! SMALLNUM          : a very small real number, used when comparing two reals
! RESIDUAL_MODULUS  : residual modulus of the material at final failure
! MSGLENGTH         : length of error message
! STAT_SUCCESS      : integer value of a successful status
! STAT_FAILURE      : integer value of a failure    status
! INTACT            : integer value for generic intact state
! COH_MAT_ONSET     : integer value for cohesive material failure onset state
! COH_MAT_FAILED    : integer value for cohesive material total failure state
use parameter_module, only : DP, ZERO, ONE, TWO, SMALLNUM, RESIDUAL_MODULUS, &
& MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, INTACT, COH_MAT_ONSET, COH_MAT_FAILED


implicit none
private

! penalty stiffness, normal and two shear directions
type,public :: cohesive_modulus
    real(DP) :: Dnn = ZERO
    real(DP) :: Dtt = ZERO
    real(DP) :: Dll = ZERO
end type

! strengths for the three directions
type,public :: cohesive_strength
    real(DP) :: tau_nc = ZERO
    real(DP) :: tau_tc = ZERO
    real(DP) :: tau_lc = ZERO
end type

! fracture toughness for the three directions
! Gnc, Gtc and Glc : toughness values, t and l dir. values may be different
! eta              : mixed-mode law component; eta_nt/nl/tl are the same; 
!                    may be different in future
type,public :: cohesive_toughness
    real(DP) :: Gnc = ZERO
    real(DP) :: Gtc = ZERO
    real(DP) :: Glc = ZERO
    real(DP) :: eta = ZERO
end type


type,public :: cohesive_material
    type(cohesive_modulus)   :: modulus
    type(cohesive_strength)  :: strength
    type(cohesive_toughness) :: toughness
end type

! cohesive material solution-dependent variables
! dm     : cohesive modulus degradation factor
! u0, uf : cohesive law parameters, initial & final failure displacements
! fstat  : failure status
type, public :: cohesive_sdv
  real(DP) :: dm    = ZERO,  u0 = ZERO,  uf = ZERO
  integer  :: fstat = INTACT
end type


interface empty
    module procedure empty_cohesive
end interface
  
interface update
    module procedure update_cohesive
end interface

interface display
  module procedure display_cohesive
end interface

interface display
    module procedure display_cohesive_sdv
end interface

interface ddsdde
    module procedure ddsdde_cohesive
end interface


public :: empty, update, display, ddsdde




contains



  pure subroutine empty_cohesive (this)
  ! Purpose:
  ! to reset this object's components into their default values (ZERO)
  
    type(cohesive_material), intent(inout) :: this
    
    this%modulus   = cohesive_modulus   (ZERO, ZERO, ZERO)
    this%strength  = cohesive_strength  (ZERO, ZERO, ZERO)
    this%toughness = cohesive_toughness (ZERO, ZERO, ZERO, ZERO)
    
  end subroutine empty_cohesive



  pure subroutine update_cohesive (this, modulus, strength, toughness)
  ! Purpose :
  ! to update this object's components
  
    type(cohesive_material),            intent(inout) :: this
    type(cohesive_modulus),   optional, intent(in)    :: modulus
    type(cohesive_strength),  optional, intent(in)    :: strength
    type(cohesive_toughness), optional, intent(in)    :: toughness
    
    if (present(modulus))   this%modulus   = modulus
    if (present(strength))  this%strength  = strength
    if (present(toughness)) this%toughness = toughness

  end subroutine update_cohesive   
  
  
  
  subroutine display_cohesive(this)
  ! Purpose:
  ! to display this cohesive object's components on cmd window
  ! this is useful for debugging
 
    type(cohesive_material), intent(in) :: this
    
    ! local variable
    character(len=20) :: display_fmt
    
    ! initialize local variable
    display_fmt = ''

    ! set display format for string and real
    ! A for string, ES for real (scientific notation)
    ! 10 is width, 3 is no. of digits aft decimal point
    ! note that for scientific real, ESw.d, w>=d+7
    display_fmt = '(1X, A, ES10.3)' 
    
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the modulus of the inquired cohesive object :'
    write(*,display_fmt) 'cohesive Dnn    is: ', this%modulus%Dnn 
    write(*,display_fmt) 'cohesive Dtt    is: ', this%modulus%Dtt
    write(*,display_fmt) 'cohesive Dll    is: ', this%modulus%Dll
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the strength of the inquired cohesive object :'
    write(*,display_fmt) 'cohesive tau_nc is: ', this%strength%tau_nc 
    write(*,display_fmt) 'cohesive tau_tc is: ', this%strength%tau_tc
    write(*,display_fmt) 'cohesive tau_lc is: ', this%strength%tau_lc
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the toughness of the inquired cohesive object :'
    write(*,display_fmt) 'cohesive Gnc    is: ', this%toughness%Gnc 
    write(*,display_fmt) 'cohesive Gtc    is: ', this%toughness%Gtc 
    write(*,display_fmt) 'cohesive Glc    is: ', this%toughness%Glc
    write(*,display_fmt) 'cohesive eta    is: ', this%toughness%eta
    write(*,'(1X, A)') ''

  end subroutine display_cohesive
  
  
  
  subroutine display_cohesive_sdv(this_sdv)
  ! Purpose:
  ! to display this cohesive_sdv's components on cmd window
  ! this is useful for debugging
  
    type(cohesive_sdv), intent(in) :: this_sdv
  
    ! local variable
    character(len=20) :: display_fmt
    
    ! initialize local variable
    display_fmt = ''

    ! set display format for string and integer
    ! A for string, I for integer, 10 for width of the number
    display_fmt = '(1X, A, I10)' 
    
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the inquired cohesive SDVs :'
    write(*,display_fmt) 'cohesive FSTAT is: ', this_sdv%FSTAT 

    ! set display format for string and real
    ! A for string, ES for real (scientific notation)
    ! 10 is width, 3 is no. of digits aft decimal point
    ! note that for scientific real, ESw.d, w>=d+7
    display_fmt = '(1X, A, ES10.3)' 
    
    write(*,display_fmt) 'cohesive dm    is: ', this_sdv%dm 
    write(*,display_fmt) 'cohesive U0    is: ', this_sdv%U0
    write(*,display_fmt) 'cohesive UF    is: ', this_sdv%UF
    write(*,'(1X, A)') ''
  
  end subroutine display_cohesive_sdv
  


  subroutine ddsdde_cohesive (this_mat, dee, stress, sdv, jump, &
  & istat, emsg, d_max)
  ! Purpose:
  ! to calculate the D matrix, stress and solution-dependent variables
  ! at an integration point of an element of lamina_material definition
  ! (restricted to 3D problems, with the standard 6 strain terms)

    ! dummy argument list:
    ! - this_mat    : material properties object,         passed-in (pass arg)
    ! - dee         : local stiffness matrix D,           to update
    ! - stress      : local stress vector,                to update
    ! - sdv         : solution-dependent variables array, to update
    ! - jump        : local displacement jump vector,     passed-in
    ! - istat       : status variable of this procedure   to output
    ! - emsg        : error message                       to output
    ! - d_max       : maximum degradation factor          passed-in (optional)
  
    type(cohesive_material),  intent(in)    :: this_mat
    real(DP),                 intent(inout) :: dee(:,:)
    real(DP),                 intent(inout) :: stress(:)
    type(cohesive_sdv),       intent(inout) :: sdv
    real(DP),                 intent(in)    :: jump(:)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    real(DP),       optional, intent(in)    :: d_max
    
    ! local variables list:
    ! - fstat     : failure status
    ! - dm        : degradation factor
    ! - u0        : cohesive law, failure onset displacement
    ! - uf        : cohesive law, total failure displacement
    ! - d_max_lcl : local copy of d_max
    integer  :: fstat
    real(DP) :: dm, u0, uf
    real(DP) :: d_max_lcl
    
    
    ! initialize intent(out) & local variables
    istat     = STAT_SUCCESS  ! default
    emsg      = ''
    fstat     = 0
    dm        = ZERO
    u0        = ZERO
    uf        = ZERO
    d_max_lcl = ZERO
    
    ! check validity of dummy arguments with intent(in) or (inout)
    
    ! check dee, stress and strain size
    if (.not. (size(dee(:,1)) == 3 .and. size(dee(1,:)) == 3)) then
      istat = STAT_FAILURE
      emsg  = 'size of dee is not supported for ddsdde_cohesive, &
      &cohesive_material_module!'
      return
    end if
    
    if (.not. (size(stress) == 3)) then
      istat = STAT_FAILURE
      emsg  = 'size of stress is not supported for ddsdde_cohesive, &
      &cohesive_material_module!'
      return
    end if
    
    if (.not. (size(jump) == 3)) then
      istat = STAT_FAILURE
      emsg  = 'size of jump is not supported for ddsdde_cohesive, &
      &cohesive_material_module!'
      return
    end if
    
    ! check sdv values
    call check_sdv (sdv, istat, emsg)
    if (istat == STAT_FAILURE) return
    
    ! check this_mat properties
    call check_mat_prop (this_mat, istat, emsg)
    if (istat == STAT_FAILURE) return
    
    ! check d_max value, must be between ZERO and ONE
    if (present(d_max)) then
      if (.not. (ZERO-SMALLNUM < d_max .and. d_max < ONE+SMALLNUM)) then
        istat = STAT_FAILURE
        emsg  = 'd_max must be within [0., 1.] for ddsdde_cohesive, &
        &cohesive_material_module!'
        return
      end if
    end if 
       
         
    ! extract max. degradation at total failure
    if(present(d_max))  then
        d_max_lcl=d_max ! input parameter, valued at the coh. elem.
    else                
        d_max_lcl=one
    end if
  


    ! --------------------------------------------------------- !
    ! reaching here, all variables are present
    ! do cohesive damage evolution to update sdv, Dee and sigma
    ! all inputs are compulsary here, no optional args
    ! --------------------------------------------------------- !


    call cohesive_law (sdv,Dee,sigma,this_mat%strength, &
    & this_mat%toughness,jump,d_max_lcl)
    
    ! update stress
    stress=sigma 
    
    ! exit program
    return
    
  end subroutine ddsdde_cohesive 
  
  
  
  
  
  
  
  subroutine cohesive_law (sdv,Dee,sigma,strength,toughness,jump,d_max)
  
    type(sdv_array),        intent(inout)   :: sdv
    real(dp),               intent(inout)   :: Dee(:,:), sigma(:)
    type(cohesive_strength), intent(in)    :: strength
    type(cohesive_toughness),intent(in)    :: toughness
    real(dp),                 intent(in)    :: jump(:), d_max
  
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
        dee=dee*(one-d_max)
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
            dm=d_max
            
        else ! fstat==intact
        ! fstat remains intact; dm, u0 and uf remain zero as initialized
            continue
        end if
    end if
    
    ! failure started
    if(fstat==cohmat_onset) then

        ! calculate dm
        if(uf <= u0 + tiny(one)) then ! brittle failure
            dm=d_max
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
        if (dm>=d_max) then
            dm=d_max
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
    
  end subroutine cohesive_law
  
  
  subroutine FailureCriterion(sigma,strength,fstat,findex)
  
    real(dp),                   intent(in) :: sigma(:)
    type(cohesive_strength),   intent(in) :: strength
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
  
  
end module cohesive_material_module
