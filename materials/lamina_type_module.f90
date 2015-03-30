module lamina_type_module
!
!  Purpose:
!    to define an object to represent a lamina material
!    
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    30/03/15  B. Y. Chen            Original code
!
use parameter_module, only : DP, ZERO, ONE, TWO, RESIDUAL_MODULUS, MSGLENGTH, &
& STAT_SUCCESS, STAT_FAILURE

implicit none
private

! elastic moduli
type, public :: lamina_modulus
  real(DP) :: E1   = ZERO, E2   = ZERO 
  real(DP) :: G12  = ZERO, G23  = ZERO
  real(DP) :: nu12 = ZERO, nu23 = ZERO
end type

! strength properties (Hashin or Pinho criterion)
type, public :: lamina_strength
  real(DP) :: Xt = ZERO, Xc = ZERO ! fibre tens. & comp. strengths
  real(DP) :: Yt = ZERO, Yc = ZERO ! matrix tens. & comp. strengths
  real(DP) :: Sl = ZERO, St = ZERO ! matrix long. & trans. shear strengths
end type

! matrix toughness properties
type, public :: lamina_matrixToughness
  real(DP) :: GmcI = ZERO, GmcII = ZERO  ! matrix toughness, mode I and II 
  real(DP) :: eta  = ZERO                ! mixed-mode law constant (BK)
end type

! fibre toughness properties, tensile and compressive
type, public :: lamina_fibreToughness
  real(DP) :: GfcT = ZERO, GfcC = ZERO 
end type

! lamina material object definition
type, public :: lamina_type
  type(lamina_modulus)          :: modulus
  type(lamina_strength)         :: strength
  type(lamina_matrixToughness)  :: matrixToughness
  type(lamina_fibreToughness)   :: fibreToughness
end type


interface empty
    module procedure empty_lamina
end interface
    
interface update
    module procedure update_lamina
end interface

interface extract
    module procedure extract_lamina
end interface

interface ddsdde
    module procedure ddsdde_lamina
end interface


public :: empty,update,extract,ddsdde




contains




  pure subroutine empty_lamina(this)
  ! Purpose:
  ! to reset this lamina object components into their default values (ZERO)
 
    type(lamina_type),intent(inout) :: this
    
    this%modulus         = lamina_modulus(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO)
    this%strength        = lamina_strength(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO)
    this%matrixToughness = lamina_matrixtoughness(ZERO,ZERO,ZERO)
    this%fibreToughness  = lamina_fibretoughness(ZERO,ZERO)

  end subroutine empty_lamina



  pure subroutine update_lamina(this, modulus, strength, fibreToughness,&
  & matrixToughness)
  ! Purpose:
  ! to update this lamina object components
    
   	type(lamina_type),                  intent(inout) :: this
    type(lamina_modulus),         optional,intent(in) :: modulus
    type(lamina_strength),        optional,intent(in) :: strength
    type(lamina_matrixToughness), optional,intent(in) :: matrixToughness
    type(lamina_fibreToughness),  optional,intent(in) :: fibreToughness
            
    if(present(modulus))          this%modulus          = modulus    
    if(present(strength))         this%strength         = strength    
    if(present(matrixToughness))  this%matrixToughness  = matrixToughness    
    if(present(fibreToughness))   this%fibreToughness   = fibreToughness

  end subroutine update_lamina  
    
    
    
  pure subroutine extract_lamina(this, modulus, strength, matrixToughness,&
  & fibreToughness)
  ! Purpose:
  ! to extract componets of this lamina object
    
   	type(lamina_type),                      intent(in) :: this
    type(lamina_modulus),         optional,intent(out) :: modulus
    type(lamina_strength),        optional,intent(out) :: strength
    type(lamina_matrixToughness), optional,intent(out) :: matrixToughness
    type(lamina_fibreToughness),  optional,intent(out) :: fibreToughness
    
    if(present(modulus))          modulus         = this%modulus    
    if(present(strength))         strength        = this%strength
    if(present(matrixToughness))  matrixtoughness = this%matrixToughness
    if(present(fibreToughness))   fibretoughness  = this%fibreToughness

  end subroutine extract_lamina 




  pure subroutine ddsdde_lamina(dee, stress, sdv, this_mat, strain, clength, &
  & maxdm)
  ! Purpose:
  ! to calculate D matrix, stress and solution-dependent variables
  ! at an integration point

    ! dummy argument list:
    ! - dee         : local stiffness matrix D, to be updated
    ! - stress      : local stress vector, to be updated
    ! - sdv         : solution-dependent variables array, to be updated
    ! - this_mat    : material properties object, passed-in
    ! - strain      : local strain vector, passed-in
    ! - clength     : elem. characteristic length, passed-in
    ! - maxdm       : maximum degradation factor
    real(DP),          intent(inout) :: dee(:,:)
    real(DP),          intent(inout) :: stress(:)
    type(sdv_array),   intent(inout) :: sdv
    type(lamina_type), intent(in)    :: this_mat
    real(DP),          intent(in)    :: strain(:)    
    real(DP),          intent(in)    :: clength
    real(DP),optional, intent(in)    :: maxdm
    
    ! local variables list:
    ! - df      : fibre degradation
    ! - cmaxdm  : local copy of maxdm
    ! - nStrain : no. strains
    ! - fstat   : failure status
    ! - ffstat  : fibre failure status
    ! - mfstat  : matrix failure status
    real(DP) :: df, cmaxdm
    integer  :: nStrain, fstat, ffstat, mfstat

    
    ! initialize local variables
    df      = ZERO 
    cmaxdm  = ZERO
    nStrain = 0 
    fstat   = 0 
    ffstat  = 0 
    mfstat  = 0
        
    ! calculate the linear elasticity stiffness matrix   
    call deemat(dee, this_mat)
            
    ! copy max. degradation (if present) into local variable
    if(present(maxdm))  then
        cmaxdm=maxdm
    else  
        ! default value is ONE              
        cmaxdm=ONE
    end if
    
    ! extract current values of failure status variable
    if(.not.allocated(sdv%i)) then 
        allocate(sdv%i(3)); sdv%i=0 ! 1st iteration
    end if
    fstat=sdv%i(1)  ! generic failure status
    ffstat=sdv%i(2) ! fibre failure status
    mfstat=sdv%i(3) ! matrix failure status
    
    
    if(ffstat==fibre_failed) then
    ! already completely failed, no need to go through coh law; update dee and calculate stress
        df=sdv%r(1)
        call deemat(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,dee,df,df,df)
        sig=matmul(dee,strain)
    else
        ! when fibres are intact, calculate stress for failure criterion check
        if(ffstat==intact) sig=matmul(dee,strain)
        ! use cohesive law, calculate damage var. (sdv%r) with strain (and sig when intact)
        call FibreCohesiveLaw(ffstat,sdv%r,sig,this_mat%strength, &
&       this_mat%fibretoughness,strain,clength,cmaxdm)
        ! update sdv
        sdv%i(2)=ffstat
        ! update ddsdde
        df=sdv%r(1) ! fibre stiffness degradation
        call deemat(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,dee,df,df,df)
        ! update stress
        sig=matmul(dee,strain)
    end if

    ! do only strength failure criterion
    if(this_mat%matrixtoughness_active) then
        write(msg_file,*)'matrix failure smeared crack model not supported! &
        & only failure criterion is assessed; use FNM cohesive subelem instead for matrix cracking'
    end if


    ! check and update mfstat
    if(mfstat<matrix_onset) then               
        ! go through failure criterion and update fstat
        call MatrixFailureCriterion(sig,this_mat%strength,mfstat)
    end if
    ! update sdv
    if(mfstat>=matrix_onset) sdv%i(3)=mfstat
    
    !~! degrade matrix stiffness if fibres have all failed
    !~if(ffstat==fibre_failed .and. mfstat>=matrix_onset) then
    !~    call deemat(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,dee,df,df,df)
    !~    ! update stress
    !~    sig=matmul(dee,strain)
    !~end if
    
    ! update stress
    if(present(stress)) stress=sig
    
    ! update fstat
    fstat=max(ffstat,mfstat)
    sdv%i(1)=fstat
    
    ! exit program
    return
  end subroutine ddsdde_lamina 




  pure subroutine deemat(dee,this_mat,df,dm2,dm3,stat,msg)
  ! Purpose:
  ! to calculate local stiffness matrix D
    
    ! dummy argument list:
    ! - dee         : local stiffness matrix
    ! - this_mat    : passed-in material object
    ! - df          : fibre degradation
    ! - dm2         : matrix degradation (in-plane, dir. 2)
    ! - dm3         : matrix degradation (out-plane, dir. 3)
    real(DP),           intent(inout) :: dee(:,:)
    type(lamina_type),  intent(in)    :: this_mat
    real(DP),                 optional, intent(in)  :: df, dm2, dm3
    integer,                  optional, intent(out) :: stat
    character(len=MSGLENGTH), optional, intent(out) :: msg
    
    !local var. list:
    ! - E1,...,G23 : elastic moduli, standard notations
    ! - del        : a common denominator for calculating dee
    ! - dm23       : matrix degradation (transverse, 2-3 plane)
    ! - nStrain    : no. of stresses/strains
    ! - istat      : local copy of stat, i.e., status variable
    ! - emsg       : local copy of msg, i.e., error message
    real(DP) :: E1,E2,E3,nu12,nu13,nu23,nu21,nu31,nu32,G12,G13,G23
    real(DP) :: del, dm23
    integer  :: nStrain
    integer  :: istat
    character(len=MSGLENGTH) :: emsg

    
    ! initialize local variables
    E1   = ZERO;  E2   = ZERO;  E3   = ZERO
    G12  = ZERO;  G13  = ZERO;  G23  = ZERO
    nu12 = ZERO;  nu13 = ZERO;  nu23 = ZERO
    nu21 = ZERO;  nu31 = ZERO;  nu32 = ZERO
    del  = ZERO;  dm23 = ZERO;  nStrain = 0
    istat = 0;    emsg = ''
    
    ! extract elastic muduli from passed-in material object
    E1   = this_mat%modulus%E1
    E2   = this_mat%modulus%E2 
    E3   = E2
    G12  = this_mat%modulus%G12
    G13  = G12
    G23  = this_mat%modulus%G23
    nu12 = this_mat%modulus%nu12
    nu13 = nu12
    nu23 = this_mat%modulus%nu23
    nu21 = (E2 / E1) * nu12
    nu31 = nu21
    nu32 = (E3 / E2) * nu23
    

    ! apply fibre degradation if present        
    if(present(df)) then    
      E1   = E1   * (ONE - df)
      nu12 = nu12 * (ONE - df)
      nu13 = nu13 * (ONE - df)      
      ! do not degrade below residual stiffness
      if (E1 < RESIDUAL_STIFFNESS) then
        E1   = RESIDUAL_STIFFNESS
        nu12 = ZERO
        nu13 = ZERO
      end if      
    end if

    
    ! apply in-plane matrix degradation if present
    if(present(dm2)) then    
      E2   = E2   * (ONE - dm2)
      nu21 = nu21 * (ONE - dm2)
      nu23 = nu23 * (ONE - dm2)
      G12  = G12  * (ONE - dm2)      
      ! do not degrade below residual stiffness
      if (E2 < RESIDUAL_STIFFNESS) then
        E2   = RESIDUAL_STIFFNESS
        nu21 = ZERO
        nu23 = ZERO
      end if      
      G12  = max(G12, RESIDUAL_STIFFNESS)      
      ! update matrix degradation for 2-3 plane
      dm23 = max(dm23, dm2)      
    end if

    
    ! apply out-plane matrix degradation if present
    if(present(dm3)) then    
      E3   = E3   * (ONE - dm3)
      nu31 = nu31 * (ONE - dm3)
      nu32 = nu32 * (ONE - dm3)
      G13  = G13  * (ONE - dm3)      
      ! do not degrade below residual stiffness
      if (E3 < RESIDUAL_STIFFNESS) then
        E3   = RESIDUAL_STIFFNESS
        nu31 = ZERO
        nu32 = ZERO
      end if      
      G13  = max(G13, RESIDUAL_STIFFNESS)      
      ! update matrix degradation for 2-3 plane
      dm23 = max(dm23, dm3)      
    end if

    
    ! apply transverse (2-3) matrix degradation 
    G23 = G23*(ONE-dm23)
    ! do not degrade below residual stiffness
    G23 = max(G23,RESIDUAL_STIFFNESS)
        
    
    ! find no. of strains
    nStrain=size(dee(:,1))
    
    ! calculate D matrix based on the no. of strains
    ! only 3D problem with the standard 6 strains is supported
    select case (nStrain)
      ! 3D problem
      case (6)  
        ! calculate D matrix terms
        dee = ZERO  ! zero all terms first
        del = ONE-nu12*nu21-nu13*nu31-nu23*nu32-TWO*nu21*nu32*nu13
        del = del/E1/E2/E3
        dee(1,1)= (ONE-nu23*nu32)/E2/E3/del
        dee(1,2)= (nu21+nu23*nu31)/E2/E3/del
        dee(1,3)= (nu31+nu21*nu32)/E2/E3/del
        dee(2,1)= dee(1,2)
        dee(2,2)= (ONE-nu13*nu31)/E1/E3/del
        dee(2,3)= (nu32+nu12*nu31)/E1/E3/del
        dee(3,1)= dee(1,3)
        dee(3,2)= dee(2,3)
        dee(3,3)= (ONE-nu12*nu21)/E1/E2/del
        dee(4,4)= G12
        dee(5,5)= G13
        dee(6,6)= G23
        ! update local status variable; no need to update error message
        istat = STAT_SUCCESS
      ! unsupported no. of strains
      case default  
        ! update local status variable and error message
        istat = STAT_FAILURE
        emsg  = 'no. strains is not supported in deemat, lamina_type_module'
    end select
    
    ! return status variable and error message, if required
    if(present(stat)) stat = istat
    if(present(msg))  msg  = emsg
    
  end subroutine deemat




  subroutine FibreCohesiveLaw(fstat,sdvr,sigma,strength,fibretoughness,strain,clength,cmaxdm)
    
    integer,                intent(inout)   :: fstat
    real(DP), allocatable,  intent(inout)   :: sdvr(:)
    real(DP),                 intent(in)    :: sigma(:),strain(:), clength, cmaxdm 
    type(lamina_strength), intent(in)       :: strength
    type(lamina_fibretoughness),intent(in)  :: fibretoughness
    
    
    ! local variables
    real(DP) :: GfcT, GfcC, dm, u0, uf, T0, u_eff, T_eff, findex, dm2
    integer  :: nStrain
    
    
    ! --------------------------------------------------------- !
    ! if already failed, calculate directly Dee and Sigma and exit
    ! --------------------------------------------------------- !
    if(fstat==fibre_failed) then    ! already failed
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
    
    ! initialize local variables
    GfcT=ZERO; GfcC=ZERO
    dm=ZERO; u0=ZERO; uf=ZERO; T0=ZERO; findex=ZERO; dm2=ZERO
    u_eff=ZERO; T_eff=ZERO
    nStrain=0
    
    
    ! extract toughness parameters
    GfcT=fibretoughness%GfcT
    GfcC=fibretoughness%GfcC

    ! extract current values of failure status variable
    if(.not.allocated(sdvr)) then
        allocate(sdvr(3)); sdvr=ZERO  ! 1st iteration
    end if
    dm=sdvr(1)
    u0=sdvr(2)
    uf=sdvr(3) 
    
    ! extract nStrain
    nStrain=size(strain)
    
    
    ! check and update fstat and damage variables
    
    ! still intact
    if(fstat==intact) then
        
        ! go through failure criterion and update fstat
        call FibreFailureCriterion(sigma,strength,fstat,findex)
        
        ! failure onset
        if(fstat==fibre_onset) then
        ! calculate u0, uf and go to next control
            
            u_eff=abs(strain(1))*clength
            T_eff=abs(sigma(1))
            
            ! effective jump and traction at failure onset, taking into acc. of overshoot
            u0=u_eff/sqrt(findex)
            T0=T_eff/sqrt(findex)
            
            ! effective jump at final failure
            if(strain(1)>ZERO) then
                uf=TWO*GfcT/T0
            else
                uf=TWO*GfcC/T0
            end if
            
            ! calculate dm to prevent overshoot of stress
            if(uf <= u0 + tiny(ONE)) then ! brittle failure
                dm=cmaxdm
                fstat=fibre_failed
            else
                dm=ONE-ONE/sqrt(findex)
            end if
            
        else if(fstat==fibre_failed) then
        ! this is the case where strength is close to ZERO
            dm=cmaxdm
            
        else ! fstat==intact
        ! fstat remains intact; dm, u0 and uf remain ZERO as initialized
            continue
        end if
    !end if
    
    ! failure started
    else if(fstat==fibre_onset) then
    
        ! calculate dm
        if(uf <= u0 + tiny(ONE)) then ! brittle failure
            dm=cmaxdm
            fstat=fibre_failed
        else
            ! effective jump and traction
            u_eff=abs(strain(1))*clength
            ! linear cohesive softening law 
            if(u_eff > tiny(ONE)) then
                dm2=uf/u_eff*(u_eff-u0)/(uf-u0)
                dm=max(dm,dm2) ! must be larger than last equilibrium dm
            end if
        end if
    
        ! check dm and update fstat
        if (dm>=cmaxdm-tiny(ONE)) then
            dm=cmaxdm
            fstat=fibre_failed
        end if
    end if
    
    
    
    ! update sdv
    sdvr(1)=dm
    sdvr(2)=u0
    sdvr(3)=uf
    
  end subroutine FibreCohesiveLaw





  subroutine FibreFailureCriterion(sigma,strength,fstat,findex)
    ! at the moment, only tensile failure is considered
    
    real(DP),                   intent(in) :: sigma(:)
    type(lamina_strength),      intent(in) :: strength
    integer,                   intent(out) :: fstat
    real(DP), optional,        intent(out) :: findex
    
    
    
    ! local variables
    real(DP)    :: Xt, Xc ! tensile, compressive strengths
    real(DP)    :: findex2

    
    ! initialize variables
    if(present(findex)) findex=ZERO                     ! intent (out)
    Xt=ZERO; Xc=ZERO; findex2=ZERO
    fstat=0
    
    
    ! extract strength parameters
    Xt=strength%Xt
    Xc=strength%Xc
    
    
    ! --------------------------------------------------------- !
    ! the following assumes max stress criterion
    ! other criteria can also be used in the future
    ! just need to put in a selection criterion and algorithm
    ! --------------------------------------------------------- ! 
    ! e.g.: select case (strength%FC)
    !           case(0)
    !               quad. stress
    !           case(1)
    !               max. stress
    !           ...
    
    if(min(Xt,abs(Xc)) > tiny(ONE)) then
        ! failure index for tensile failure; matrix crack perpendicular to shell plane, no z-dir stress components
        findex2=max(sigma(1),ZERO)/Xt + abs(min(sigma(1),ZERO)/Xc)
        
        if(findex2>=ONE) fstat=fibre_onset
    else
        ! strength close to ZERO == already failed
        fstat=fibre_failed
    end if
    
    if(present(findex)) findex=findex2
    
  end subroutine FibreFailureCriterion
    



    
  subroutine MatrixFailureCriterion(sigma,strength,fstat,findex)
    ! at the moment, only tensile failure is considered
    
    real(DP),                   intent(in) :: sigma(:)
    type(lamina_strength),      intent(in) :: strength
    integer,                   intent(out) :: fstat
    real(DP), optional,        intent(out) :: findex
    
    
    
    ! local variables
    integer     :: nStrain
    real(DP)    :: Yt,Yc,Sl,St ! tensile, compressive, and shear strengths
    real(DP)    :: findex2

    
    ! initialize variables
    if(present(findex)) findex=ZERO                     ! intent (out)
    Yt=ZERO; Yc=ZERO; Sl=ZERO; St=ZERO; findex2=ZERO
    nStrain=0; fstat=0
    
    ! extract nStrain
    nStrain=size(sigma)
    
    ! extract strength parameters
    Yt=strength%Yt
    Yc=strength%Yc
    Sl=strength%Sl
    St=strength%St
    
    
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
    
    if(min(Yt,St,Sl) > tiny(ONE)) then
        ! failure index for tensile failure; matrix crack perpendicular to shell plane, no z-dir stress components
        if(nStrain==3) then 
        findex2=sqrt((max(sigma(2),ZERO)/Yt)**2 + (sigma(3)/Sl)**2)
        else
        findex2=sqrt((max(sigma(2),ZERO)/Yt)**2 + (sigma(4)/Sl)**2)
        end if
        
        if(findex2>=ONE) fstat=matrix_onset
    else
        ! strength close to ZERO == already failed
        fstat=matrix_onset
    end if
    
    if(present(findex)) findex=findex2
    
  end subroutine MatrixFailureCriterion
    
    
    
end module lamina_type_module
