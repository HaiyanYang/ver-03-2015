module lamina_type_module
!
!  Purpose:
!    define an object to represent a lamina material
!    with the associated procedures to empty, update, extract
!    its components, to integrate local stiffness matrix, and
!    to update local solution-dependent variables (damage variables)
!    
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    30/03/15  B. Y. Chen            Original code
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
! MATRIX_ONSET      : integer value for matrix failure onset state
! MATRIX_FAILED     : integer value for matrix total failure state
! FIBRE_ONSET       : integer value for fibre  failure onset state
! FIBRE_FAILED      : integer value for fibre  total failure state
! LAMINA_SDV        : a der. type to group specific variables of lamina material
use parameter_module, only : DP, ZERO, ONE, TWO, SMALLNUM, RESIDUAL_MODULUS, &
& MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, INTACT, MATRIX_ONSET, MATRIX_FAILED,&
& FIBRE_ONSET, FIBRE_FAILED, LAMINA_SDV 

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
    this%matrixToughness = lamina_matrixToughness(ZERO,ZERO,ZERO)
    this%fibreToughness  = lamina_fibreToughness(ZERO,ZERO)

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
    if(present(matrixToughness))  matrixToughness = this%matrixToughness
    if(present(fibreToughness))   fibreToughness  = this%fibreToughness

  end subroutine extract_lamina 
  
  

  pure subroutine ddsdde_lamina(dee, stress, sdv, this_mat, strain, clength, &
  & istat, emsg, maxdm)
  ! Purpose:
  ! to calculate the D matrix, stress and solution-dependent variables
  ! at an integration point of an element with lamina_type material definition.
  ! (restricted to 3D problems, with the standard 6 strain terms)

    ! dummy argument list:
    ! - dee         : local stiffness matrix D,           to update
    ! - stress      : local stress vector,                to update
    ! - sdv         : solution-dependent variables array, to update
    ! - this_mat    : material properties object,         passed-in
    ! - strain      : local strain vector,                passed-in
    ! - clength     : elem. characteristic length,        passed-in
    ! - istat       : status variable of this procedure   to output
    ! - emsg        : error message                       to output
    ! - maxdm       : maximum degradation factor          passed-in (optional)
    real(DP),                 intent(inout) :: dee(:,:)
    real(DP),                 intent(inout) :: stress(:)
    type(LAMINA_SDV),         intent(inout) :: sdv
    type(lamina_type),        intent(in)    :: this_mat
    real(DP),                 intent(in)    :: strain(:)    
    real(DP),                 intent(in)    :: clength
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    real(DP),       optional, intent(in)    :: maxdm
    
    ! local variables list:
    ! - fstat     : generic failure status
    ! - ffstat    : fibre   failure status
    ! - mfstat    : matrix  failure status
    ! - df        : fibre degradation
    ! - u0        : fibre cohesive law, failure onset displacement
    ! - uf        : fibre cohesive law, total failure displacement
    ! - maxdm_lcl : local copy of maxdm
    integer  :: fstat, ffstat, mfstat
    real(DP) :: df, u0, uf
    real(DP) :: maxdm_lcl

    
    ! initialize intent(out) & local variables
    istat     = STAT_SUCCESS  ! default
    emsg      = ''
    fstat     = 0 
    ffstat    = 0 
    mfstat    = 0
    df        = ZERO 
    u0        = ZERO 
    uf        = ZERO 
    maxdm_lcl = ZERO
    
    ! check validity of dummy arguments with intent(in) or (inout)
    
    ! check dee, stress and strain size
    if (.not. (size(dee(:,1)) == 6 .and. size(dee(1,:)) == 6)) then
      istat = STAT_FAILURE
      emsg  = 'size of dee is not supported for ddsdde_lamina, &
      &lamina_type_module!'
      return
    end if
    
    if (.not. (size(stress) == 6)) then
      istat = STAT_FAILURE
      emsg  = 'size of stress is not supported for ddsdde_lamina, &
      &lamina_type_module!'
      return
    end if
    
    if (.not. (size(strain) == 6)) then
      istat = STAT_FAILURE
      emsg  = 'size of strain is not supported for ddsdde_lamina, &
      &lamina_type_module!'
      return
    end if
    
    ! check sdv values
    call check_sdv (sdv, istat, emsg)
    if (istat == STAT_FAILURE) return
    
    ! check this_mat properties
    call check_mat_prop (this_mat, istat, emsg)
    if (istat == STAT_FAILURE) return
    
    ! check clength value
    if (.not. (clength > SMALLNUM)) then 
      istat = STAT_FAILURE
      emsg  = 'clength must be larger than zero for ddsdde_lamina, &
      &lamina_type_module!'
      return
    end if
    
    ! check maxdm value, must be between ZERO and ONE
    if (present(maxdm)) then
      if (.not. (ZERO-SMALLNUM < maxdm .and. maxdm < ONE+SMALLNUM)) then
        istat = STAT_FAILURE
        emsg  = 'maxdm must be within [0., 1.] for ddsdde_lamina, &
        &lamina_type_module!'
        return
      end if
    end if 
        

    ! extract values of failure status and cohesive law variables
    fstat  = sdv%FSTAT 
    ffstat = sdv%FFSTAT
    mfstat = sdv%MFSTAT
    df     = sdv%DF
    u0     = sdv%U0
    uf     = sdv%UF
         
    ! copy max. degradation (if present) into local variable
    if(present(maxdm))  then
      maxdm_lcl = maxdm
    else  
      ! default value is ONE              
      maxdm_lcl = ONE
    end if
    
    !---------------------------------------------------------------------------
    ! fstat is always equal to ffstat, 
    ! unless ffstat = INTACT and mfstat > INTACT, then fstat = mfstat
    !
    ! Possible changes of ffstat, mfstat and fstat after the failure procedure
    !
    ! From:  ffstat = FIBRE_FAILED; mfstat = * ; fstat = ffstat
    ! To  :  ffstat = FIBRE_FAILED; mfstat = * ; fstat = ffstat (NO CHANGE)
    !
    ! From:  ffstat = FIBRE_ONSET ; mfstat = * ; fstat = ffstat
    ! To 1:  ffstat = FIBRE_ONSET ; mfstat = * ; fstat = ffstat
    ! To 2:  ffstat = FIBRE_FAILED; mfstat = * ; fstat = ffstat
    !
    ! From:  ffstat = INTACT      ; mfstat = * ; fstat = ffstat
    ! To 1:  ffstat = INTACT      ; mfstat = * ; fstat = mfstat (mfstat>INTACT)
    ! To 2:  ffstat = FIBRE_ONSET ; mfstat = * ; fstat = ffstat
    ! To 3:  ffstat = FIBRE_FAILED; mfstat = * ; fstat = ffstat
    !---------------------------------------------------------------------------
    
    ! ffstat is the key variable, select what to do next based on ffstat
    ffstat_bfr_cohlaw: select case (ffstat)
    
      ! if fibres are already FAILED, just calculate stress and return;
      ! no need to go through failure criterion, coh law and update sdvs
      case (FIBRE_FAILED)
        ! calculate dee; degrade both fibre and matrix stiffness properties
        call deemat3d (dee, this_mat, df=df, dm2=df, dm3=df)
        ! calculate stress
        stress = matmul(dee, strain)
        ! exit program
        return
        
      ! when fibres are DAMAGED, no need to calculate stress; 
      ! only strain is needed for subsequent cohesive law calculations;
      ! so do NOTHING  
      case (FIBRE_ONSET)
        continue
        
      ! when fibres are INTACT, calculate stress for failure criterion check
      case (INTACT)
        ! calculate dee using original material properties only
        call deemat3d (dee, this_mat)
        ! calculate stress
        stress = matmul(dee, strain)
        
      ! unexpected value, update istat and emsg, then return
      case default
        istat = STAT_FAILURE
        emsg  = 'ffstat value is not recognized in ffstat_bfr_cohlaw, &
        & lamina_ddsdde, lamina_type_module!'
        return
        
    end select ffstat_bfr_cohlaw
      
      
    ! call fibre cohesive law, calculate fibre damage var.
    call fibre_cohesive_law (ffstat, df, u0, uf, stress, strain, &
  & clength, this_mat%strength, this_mat%fibreToughness,         &
  & maxdm_lcl, istat, emsg)
  
    ! check if the cohesive law is run successfully
    if(istat == STAT_FAILURE) return

    !---------------------------------------------------------------------------
    ! going through cohesive law, the possible outcomes of ffstat are:
    ! 1. intact -> intact : check matrix status (update fstat if mfstat changes)
    ! 2. intact -> onset  : update fstat, D, stress
    ! 3. onset  -> onset  : update fstat, D, stress
    ! 4. onset  -> failed : update fstat, D, stress
    !---------------------------------------------------------------------------
  
    ! ffstat is the key variable, select what to do next based on ffstat
    ffstat_aft_cohlaw: select case (ffstat)
    
      ! if fibres are DAMAGED/FAILED after the cohesive law, 
      ! update fstat (and to sdv), deemat and stress
      case (FIBRE_ONSET, FIBRE_FAILED)
        ! update fstat
        fstat = ffstat
        ! update to sdv (everything except mfstat)
        sdv%FSTAT  = fstat
        sdv%FFSTAT = ffstat
        sdv%DF     = df
        sdv%U0     = u0
        sdv%UF     = uf
        ! update D matrix
        call deemat3d (dee, this_mat, df=df, dm2=df, dm3=df) 
        ! update stress
        stress = matmul(dee, strain) 
        
      ! if fibres are STILL INTACT, check for matrix status;
      ! if matrix failure onset, update fstat (and sdv)
      case (INTACT)
        ! check for matrix failure
        if (mfstat == INTACT) then       
          ! go through failure criterion and update fstat
          call matrix_failure_criterion3d (mfstat, stress, this_mat%strength)
          ! update fstat if matrix reached failure onset
          if(mfstat == MATRIX_ONSET) then
            fstat = mfstat 
            ! update to sdv (only fstat and mfstat)
            sdv%FSTAT  = fstat
            sdv%MFSTAT = mfstat
            ! no need to update D matrix and stress
          end if
        end if
        
      ! unexpected value, update istat and emsg, then return
      case default
        istat = STAT_FAILURE
        emsg  = 'ffstat value is not recognized in ffstat_aft_cohlaw, &
        & lamina_ddsdde, lamina_type_module!'
        return
 
    end select ffstat_aft_cohlaw
   
    
  end subroutine ddsdde_lamina 



! the rest are private procedures used in ddsdde_lamina subroutine
! they can be considered as internal procedures of ddsdde_lamina

! for internal/private procedures, the checking of the validity of 
! intent(in) and (inout) dummy arguments can be spared

! but any logical construct in complex procedures should be complete 
! with istat and emsg. This is the case with fibre_cohesive_law


  pure subroutine check_sdv (sdv, istat, emsg)
  ! Purpose:
  ! to check the validity of the input solution-dependent variables
  
    type(LAMINA_SDV),         intent(in)  :: sdv
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    
    ! initialize intent out variables
    istat = STAT_SUCCESS
    emsg  = ''
    
    if (.not. (ZERO-SMALLNUM < sdv%DF .and. sdv%DF < ONE+SMALLNUM)) then
      istat = STAT_FAILURE
      emsg  = 'lamina DF must be between [0., 1.], lamina_type_module'
      return
    end if
    
    if (.not. (ZERO-SMALLNUM < sdv%U0)) then
      istat = STAT_FAILURE
      emsg  = 'lamina U0 must be >= 0., lamina_type_module'
      return
    end if
    
    if (.not. (ZERO-SMALLNUM < sdv%UF)) then
      istat = STAT_FAILURE
      emsg  = 'lamina UF must be >= 0., lamina_type_module'
      return
    end if
    
    select case (sdv%FSTAT)
      case (INTACT, MATRIX_ONSET, MATRIX_FAILED, FIBRE_ONSET, FIBRE_FAILED)
        continue
      case default
        istat = STAT_FAILURE
        emsg  = 'lamina FSTAT value is incorrect, lamina_type_module'
        return
    end select
    
    select case (sdv%FFSTAT)
      case (INTACT, FIBRE_ONSET, FIBRE_FAILED)
        continue
      case default
        istat = STAT_FAILURE
        emsg  = 'lamina FFSTAT value is incorrect, lamina_type_module'
        return
    end select
    
    select case (sdv%MFSTAT)
      case (INTACT, MATRIX_ONSET, MATRIX_FAILED)
        continue
      case default
        istat = STAT_FAILURE
        emsg  = 'lamina MFSTAT value is incorrect, lamina_type_module'
        return
    end select
    
  end subroutine check_sdv



  pure subroutine check_mat_prop (this, istat, emsg)
  ! Purpose:
  ! to check the validity of the input material properties
 
    type(lamina_type),        intent(in)  :: this
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    
    ! initialize intent out variables
    istat = STAT_SUCCESS
    emsg  = ''
    
    ! elastic moduli must be positive non-zero
    if (this%modulus%E1 < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina E1 must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%modulus%E2 < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina E2 must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%modulus%G12 < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina G12 must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%modulus%G23 < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina G23 must be greater than zero, lamina_type_module'
      return
    end if
  
    ! check on nu12 and nu23 are omitted; they can take on both + and - values
    
    ! strengths must be positive non-zero
    if (this%strength%Xt < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Xt must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%strength%Xc < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Xc must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%strength%Yt < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Yt must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%strength%Yc < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Yc must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%strength%Sl < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Sl must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%strength%St < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina St must be greater than zero, lamina_type_module'
      return
    end if

    
    ! matrix toughnesses and mixed-mode ratio must be positive non-zero
    if (this%matrixToughness%GmcI < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina GmcI must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%matrixToughness%GmcII < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina GmcII must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%matrixToughness%eta < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina eta must be greater than zero, lamina_type_module'
      return
    end if


    ! fibre toughnesses must be positive non-zero
    if (this%fibreToughness%GfcT < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina GfcT must be greater than zero, lamina_type_module'
      return
    end if
    
    if (this%fibreToughness%GfcC < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina GfcC must be greater than zero, lamina_type_module'
      return
    end if


  end subroutine check_mat_prop



  pure subroutine deemat3d (dee, this_mat, df, dm2, dm3)
  ! Purpose:
  ! to calculate local stiffness matrix D
  ! for 3D problem with the standard 6 strains
    
    ! dummy argument list:
    ! - dee         : local stiffness matrix,                 to update
    ! - this_mat    : passed-in material object,              passed-in
    ! - df          : fibre degradation,                      passed-in
    ! - dm2         : matrix degradation (in-plane, dir. 2),  passed-in
    ! - dm3         : matrix degradation (out-plane, dir. 3), passed-in
    real(DP),           intent(inout) :: dee(:,:)
    type(lamina_type),  intent(in)    :: this_mat
    real(DP),   optional, intent(in)  :: df, dm2, dm3

    
    !local var. list:
    ! - E1,...,G23 : elastic moduli, standard notations
    ! - del        : a common denominator for calculating dee
    ! - dm23       : matrix degradation (transverse, 2-3 plane)
    real(DP) :: E1, E2, E3, nu12, nu13, nu23, nu21, nu31, nu32, G12, G13, G23
    real(DP) :: del, dm23

    
    ! initialize local variables
    E1   = ZERO;  E2   = ZERO;  E3   = ZERO
    G12  = ZERO;  G13  = ZERO;  G23  = ZERO
    nu12 = ZERO;  nu13 = ZERO;  nu23 = ZERO
    nu21 = ZERO;  nu31 = ZERO;  nu32 = ZERO
    del  = ZERO;  dm23 = ZERO
    
    ! for private procedures used in the main procedure (ddsdde_lamina),
    ! the values of the input arguments are not checked. 
    ! they are assumed to be of valid values.
      
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
    G23 = G23 * (ONE-dm23)
    ! do not degrade below residual stiffness
    G23 = max(G23,RESIDUAL_STIFFNESS)
        
    
    ! calculate D matrix terms
    ! zero all terms first
    dee = ZERO
  
    del = ONE-nu12*nu21-nu13*nu31-nu23*nu32-TWO*nu21*nu32*nu13
    del = del/E1/E2/E3
    
    dee(1,1) = (ONE  - nu23 * nu32)/E2/E3/del
    dee(1,2) = (nu21 + nu23 * nu31)/E2/E3/del
    dee(1,3) = (nu31 + nu21 * nu32)/E2/E3/del
    dee(2,1) = dee(1,2)
    dee(2,2) = (ONE  - nu13 * nu31)/E1/E3/del
    dee(2,3) = (nu32 + nu12 * nu31)/E1/E3/del
    dee(3,1) = dee(1,3)
    dee(3,2) = dee(2,3)
    dee(3,3) = (ONE  - nu12 * nu21)/E1/E2/del
    dee(4,4) = G12
    dee(5,5) = G13
    dee(6,6) = G23
    
    
  end subroutine deemat


    
  pure subroutine fibre_cohesive_law(ffstat, df, u0, uf, stress, strain, &
  & clength, strength, fibreToughness, maxdm, istat, emsg) 
  ! Purpose:
  ! to update fibre failure status, stiffness degradation factor and 
  ! cohesive law variables according to a linear cohesive softening law;
  ! this is a complex subroutine, istat and emsg are required to flag an
  ! error when an unexpected logical case is met
  
    ! dummy argument list:
    ! - ffstat          : fibre failure status                  to update
    ! - df              : fibre stiffness degradation factor    to update
    ! - u0              : fibre cohesive law, u0                to update
    ! - uf              : fibre cohesive law, uf                to update
    ! - stress          : stress vector                         passed-in
    ! - strain          : strain vector                         passed-in
    ! - clength         : element characteristic length         passed-in
    ! - strength        : strength parameters of lamina         passed-in
    ! - fibreToughness  : fibre-direction toughness parameters  passed-in
    ! - maxdm           : maximum degradation factor            passed-in
    ! - istat           : status variable of this procedure     to output
    ! - emsg            : error message                         to output
    integer,                      intent(inout) :: ffstat
    real(DP),                     intent(inout) :: df, u0, uf
    real(DP),                     intent(in)    :: stress(:), strain(:)
    real(DP),                     intent(in)    :: clength
    type(lamina_strength),        intent(in)    :: strength
    type(lamina_fibreToughness),  intent(in)    :: fibreToughness
    real(DP),                     intent(in)    :: maxdm
    integer,                      intent(out)   :: istat
    character(len=MSGLENGTH),     intent(out)   :: emsg
    
    
    ! local variables:
    ! Xt      : fibre tensile     strength
    ! Xc      : fibre compressive strength
    ! GfcT    : fibre tensile     fracture toughness
    ! GfcC    : fibre compressive fracture toughness
    ! findex  : failure index of the failure criterion
    ! T0      : traction at failure onset
    ! u_eff   : effective displacement for cohesive law
    ! T_eff   : effective traction     for cohesive law
    ! df_tmp  : temporary df
    real(DP) :: Xt, Xc, GfcT, GfcC, findex, T0, u_eff, T_eff, df_tmp
    
    
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
    
    ! initialize intent(out) & local variables
    istat  = STAT_SUCCESS
    emsg   = ''
    Xt     = ZERO;   Xc    = ZERO
    GfcT   = ZERO;   GfcC  = ZERO
    findex = ZERO;   T0    = ZERO
    u_eff  = ZERO;   T_eff = ZERO
    df_tmp = ZERO
    
    ! for private procedures used in the main procedure (ddsdde_lamina),
    ! the values of the input arguments are not checked. 
    ! they are assumed to be of valid values.
    
    ! extract strength parameters
    Xt   = strength%Xt
    Xc   = strength%Xc
    
    ! extract toughness parameters
    GfcT = fibreToughness%GfcT
    GfcC = fibreToughness%GfcC 
    
    
    ! check and update ffstat and damage variables
    !
    ! possible changes of ffstat: 
    !
    ! INTACT -> failure criterion -> INTACT : no change
    !                             -> ONSET  : update ffstat, u0, uf and df
    !                             -> FAILED : update ffstat, u0, uf and df
    !
    ! ONSET  -> cohesive law      -> ONSET  : update ffstat and df
    !                             -> FAILED : update ffstat and df
    !
    ! so basically, select what to do based in ffstat: 
    ! - if it is still INTACT, go to failure criterion, calculate
    !   u0, uf and df if failure criterion is met;
    ! - if it is FAILURE ONSET, go to cohesive law, and update df
    ! - if it is other value, flag error stat & msg and return
    
    ffstat_select: select case (ffstat)
      
      ! if ffstat is still INTACT, then df, u0, uf have never been updated;
      ! failure criterion needs to be firstly applied on the stress.
      ! if failure onset is reached, calculate the cohesive law var. u0 and uf,
      ! then update ffstat, calculate df (update ffstat again if brittle);
      ! if failure onset is not reached, do nothing, exit
      case (INTACT) ffstat_select
      
        ! apply failure criterion and calculate the failure index
        findex = max(stress(1), ZERO) / Xt + abs( min(stress(1), ZERO) / Xc )
        
        ! check to see if the fibre failure onset criterion is reached;
        ! if so, update ffstat, calculate u0, uf and df;
        ! if not, do nothing
        findex_if: if (findex > ONE-SMALLNUM) then 
          ! update ffstat
          ffstat = FIBRE_ONSET
          ! calculate effective jump and traction at failure onset
          u_eff  = abs(strain(1)) * clength
          T_eff  = abs(stress(1))    
          ! adjust    effective jump and traction at failure onset, 
          ! with overshoot of failure index (ideally, findex = ONE)
          u0     = u_eff / findex
          T0     = T_eff / findex    
          ! calculate effective jump at final failure
          if (strain(1) > ZERO) then
            ! use tensile     fracture toughness
            uf   = TWO * GfcT / T0
          else
            ! use compressive fracture toughness
            uf   = TWO * GfcC / T0
          end if   
          ! calculate df
          if(uf < u0 + SMALLNUM) then 
            ! GfcT/GfcC too small, brittle failure
            ! update df and ffstat
            df     = maxdm
            ffstat = FIBRE_FAILED
          else
            ! df is calculated such that 
            ! the updated stress(1) will be Xt or -Xc
            df     = ONE - ONE / findex
          end if    
        end if findex_if
 
 
      ! if ffstat is already FIBRE_ONSET, then cohesive law var. u0 and uf must
      ! have already been defined. 
      ! the only thing to do here is to update df with respect to displacement,
      ! according to the already-defined cohesive law. 
      ! if df reaches maxdm, update ffstat to FIBRE_FAILED
      case (FIBRE_ONSET) ffstat_select
    
        ! calculate effective displacement
        u_eff = abs(strain(1)) * clength
        ! go to the cohesive law ONLY when u_eff is NONZERO
        if (u_eff > SMALLNUM) then
          ! use the defined linear cohesive softening law var. u0 and uf
          ! to calculate df_tmp
          df_tmp = (uf / u_eff) * (u_eff-u0) / (uf-u0)
          ! update df only if df_tmp is larger
          if (df_tmp > df) df = df_tmp
        end if
    
        ! check df and update ffstat if df reaches the max
        if (df > maxdm-SMALLNUM) then
            df     = maxdm
            ffstat = FIBRE_FAILED
        end if
      
      case default ffstat_select
        ! this case should never be reached; flag error status and msg
        istat = STAT_FAILURE
        emsg  = 'unexpected ffstat value in fibre_cohesive_law, & 
        &lamina_type_module!'
        return
      
    end select ffstat_select
    
  end subroutine fibre_cohesive_law


    
  pure subroutine matrix_failure_criterion3d (mfstat, stress, strength)
  ! Purpose:
  ! to implement a matrix failure criterion based on the stress and strength
  ! for 3D problems with standard 6 strains.
  ! at the moment, only in-plane tensile failure is considered;
  ! only in-plane stress components are used in the failure criterion.
  ! more complicated failure criterion can be implemented here in the future
  ! to include compressive failure and out-plane stress components, and  
  ! to output matrix crack angle w.r.t shell plane
    
    ! list of dummy arguments:
    ! mfstat    : matrix failure status variable      to update
    ! stress    : stress array                        passed-in
    ! strength  : strength parameters                 passed-in
    integer,               intent(inout) :: mfstat
    real(DP),              intent(in)    :: stress(:)
    type(lamina_strength), intent(in)    :: strength
    
    ! local variables list:
    ! Yt, Yc, Sl, St  : lamina strength parameters
    ! findex          : failure index of the failure criterion
    real(DP)  :: Yt, Yc, Sl, St
    real(DP)  :: findex

    
    ! initialize local variables
    Yt      = ZERO 
    Yc      = ZERO 
    Sl      = ZERO 
    St      = ZERO
    findex  = ZERO
    
    ! for private procedures used in the main procedure (ddsdde_lamina),
    ! the values of the input arguments are not checked. 
    ! they are assumed to be of valid values.
    
    ! extract strength parameters
    Yt = strength%Yt
    Yc = strength%Yc
    Sl = strength%Sl
    St = strength%St
    
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
    
    ! calculate the failure index for tensile failure
    ! matrix crack is assumed to be perpendicular to shell plane
    ! no out-plane stress components are considered

    findex = sqrt( (max(stress(2),ZERO)/Yt)**2 + (stress(4)/Sl)**2 )
    
    if(findex > ONE-SMALLNUM) mfstat = MATRIX_ONSET
    
    
  end subroutine matrix_failure_criterion
    
    
    
end module lamina_type_module
