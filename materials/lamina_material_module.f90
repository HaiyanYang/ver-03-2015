module lamina_material_module
!
!  Purpose:
!    define an object to represent a lamina material
!    with the associated procedures to empty, update, extract and display
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
use parameter_module, only : DP, ZERO, ONE, TWO, SMALLNUM, RESIDUAL_MODULUS, &
                           & MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, INTACT,  &
                           & MATRIX_ONSET, FIBRE_ONSET, FIBRE_FAILED

implicit none
private

! elastic moduli, standard notations
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

! fibre toughness properties
! mode I toughness, tensile (GfcT) and compressive (GfcC)
type, public :: lamina_fibreToughness
  real(DP) :: GfcT = ZERO, GfcC = ZERO 
end type

! lamina material object definition
type, public :: lamina_material
  type(lamina_modulus)          :: modulus
  type(lamina_strength)         :: strength
  type(lamina_fibreToughness)   :: fibreToughness
end type

! lamina material solution-dependent variables
! df     : fibre modulus degradation factor
! u0, uf : cohesive law parameters, initial & final failure displacements
! fstat  : generic failure status (= ffstat/mfstat, whichever is more severe)
! ffstat : fibre   failure status
! mfstat : matrix  failure status
type, public :: lamina_sdv
  real(DP) :: df    = ZERO,   u0     = ZERO,   uf     = ZERO
  integer  :: fstat = INTACT, ffstat = INTACT, mfstat = INTACT
end type

interface empty
  module procedure empty_lamina
end interface
    
interface update
  module procedure update_lamina
end interface

interface display
  module procedure display_lamina
end interface

interface display
    module procedure display_lamina_sdv
end interface

interface ddsdde
    module procedure ddsdde_lamina
end interface


public :: empty, update, display, ddsdde




contains




  pure subroutine empty_lamina(this)
  ! Purpose:
  ! to reset this lamina object's components into their default values (ZERO)
 
    type(lamina_material),intent(inout) :: this
    
    this%modulus         = lamina_modulus(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO)
    this%strength        = lamina_strength(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO)
    this%fibreToughness  = lamina_fibreToughness(ZERO,ZERO)

  end subroutine empty_lamina



  pure subroutine update_lamina(this, modulus, strength, fibreToughness)
  ! Purpose:
  ! to update this lamina object's components
    
   	type(lamina_material),              intent(inout) :: this
    type(lamina_modulus),         optional,intent(in) :: modulus
    type(lamina_strength),        optional,intent(in) :: strength
    type(lamina_fibreToughness),  optional,intent(in) :: fibreToughness
            
    if(present(modulus))          this%modulus          = modulus    
    if(present(strength))         this%strength         = strength       
    if(present(fibreToughness))   this%fibreToughness   = fibreToughness

  end subroutine update_lamina  
    
   
  
  subroutine display_lamina(this)
  ! Purpose:
  ! to display this lamina object's components on cmd window
  ! this is useful for debugging
 
    type(lamina_material), intent(in) :: this
    
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
    write(*,'(1X, A)') 'Display the modulus of the inquired lamina object :'
    write(*,display_fmt) 'lamina E1   is: ', this%modulus%E1 
    write(*,display_fmt) 'lamina E2   is: ', this%modulus%E2
    write(*,display_fmt) 'lamina G12  is: ', this%modulus%G12
    write(*,display_fmt) 'lamina G23  is: ', this%modulus%G23
    write(*,display_fmt) 'lamina nu12 is: ', this%modulus%nu12
    write(*,display_fmt) 'lamina nu23 is: ', this%modulus%nu23
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the strength of the inquired lamina object :'
    write(*,display_fmt) 'lamina Xt   is: ', this%strength%Xt 
    write(*,display_fmt) 'lamina Xc   is: ', this%strength%Xc
    write(*,display_fmt) 'lamina Yt   is: ', this%strength%Yt
    write(*,display_fmt) 'lamina Yc   is: ', this%strength%Yc
    write(*,display_fmt) 'lamina St   is: ', this%strength%St
    write(*,display_fmt) 'lamina Sl   is: ', this%strength%Sl
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the fibre toughness of the inquired lamina &
    &object :'
    write(*,display_fmt) 'lamina GfcT  is: ', this%fibreToughness%GfcT 
    write(*,display_fmt) 'lamina GfcC  is: ', this%fibreToughness%GfcC
    write(*,'(1X, A)') ''

  end subroutine display_lamina
  
  
  
  subroutine display_lamina_sdv(this_sdv)
  ! Purpose:
  ! to display this lamina_sdv's components on cmd window
  ! this is useful for debugging
  
    type(lamina_sdv), intent(in) :: this_sdv
  
    ! local variable
    character(len=20) :: display_fmt
    
    ! initialize local variable
    display_fmt = ''

    ! set display format for string and integer
    ! A for string, I for integer, 10 for width of the number
    display_fmt = '(1X, A, I10)' 
    
    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the inquired lamina SDVs :'
    write(*,display_fmt) 'lamina FSTAT  is: ', this_sdv%FSTAT 
    write(*,display_fmt) 'lamina FFSTAT is: ', this_sdv%FFSTAT 
    write(*,display_fmt) 'lamina MFSTAT is: ', this_sdv%MFSTAT 

    ! set display format for string and real
    ! A for string, ES for real (scientific notation)
    ! 10 is width, 3 is no. of digits aft decimal point
    ! note that for scientific real, ESw.d, w>=d+7
    display_fmt = '(1X, A, ES10.3)' 
    
    write(*,display_fmt) 'lamina DF     is: ', this_sdv%DF 
    write(*,display_fmt) 'lamina U0     is: ', this_sdv%U0
    write(*,display_fmt) 'lamina UF     is: ', this_sdv%UF
    write(*,'(1X, A)') ''
  
  end subroutine display_lamina_sdv
  
  

  pure subroutine ddsdde_lamina(this_mat, dee, stress, sdv, strain, clength, &
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
    ! - strain      : local strain vector,                passed-in
    ! - clength     : elem. characteristic length,        passed-in
    ! - istat       : status variable of this procedure   to output
    ! - emsg        : error message                       to output
    ! - d_max       : maximum degradation factor          passed-in (optional)
    type(lamina_material),    intent(in)    :: this_mat
    real(DP),                 intent(inout) :: dee(:,:)
    real(DP),                 intent(inout) :: stress(:)
    type(lamina_sdv),         intent(inout) :: sdv
    real(DP),                 intent(in)    :: strain(:)    
    real(DP),                 intent(in)    :: clength
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    real(DP),       optional, intent(in)    :: d_max
    
    ! local variables list:
    ! - fstat     : generic failure status
    ! - ffstat    : fibre   failure status
    ! - mfstat    : matrix  failure status
    ! - df        : fibre degradation
    ! - u0        : fibre cohesive law, failure onset displacement
    ! - uf        : fibre cohesive law, total failure displacement
    ! - d_max_lcl : local copy of d_max
    ! - findex    : failure index of a failure criterion
    integer  :: fstat, ffstat, mfstat
    real(DP) :: df, u0, uf
    real(DP) :: d_max_lcl
    real(DP) :: findex

    
    ! initialize intent(out) & local variables
    istat     = STAT_SUCCESS  ! default
    emsg      = ''
    fstat     = 0 
    ffstat    = 0 
    mfstat    = 0
    df        = ZERO 
    u0        = ZERO 
    uf        = ZERO 
    d_max_lcl = ZERO
    findex    = ZERO
    
    ! check validity of dummy arguments with intent(in) or (inout)
    
    ! check this_mat properties
    call check_mat_prop (this_mat, istat, emsg)
    if (istat == STAT_FAILURE) return
    
    ! check dee, stress and strain size
    if (.not. (size(dee(:,1)) == 6 .and. size(dee(1,:)) == 6)) then
      istat = STAT_FAILURE
      emsg  = 'size of dee is not supported for ddsdde_lamina, &
      &lamina_material_module!'
      return
    end if
    
    if (.not. (size(stress) == 6)) then
      istat = STAT_FAILURE
      emsg  = 'size of stress is not supported for ddsdde_lamina, &
      &lamina_material_module!'
      return
    end if
    
    if (.not. (size(strain) == 6)) then
      istat = STAT_FAILURE
      emsg  = 'size of strain is not supported for ddsdde_lamina, &
      &lamina_material_module!'
      return
    end if
    
    ! check sdv values
    call check_sdv (sdv, istat, emsg)
    if (istat == STAT_FAILURE) return
    
    ! check clength value
    if (.not. (clength > SMALLNUM)) then 
      istat = STAT_FAILURE
      emsg  = 'clength must be larger than zero for ddsdde_lamina, &
      &lamina_material_module!'
      return
    end if
    
    ! check d_max value, must be between ZERO and ONE
    if (present(d_max)) then
      if (.not. (ZERO-SMALLNUM < d_max .and. d_max < ONE+SMALLNUM)) then
        istat = STAT_FAILURE
        emsg  = 'd_max must be within [0., 1.] for ddsdde_lamina, &
        &lamina_material_module!'
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
    if (present(d_max))  then
      d_max_lcl = d_max
    else  
      ! default value is ONE              
      d_max_lcl = ONE
    end if
    
    !---------------------------------------------------------------------------
    ! fstat is always equal to ffstat, 
    ! unless ffstat = INTACT and mfstat > INTACT, then fstat = mfstat
    !
    ! Possible changes of ffstat, mfstat and fstat before&after cohesive law
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
    ! for calculations before the fibre cohesive law:
    ! - if the fibres are already failed, just calculate D and stress, 
    !   then exit the problem;
    ! - if the fibres are already damaged (but not yet failed), then the fibre
    !   cohesive law is needed to update the damage variables, D and stress, 
    !   based on the current strain values;
    ! - if the fibres are still intact, then the fibre & matrix failure criteria
    !   are needed to check if failure criteria are met based on the current 
    !   stress values; if so, update the damage variables, D and stress 
    !   according to the fibre cohesive law
    ffstat_bfr_cohlaw: select case (ffstat)
    
      ! if fibres are already FAILED, just calculate stress and return;
      ! no need to go through failure criterion, coh law and update sdvs
      case (FIBRE_FAILED)
        ! calculate dee; degrade both fibre and matrix stiffness properties
        call deemat_3d (this_mat, dee, df=df, dm2=df, dm3=df)
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
        call deemat_3d (this_mat, dee)
        ! calculate stress
        stress = matmul(dee, strain)
        
      ! unexpected value, update istat and emsg, then return
      case default
        istat = STAT_FAILURE
        emsg  = 'ffstat value is not recognized in ffstat_bfr_cohlaw, &
        & lamina_ddsdde, lamina_material_module!'
        return
        
    end select ffstat_bfr_cohlaw
      
      
    ! call fibre cohesive law, calculate fibre damage variables (ffstat, df, u0,
    ! uf) based on this material's properties, current stress and strain states
    ! and the element's characteristic length
    call fibre_cohesive_law (this_mat, ffstat, df, u0, uf, &
    & stress, strain, clength, d_max_lcl, istat, emsg)
  
    ! check if the cohesive law is run successfully;
    ! if not, return (in this case, emsg will show the reason of error)
    if (istat == STAT_FAILURE) return

    !---------------------------------------------------------------------------
    ! going through the cohesive law, the possible outcomes of ffstat are:
    ! 1. intact -> intact : check matrix status, and update fstat, sdv
    ! 2. intact -> onset  : update fstat, sdv, D, stress
    ! 3. onset  -> onset  : update fstat, sdv, D, stress
    ! 4. onset  -> failed : update fstat, sdv, D, stress
    !---------------------------------------------------------------------------
  
    ! ffstat is the key variable, select what to do next based on ffstat for
    ! the calculations after the fibre cohesive law:
    ! - if fibres are damaged or failed, regardless of its prior status 
    !   (intact or damaged), damage variables, D and stress should be updated;
    ! - if fibres remain intact, then matrix failure criterion needs to be
    !   checked, if it is met, fstat would be updated as mfstat; no need to 
    !   update any damage variable or D or stress, because matrix failure here
    !   does not lead to softening (note that in FNM, matrix failure is
    !   explicitly represented by a cohesive sub-element; only the matrix
    !   failure status is interested here)
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
        call deemat_3d (this_mat, dee, df=df, dm2=df, dm3=df) 
        ! update stress
        stress = matmul(dee, strain) 
        
      ! if fibres are STILL INTACT, check for matrix status;
      ! if matrix failure onset, update fstat (and sdv)
      case (INTACT)
        ! check for matrix failure
        if (mfstat == INTACT) then       
          ! go through failure criterion and update fstat
          call matrix_failure_criterion_3d (this_mat, stress, findex)
          ! update fstat if matrix reached failure onset (findex >= 1.)
          if (findex > ONE - SMALLNUM) then
            mfstat = MATRIX_ONSET
            fstat  = mfstat 
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
        & lamina_ddsdde, lamina_material_module!'
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
  
    type(lamina_sdv),         intent(in)  :: sdv
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    
    ! initialize intent out variables
    istat = STAT_SUCCESS
    emsg  = ''
    
    if (.not. (ZERO-SMALLNUM < sdv%DF .and. sdv%DF < ONE+SMALLNUM)) then
      istat = STAT_FAILURE
      emsg  = 'lamina DF must be between [0., 1.], lamina_material_module'
      return
    end if
    
    if (.not. (ZERO-SMALLNUM < sdv%U0)) then
      istat = STAT_FAILURE
      emsg  = 'lamina U0 must be >= 0., lamina_material_module'
      return
    end if
    
    if (.not. (ZERO-SMALLNUM < sdv%UF)) then
      istat = STAT_FAILURE
      emsg  = 'lamina UF must be >= 0., lamina_material_module'
      return
    end if
    
    select case (sdv%FSTAT)
      case (INTACT, MATRIX_ONSET, FIBRE_ONSET, FIBRE_FAILED)
        continue
      case default
        istat = STAT_FAILURE
        emsg  = 'lamina FSTAT value is incorrect, lamina_material_module'
        return
    end select
    
    select case (sdv%FFSTAT)
      case (INTACT, FIBRE_ONSET, FIBRE_FAILED)
        continue
      case default
        istat = STAT_FAILURE
        emsg  = 'lamina FFSTAT value is incorrect, lamina_material_module'
        return
    end select
    
    select case (sdv%MFSTAT)
      case (INTACT, MATRIX_ONSET)
        continue
      case default
        istat = STAT_FAILURE
        emsg  = 'lamina MFSTAT value is incorrect, lamina_material_module'
        return
    end select
    
  end subroutine check_sdv



  pure subroutine check_mat_prop (this, istat, emsg)
  ! Purpose:
  ! to check the validity of the input material properties
 
    type(lamina_material),    intent(in)  :: this
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    
    ! initialize intent out variables
    istat = STAT_SUCCESS
    emsg  = ''
    
    ! elastic moduli must be positive non-zero
    if (this%modulus%E1 < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina E1 must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%modulus%E2 < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina E2 must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%modulus%G12 < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina G12 must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%modulus%G23 < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina G23 must be greater than zero, lamina_material_module'
      return
    end if
  
    ! check on nu12 and nu23 are omitted; they can take on both + and - values
    
    ! strengths must be positive non-zero
    if (this%strength%Xt < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Xt must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%strength%Xc < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Xc must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%strength%Yt < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Yt must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%strength%Yc < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Yc must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%strength%Sl < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina Sl must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%strength%St < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina St must be greater than zero, lamina_material_module'
      return
    end if

    ! fibre toughnesses must be positive non-zero
    if (this%fibreToughness%GfcT < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina GfcT must be greater than zero, lamina_material_module'
      return
    end if
    
    if (this%fibreToughness%GfcC < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'lamina GfcC must be greater than zero, lamina_material_module'
      return
    end if


  end subroutine check_mat_prop



  pure subroutine deemat_3d (this_mat, dee, df, dm2, dm3)
  ! Purpose:
  ! to calculate local stiffness matrix D
  ! for 3D problem with the standard 6 strains
    
    ! dummy argument list:
    ! - this_mat    : passed-in material object,              pass arg.
    ! - dee         : local stiffness matrix,                 to update
    ! - df          : fibre degradation,                      passed-in
    ! - dm2         : matrix degradation (in-plane, dir. 2),  passed-in
    ! - dm3         : matrix degradation (out-plane, dir. 3), passed-in
    type(lamina_material), intent(in)    :: this_mat
    real(DP),              intent(inout) :: dee(:,:)
    real(DP),    optional, intent(in)    :: df, dm2, dm3

    
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
    if (present(df)) then    
      E1   = E1   * (ONE - df)
      nu12 = nu12 * (ONE - df)
      nu13 = nu13 * (ONE - df)      
      ! do not degrade below residual stiffness
      if (E1 < RESIDUAL_MODULUS + SMALLNUM) then
        E1   = RESIDUAL_MODULUS
        nu12 = ZERO
        nu13 = ZERO
      end if      
    end if

    
    ! apply in-plane matrix degradation if present
    if (present(dm2)) then    
      E2   = E2   * (ONE - dm2)
      nu21 = nu21 * (ONE - dm2)
      nu23 = nu23 * (ONE - dm2)
      G12  = G12  * (ONE - dm2)      
      ! do not degrade below residual stiffness
      if (E2 < RESIDUAL_MODULUS + SMALLNUM) then
        E2   = RESIDUAL_MODULUS
        nu21 = ZERO
        nu23 = ZERO
      end if      
      G12  = max(G12, RESIDUAL_MODULUS)      
      ! update matrix degradation for 2-3 plane
      dm23 = max(dm23, dm2)      
    end if

    
    ! apply out-plane matrix degradation if present
    if (present(dm3)) then    
      E3   = E3   * (ONE - dm3)
      nu31 = nu31 * (ONE - dm3)
      nu32 = nu32 * (ONE - dm3)
      G13  = G13  * (ONE - dm3)      
      ! do not degrade below residual stiffness
      if (E3 < RESIDUAL_MODULUS + SMALLNUM) then
        E3   = RESIDUAL_MODULUS
        nu31 = ZERO
        nu32 = ZERO
      end if      
      G13  = max(G13, RESIDUAL_MODULUS)      
      ! update matrix degradation for 2-3 plane
      dm23 = max(dm23, dm3)      
    end if

    
    ! apply transverse (2-3) matrix degradation 
    G23 = G23 * (ONE-dm23)
    ! do not degrade below residual stiffness
    G23 = max(G23,RESIDUAL_MODULUS)
        
    
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
    
    
  end subroutine deemat_3d


    
  pure subroutine fibre_cohesive_law(this_mat, ffstat, df, u0, uf, &
  & stress, strain, clength, d_max, istat, emsg) 
  ! Purpose:
  ! to update fibre failure status, stiffness degradation factor and 
  ! cohesive law variables according to a linear cohesive softening law;
  ! this is a complex subroutine, istat and emsg are required to flag an
  ! error when an unexpected logical case is met
  
    ! dummy argument list:
    ! - this_mat        : lamina material object                pass arg.
    ! - ffstat          : fibre failure status                  to update
    ! - df              : fibre stiffness degradation factor    to update
    ! - u0              : fibre cohesive law, u0                to update
    ! - uf              : fibre cohesive law, uf                to update
    ! - stress          : stress vector                         passed-in
    ! - strain          : strain vector                         passed-in
    ! - clength         : element characteristic length         passed-in
    ! - d_max           : maximum degradation factor            passed-in
    ! - istat           : status variable of this procedure     to output
    ! - emsg            : error message                         to output
    type(lamina_material),        intent(in)    :: this_mat
    integer,                      intent(inout) :: ffstat
    real(DP),                     intent(inout) :: df, u0, uf
    real(DP),                     intent(in)    :: stress(:), strain(:)
    real(DP),                     intent(in)    :: clength
    real(DP),                     intent(in)    :: d_max
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
    Xt   = this_mat%strength%Xt
    Xc   = this_mat%strength%Xc
    
    ! extract toughness parameters
    GfcT = this_mat%fibreToughness%GfcT
    GfcC = this_mat%fibreToughness%GfcC 
    
    
    ! check and update ffstat and damage variables
    !
    ! possible changes of ffstat: 
    !
    ! INTACT -> failure criterion -> INTACT : no change
    !                             -> ONSET  : update ffstat, u0, uf and df
    !                             -> FAILED : update ffstat, u0, uf and df
    !
    ! ONSET  -> cohesive law      -> ONSET  : update df
    !                             -> FAILED : update ffstat and df
    !
    ! so basically, select what to do based in ffstat: 
    ! - if it is still INTACT, go to failure criterion, calculate
    !   u0, uf and df if failure criterion is met;
    ! - if it is FAILURE ONSET, go to cohesive law, and update df
    !     - if df is close to ONE, update ffstat to be FAILED
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
          if (uf < u0 + SMALLNUM) then 
            ! GfcT/GfcC too small, brittle failure
            ! update df and ffstat
            df     = d_max
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
      ! if df reaches d_max, update ffstat to FIBRE_FAILED
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
        if (df > d_max-SMALLNUM) then
            df     = d_max
            ffstat = FIBRE_FAILED
        end if
      
      case default ffstat_select
        ! this case should never be reached; flag error status and msg
        istat = STAT_FAILURE
        emsg  = 'unexpected ffstat value in fibre_cohesive_law, & 
        &lamina_material_module!'
        return
      
    end select ffstat_select
    
  end subroutine fibre_cohesive_law


    
  pure subroutine matrix_failure_criterion_3d (this_mat, stress, findex)
  ! Purpose:
  ! to implement a matrix failure criterion based on the stress and strength
  ! for 3D problems with standard 6 strains.
  ! at the moment, only in-plane tensile failure is considered;
  ! only in-plane stress components are used in the failure criterion.
  ! more complicated failure criterion can be implemented here in the future
  ! to include compressive failure and out-plane stress components, and  
  ! to output matrix crack angle w.r.t shell plane
    
    ! list of dummy arguments:
    ! - this_mat  : lamina object       pass arg.
    ! - stress    : stress array        passed-in
    ! - findex    : failure index       to output
    type(lamina_material), intent(in)    :: this_mat
    real(DP),              intent(in)    :: stress(:)
    real(DP),              intent(out)   :: findex
    
    ! local variables list:
    ! - Yt, Yc, Sl, St      : lamina strength parameters
    ! - tau_n, tau_t, tau_l : normal and two shear tractions on the potential
    !                         matrix crack surface
    ! - sigma1/2/3          : standard normal stress components
    ! - tau12/13/23         : standard shear stress components
    real(DP)  :: Yt, Yc, Sl, St
    real(DP)  :: tau_n, tau_t, tau_l
    real(DP)  :: sigma1, sigma2, sigma3, tau12, tau13, tau23

    
    ! initialize local & intent(out) variables
    findex  = ZERO
    Yt      = ZERO 
    Yc      = ZERO 
    Sl      = ZERO 
    St      = ZERO
    tau_n   = ZERO
    tau_t   = ZERO
    tau_l   = ZERO
    sigma1  = ZERO
    sigma2  = ZERO
    sigma3  = ZERO
    tau_12  = ZERO
    tau_13  = ZERO
    tau_23  = ZERO
    
    
    ! for private procedures used in the main procedure (ddsdde_lamina),
    ! the values of the input arguments are not checked. 
    ! they are assumed to be of valid values.
    
    ! extract strength parameters
    Yt = this_mat%strength%Yt
    Yc = this_mat%strength%Yc
    Sl = this_mat%strength%Sl
    St = this_mat%strength%St
    
    ! extract stress components
    sigma1  = stress(1)
    sigma2  = stress(2)
    sigma3  = stress(3)
    tau_12  = stress(4)
    tau_13  = stress(5)
    tau_23  = stress(6)
    
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
    
    tau_n = max(sigma2, ZERO)
    tau_t = ZERO
    tau_l = tau_12

    findex = sqrt( (tau_n/Yt)**2 + (tau_t/St)**2 + (tau_l/Sl)**2 )
    
    
  end subroutine matrix_failure_criterion_3d
    
    
    
end module lamina_material_module
