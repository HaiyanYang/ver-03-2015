module cohesive_material_module
!
!  Purpose:
!    for 3D problems only
!
!    define a cohesive sdv object to represent the solution dependent variables
!    of a cohesive material constitutive law, with the associated procedures to
!    display and extract its components. All updates are done within this module
!
!    define a cohesive material object to represent a cohesive material
!    with the associated procedures to empty, set, and display
!    its components, to integrate local stiffness matrix, and
!    to update cohesive sdv
!
!    define an integration point object of this material which stores all the
!    necessary information for calculation and output,
!    with associated procedures to empty, update, extract and display
!
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    05/04/15  B. Y. Chen            Original code
!    11/04/15  B. Y. Chen            Added cohesive_ig_point object and its
!                                    assoc. procedures
!

! list out ALL the parameter used, no matter how long:
! DP                : precision of real number, Double Precision
! ZERO, ONE, TWO,
! HALF              : self explanatory
! SMALLNUM          : a very small real number, used when comparing two reals
! RESIDUAL_MODULUS  : residual modulus of the material at final failure
! MSGLENGTH         : length of error message
! STAT_SUCCESS      : integer value of a successful status
! STAT_FAILURE      : integer value of a failure    status
! INTACT            : integer value for generic intact state
! COH_MAT_ONSET     : integer value for cohesive material failure onset state
! COH_MAT_FAILED    : integer value for cohesive material total failure state
! NDIM              : no. of dimensions
! NST               : no. of stress (traction) terms
use parameter_module, only : DP, ZERO, ONE, TWO, HALF,              &
                           & SMALLNUM, RESIDUAL_MODULUS,            &
                           & MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                           & INTACT, COH_MAT_ONSET, COH_MAT_FAILED, &
                           & NDIM, NST => NST_COHESIVE


implicit none
private

!**** auxiliary objects ****
! no encapsulation, no type-bound procedures
! these objects are created to facilitate the definition of the
! cohesive material object
! these objects are NOT meant to be used outside of this module

! penalty stiffness, normal and two shear directions
type, public :: cohesive_modulus
  real(DP) :: Dnn = ZERO
  real(DP) :: Dtt = ZERO
  real(DP) :: Dll = ZERO
end type cohesive_modulus

! strengths for the three directions
type, public :: cohesive_strength
  real(DP) :: tau_nc = ZERO
  real(DP) :: tau_tc = ZERO
  real(DP) :: tau_lc = ZERO
end type cohesive_strength

! fracture toughness for the three directions
! Gnc, Gtc and Glc : toughness values, t and l dir. values may be different
! alpha            : mixed-mode law coefficient (power law)
type, public :: cohesive_toughness
  real(DP) :: Gnc   = ZERO
  real(DP) :: Gtc   = ZERO
  real(DP) :: Glc   = ZERO
  real(DP) :: alpha = ZERO
end type cohesive_toughness

!**** end auxiliary objects ****


!**** the rest are actual objects of the module ****
! with encapsulation and allowed procedures

! cohesive material solution-dependent variables
! dm     : cohesive modulus degradation factor
! u0, uf : cohesive law parameters, initial & final failure separations
! fstat  : failure status
type, public :: cohesive_sdv
  private
  real(DP) :: dm    = ZERO,  u0 = ZERO,  uf = ZERO
  integer  :: fstat = INTACT
end type cohesive_sdv
! associated procedures: extract, display
! all other operations on the object are performed within this module
! note: this means that sdv cannot be defined or modified outside of
! this module

! cohesive material object definition
type,public :: cohesive_material
  private
  type(cohesive_modulus)   :: modulus
  type(cohesive_strength)  :: strength
  type(cohesive_toughness) :: toughness
end type cohesive_material
! associated procedures: empty, set, display, ddsdde
! all other operations on the object are performed within this module

! cohesive material integration point object
! stores everything needed for the integration of cohesive material in elements
type, public :: cohesive_ig_point
  private
  real(DP), allocatable :: x(:)           ! physical coordinates
  real(DP), allocatable :: u(:)           ! displacement
  real(DP), allocatable :: traction(:)    ! traction for output
  real(DP), allocatable :: separation(:)  ! separation for output
  type(cohesive_sdv) :: converged_sdv ! sdv of last converged increment
  type(cohesive_sdv) :: iterating_sdv ! sdv of current iteration
end type cohesive_ig_point
! associated procedures: empty, display, update, extract
! all other operations on the object are performed within this module


interface empty
  module procedure empty_cohesive
  module procedure empty_cohesive_ig_point
end interface empty

interface set
  module procedure set_cohesive
end interface set

interface display
  module procedure display_cohesive
  module procedure display_cohesive_sdv
  module procedure display_cohesive_ig_point
end interface display

interface ddsdde
  module procedure ddsdde_cohesive
  module procedure ddsdde_cohesive_intact
end interface ddsdde

interface update
  module procedure update_cohesive_ig_point
end interface update

interface extract
  module procedure extract_cohesive_sdv
  module procedure extract_cohesive_ig_point
end interface extract


public :: empty, set, display, ddsdde, update, extract




contains




  pure subroutine extract_cohesive_sdv (sdv, dm, u0, uf, fstat)
  ! Purpose:
  ! to extract components of the cohesive sdv object

    type(cohesive_sdv), intent(in)  :: sdv
    real(DP), optional, intent(out) :: dm
    real(DP), optional, intent(out) :: u0
    real(DP), optional, intent(out) :: uf
    integer,  optional, intent(out) :: fstat

    if (present(dm))    dm    = sdv%dm
    if (present(u0))    u0    = sdv%u0
    if (present(uf))    uf    = sdv%uf
    if (present(fstat)) fstat = sdv%fstat

  end subroutine extract_cohesive_sdv



  subroutine display_cohesive_sdv (this_sdv)
  ! Purpose:
  ! to display this cohesive_sdv's components on cmd window
  ! this is useful for debugging

    type(cohesive_sdv), intent(in) :: this_sdv

    ! local variable to set the output format
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

    write(*,display_fmt) 'cohesive DM    is: ', this_sdv%dm
    write(*,display_fmt) 'cohesive U0    is: ', this_sdv%U0
    write(*,display_fmt) 'cohesive UF    is: ', this_sdv%UF
    write(*,'(1X, A)') ''

  end subroutine display_cohesive_sdv



  pure subroutine empty_cohesive (this)
  ! Purpose:
  ! to reset this object's components into their default values (ZERO)

    type(cohesive_material), intent(inout) :: this

    this%modulus   = cohesive_modulus   (ZERO, ZERO, ZERO)
    this%strength  = cohesive_strength  (ZERO, ZERO, ZERO)
    this%toughness = cohesive_toughness (ZERO, ZERO, ZERO, ZERO)

  end subroutine empty_cohesive



  pure subroutine set_cohesive (this, modulus, strength, toughness, &
  & istat, emsg)
  ! Purpose :
  ! to set this object's components during preprocessing before analysis;
  ! property check subroutine, error status and message are needed to
  ! flag an error when the input properties are unacceptable

    ! - istat       : status variable of this procedure   to output
    ! - emsg        : error message                       to output
    type(cohesive_material),  intent(inout) :: this
    type(cohesive_modulus),   intent(in)    :: modulus
    type(cohesive_strength),  intent(in)    :: strength
    type(cohesive_toughness), intent(in)    :: toughness
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg

    ! local copy of intent(inout) variable
    type(cohesive_material) :: this_lcl

    ! initialize intent(out) & local variables
    istat = STAT_SUCCESS  ! default
    emsg  = ''

    this_lcl%modulus   = modulus
    this_lcl%strength  = strength
    this_lcl%toughness = toughness

    ! check this_mat properties
    call check_mat_prop (this_lcl, istat, emsg)
    if (istat == STAT_FAILURE) return

    ! update to dummy arg if inputs are valid
    this = this_lcl


    contains


      pure subroutine check_mat_prop (this, istat, emsg)
      ! Purpose:
      ! to check the validity of the input material properties

        type(cohesive_material),  intent(in)  :: this
        integer,                  intent(out) :: istat
        character(len=MSGLENGTH), intent(out) :: emsg

        ! initialize intent out variables
        istat = STAT_SUCCESS
        emsg  = ''

        ! elastic moduli must be positive non-zero
        if (this%modulus%Dnn < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive Dnn must be greater than zero, cohesive_material_module'
          return
        end if

        if (this%modulus%Dtt < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive Dtt must be greater than zero, cohesive_material_module'
          return
        end if

        if (this%modulus%Dll < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive Dll must be greater than zero, cohesive_material_module'
          return
        end if


        ! strengths must be positive non-zero
        if (this%strength%tau_nc < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive tau_nc must be greater than zero, &
          &cohesive_material_module'
          return
        end if

        if (this%strength%tau_tc < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive tau_tc must be greater than zero, &
          &cohesive_material_module'
          return
        end if

        if (this%strength%tau_lc < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive tau_lc must be greater than zero, &
          &cohesive_material_module'
          return
        end if


        ! matrix toughnesses and mixed-mode ratio must be positive non-zero
        if (this%toughness%Gnc < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive Gnc must be greater than zero, cohesive_material_module'
          return
        end if

        if (this%toughness%Gtc < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive Gtc must be greater than zero, cohesive_material_module'
          return
        end if

        if (this%toughness%Glc < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive Glc must be greater than zero, cohesive_material_module'
          return
        end if

        if (this%toughness%alpha < SMALLNUM) then
          istat = STAT_FAILURE
          emsg  = 'cohesive alpha must be greater than zero, cohesive_material_module'
          return
        end if


      end subroutine check_mat_prop



  end subroutine set_cohesive



  subroutine display_cohesive (this)
  ! Purpose:
  ! to display this cohesive object's components on cmd window
  ! this is useful for debugging

    type(cohesive_material), intent(in) :: this

    ! local variable to set the output format
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
    write(*,display_fmt) 'cohesive alpha  is: ', this%toughness%alpha
    write(*,'(1X, A)') ''

  end subroutine display_cohesive



  pure subroutine ddsdde_cohesive_intact (this_mat, dee, traction, separation)
  ! Purpose:
  ! to use intact material properties to calculate the D matrix and traction
  ! at an integration point of an element of cohesive_material definition
  ! (restricted to 3D problems, with the standard 3 separation terms)

    ! dummy argument list:
    ! - this_mat    : material properties object,         passed-in (pass arg)
    ! - dee         : local stiffness matrix D,           to update
    ! - traction    : local traction vector,              to update
    ! - separation  : local separation vector,            passed-in
    type(cohesive_material),  intent(in)    :: this_mat
    real(DP),                 intent(inout) :: dee(NST,NST)
    real(DP),                 intent(inout) :: traction(NST)
    real(DP),                 intent(in)    :: separation(NST)

    !**** check validity of non-pass dummy arguments with intent(in/inout) ****
    ! they include : dee, traction and separation

    ! dee and traction input values are not used; they can be intent(out).
    ! they are defined as intent(inout) to avoid any potential memory leak.
    ! so no need to check their input values

    ! separation components can take any real value, nothing to check

    !**** MAIN CALCULATIONS ****

    ! calculate dee using original material properties only
    call deemat_3d (this_mat, dee)

    ! calculate traction
    traction = matmul(dee, separation)

  end subroutine ddsdde_cohesive_intact



  pure subroutine ddsdde_cohesive (this_mat, dee, traction, sdv, separation, &
  & istat, emsg, d_max)
  ! Purpose:
  ! to calculate the D matrix, traction and solution-dependent variables
  ! at an integration point of an element of cohesive_material definition
  ! (restricted to 3D problems, with the standard 3 separation terms)

  ! note: this is a complex outbound procedure with istat and emsg,
  ! local copies of intent(inout) dummy args are used for calculations

    ! dummy argument list:
    ! - this_mat    : material properties object,         passed-in (pass arg)
    ! - dee         : local stiffness matrix D,           to update
    ! - traction    : local traction vector,              to update
    ! - sdv         : solution-dependent variables array, to update
    ! - separation  : local separation vector,            passed-in
    ! - istat       : status variable of this procedure   to output
    ! - emsg        : error message                       to output
    ! - d_max       : maximum degradation factor          passed-in (optional)

    type(cohesive_material),  intent(in)    :: this_mat
    real(DP),                 intent(inout) :: dee(NST,NST)
    real(DP),                 intent(inout) :: traction(NST)
    type(cohesive_sdv),       intent(inout) :: sdv
    real(DP),                 intent(in)    :: separation(NST)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    real(DP),       optional, intent(in)    :: d_max

    ! local variables list:
    !** local copies of intent(inout) dummy args:
    ! - dee_lcl         : local copy of dee
    ! - traction_lcl    : local copy of traction
    ! - fstat           : local copy of sdv failure status
    ! - dm              : local copy of sdv degradation factor
    ! - u0              : local copy of sdv failure onset displacement
    ! - uf              : local copy of sdv total failure displacement
    !** local copies of optional dummy args
    ! - d_max_lcl       : local copy of d_max
    !** purely local variables:
    ! - is_closed_crack : logical variable, true if the cohesive material is
    !                     under compression
    real(DP) :: dee_lcl(NST,NST)
    real(DP) :: traction_lcl(NST)
    integer  :: fstat
    real(DP) :: dm, u0, uf
    real(DP) :: d_max_lcl
    logical  :: is_closed_crack


    !**** initialize intent(out) & local variables ****
    istat     = STAT_SUCCESS  ! default
    emsg      = ''
    dee_lcl   = ZERO
    traction_lcl = ZERO
    fstat     = 0
    dm        = ZERO
    u0        = ZERO
    uf        = ZERO
    d_max_lcl = ZERO
    is_closed_crack = .false. ! default

    !**** check validity of non-pass dummy arguments with intent(in/inout) ****
    ! they include : dee, traction, sdv, separation, d_max

    ! dee and traction input values are not used; they can be intent(out).
    ! they are defined as intent(inout) to avoid any potential memory leak.
    ! so no need to check their input values

    ! sdv objects are only modified within this module, so no need to check

    ! separation components can take any real value, nothing to check

    ! check d_max value, must be between ZERO and ONE
    if (present(d_max)) then
      if (.not. (ZERO-SMALLNUM < d_max .and. d_max < ONE+SMALLNUM)) then
        istat = STAT_FAILURE
        emsg  = 'd_max must be within [0., 1.] for ddsdde_cohesive, &
        &cohesive_material_module!'
        return
      end if
    end if


    !**** assign values to local variables ****

    ! dee and traction's passed-in values are not needed; no need to pass
    ! them to their local copies

    ! extract values of sdv's failure status and cohesive law variables
    ! and store them in their respective local copies
    fstat  = sdv%FSTAT
    dm     = sdv%DM
    u0     = sdv%U0
    uf     = sdv%UF

    ! copy max. degradation (if present) into local variable
    if(present(d_max))  then
      d_max_lcl = d_max
    else
      ! default value is ONE
      d_max_lcl = ONE
    end if

    ! if normal separation is less than ZERO, then the crack is closed
    is_closed_crack = separation(1) < (ZERO - SMALLNUM)


    !**** MAIN CALCULATIONS ****

    !---------------------------------------------------------------------------
    ! Possible changes of fstat before&after cohesive law
    !
    ! From:  fstat = COH_MAT_FAILED
    ! To  :  fstat = COH_MAT_FAILED (NO CHANGE)
    !
    ! From:  fstat = COH_MAT_ONSET
    ! To 1:  fstat = COH_MAT_ONSET
    ! To 2:  fstat = COH_MAT_FAILED
    !
    ! From:  fstat = INTACT
    ! To 1:  fstat = INTACT
    ! To 2:  fstat = COH_MAT_ONSET
    ! To 3:  fstat = COH_MAT_FAILED
    !---------------------------------------------------------------------------

    ! fstat is the key variable, select what to do next based on fstat
    ! for calculations before the cohesive law:
    ! - if the cohesive material is already failed, just calculate D and
    !   traction, then update inout dummy args and exit the program;
    ! - if the cohesive material is already damaged but not yet failed, then the
    !   cohesive law is needed to update the damage variables, D and traction,
    !   based on the current separation vector;
    ! - if the cohesive material is still intact, then the cohesive failure
    !   criterion is needed to check if failure onset is met based on the
    !   current traction values; if so, update the damage variables, D and
    !   traction according to the cohesive law

    fstat_bfr_cohlaw: select case (fstat)

      ! if the coh mat is already FAILED, just calculate D and traction, then
      ! update inout dummy args and exit the program;
      ! no need to go through failure criterion, coh law and update sdvs
      case (COH_MAT_FAILED)
        ! calculate dee; degrade stiffness
        call deemat_3d (this_mat, dee_lcl, dm, is_closed_crack)
        ! calculate traction
        traction_lcl = matmul(dee_lcl, separation)

        !**** update intent(inout) dummy args before successful return ****
        dee      = dee_lcl
        traction = traction_lcl
        ! sdv remains unchanged
        ! exit program
        return

      ! when the coh mat is DAMAGED, no need to calculate traction;
      ! only separation is needed for subsequent cohesive law calculations;
      ! so do NOTHING
      case (COH_MAT_ONSET)
        continue

      ! when the coh mat is INTACT, calculate traction for failure criterion
      case (INTACT)
        ! calculate dee using original material properties only
        call deemat_3d (this_mat, dee_lcl)
        ! calculate traction
        traction_lcl = matmul(dee_lcl, separation)

      ! unexpected value, update istat and emsg and intent(out) dummy args
      ! then return
      case default
        istat = STAT_FAILURE
        emsg  = 'fstat value is not recognized in fstat_bfr_cohlaw, &
        & cohesive_ddsdde, cohesive_material_module!'
        return

    end select fstat_bfr_cohlaw



    ! call cohesive law, calculate damage variables (fstat, dm, u0, uf) based
    ! on this material's properties and current traction and separation vectors
    call cohesive_law (this_mat, fstat, dm, u0, uf, traction_lcl, separation, &
    & d_max_lcl, istat, emsg)

    ! check if the cohesive law is run successfully;
    ! if not, return (in this case, emsg will show the cause of error)
    ! note that if other intent(out) dummy args are present, they need to be
    ! assigned error-exit values before return
    if (istat == STAT_FAILURE) return


    !---------------------------------------------------------------------------
    ! after the cohesive law, the damage variables are already updated.
    ! the possible actions based on the change of fstat are:
    ! 1. intact -> intact : no change, do nothing
    ! 2. intact -> onset  : update D, traction
    ! 3. onset  -> onset  : update D, traction
    ! 4. onset  -> failed : update D, traction
    ! then update inout dummy args and exit the program
    !---------------------------------------------------------------------------

    ! fstat is the key variable, select what to do next based on fstat for
    ! the calculations after the cohesive law:
    ! - if the coh mat is damaged or failed, regardless of its prior status
    !   (intact or damaged), D and traction should be updated; then update
    !   inout dummy args and exit the program
    ! - if the coh mat remains intact, then the D and traction have already been
    !   calculated before the cohesive law, and sdv remain unchanged;
    !   so just update inout dummy args and exit the program
    fstat_aft_cohlaw: select case (fstat)

      ! if the coh mat is DAMAGED/FAILED after the cohesive law,
      ! update deemat and traction, then update all inout dummy args and exit
      case (COH_MAT_ONSET, COH_MAT_FAILED)
        ! update D matrix
        call deemat_3d (this_mat, dee_lcl, dm, is_closed_crack)
        ! update traction
        traction_lcl = matmul(dee_lcl, separation)

        !**** update intent(inout) dummy args before successful return ****
        dee       = dee_lcl
        traction  = traction_lcl
        sdv%FSTAT = fstat
        sdv%DM    = dm
        sdv%U0    = u0
        sdv%UF    = uf
        ! exit the program
        return

      ! if the coh mat is STILL INTACT, D and traction has already
      ! been calculated before the cohesive law and damage var. remain unchanged
      ! just update inout dummy args and exit the program
      case (INTACT)
        !**** update intent(inout) dummy args before successful return ****
        dee       = dee_lcl
        traction  = traction_lcl
        ! sdv remain unchanged
        ! exit the program
        return

      ! unexpected value, update istat and emsg and intent(out) dummy args
      ! then return
      case default
        istat = STAT_FAILURE
        emsg  = 'fstat value is not recognized in fstat_aft_cohlaw, &
        & cohesive_ddsdde, cohesive_material_module!'
        return

    end select fstat_aft_cohlaw


  end subroutine ddsdde_cohesive


! the rest are private procedures used in ddsdde_cohesive subroutine
! they can be considered as internal procedures of ddsdde_cohesive

! for internal/private procedures, the checking of the validity of
! intent(in) and (inout) dummy arguments can be spared

! but any logical construct in complex procedures should be complete
! with istat and emsg. This is the case with cohesive_law


  pure subroutine deemat_3d (this_mat, dee, dm, is_closed_crack)
  ! Purpose:
  ! to calculate local stiffness matrix D
  ! for 3D problem with the standard 3 displacement separations

    ! dummy argument list:
    ! - this_mat        : passed-in material object,     pass arg.
    ! - dee             : local stiffness matrix,        to update
    ! - dm              : degradation factor,            passed-in (optional)
    ! - is_closed_crack : true if the crack is closed    passed-in (optional)
    type(cohesive_material), intent(in)    :: this_mat
    real(DP),                intent(inout) :: dee(:,:)
    real(DP),      optional, intent(in)    :: dm
    logical,       optional, intent(in)    :: is_closed_crack


    !local var. list:
    ! - Dnn, Dtt, Dll : penalty stiffness, standard notations
    ! - crack_closed  : local copy of is_closed_crack
    real(DP) :: Dnn, Dtt, Dll
    logical  :: crack_closed


    ! initialize local variables
    Dnn  = ZERO
    Dtt  = ZERO
    Dll  = ZERO
    crack_closed = .false.

    ! for private procedures used in the main procedure (ddsdde_cohesive),
    ! the values of the input arguments are not checked.
    ! they are assumed to be of valid values.

    ! extract elastic muduli from passed-in material object
    Dnn  = this_mat%modulus%Dnn
    Dtt  = this_mat%modulus%Dtt
    Dll  = this_mat%modulus%Dll

    ! extract info about if the crack is closed
    if (present(is_closed_crack)) crack_closed = is_closed_crack

    ! apply fibre degradation if dm is present
    if (present(dm)) then
      ! degrade normal stiffness if crack is NOT closed
      if (.not. crack_closed) Dnn = Dnn * (ONE - dm)
      Dtt = Dtt * (ONE - dm)
      Dll = Dll * (ONE - dm)
      ! do not degrade below residual stiffness
      if (Dnn < RESIDUAL_MODULUS + SMALLNUM) Dnn = RESIDUAL_MODULUS
      if (Dtt < RESIDUAL_MODULUS + SMALLNUM) Dtt = RESIDUAL_MODULUS
      if (Dll < RESIDUAL_MODULUS + SMALLNUM) Dll = RESIDUAL_MODULUS
    end if

    ! calculate D matrix terms
    ! zero all terms first
    dee = ZERO

    dee(1,1) = Dnn
    dee(2,2) = Dtt
    dee(3,3) = Dll


  end subroutine deemat_3d



  pure subroutine cohesive_law (this_mat, fstat, dm, u0, uf, traction, &
  & separation, d_max, istat, emsg)
  ! Purpose:
  ! to update cohesive failure status, stiffness degradation factor and
  ! cohesive law variables according to a linear cohesive softening law;
  ! this is a complex subroutine, istat and emsg are required to flag an
  ! error when an unexpected logical case is met

    ! dummy argument list:
    ! - this_mat        : cohesive material object              pass arg.
    ! - fstat           : failure status                        to update
    ! - dm              : stiffness degradation factor          to update
    ! - u0              : cohesive law, u0                      to update
    ! - uf              : cohesive law, uf                      to update
    ! - traction        : traction   vector                     passed-in
    ! - separation      : separation vector                     passed-in
    ! - d_max           : maximum   degradation factor          passed-in
    ! - istat           : status variable of this procedure     to output
    ! - emsg            : error message                         to output
    type(cohesive_material),  intent(in)    :: this_mat
    integer,                  intent(inout) :: fstat
    real(DP),                 intent(inout) :: dm, u0, uf
    real(DP),                 intent(in)    :: traction(:), separation(:)
    real(DP),                 intent(in)    :: d_max
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg

    ! local variables list:
    ! - tau_n/t/l                 : normal and two shear tractions
    ! - delta_n/t/l               : normal and two shear separations
    ! - findex                    : failure index
    ! - Gnc/tc/lc, alpha          : cohesive toughnesses and power law exponent
    ! - Gn/t/l                    : normal and two shear strain energy density
    ! - lambda_n/t/l              : normal and two shear mode ratios
    ! - Gmc                       : mixed-mode fracture toughness
    ! - u_eff, T_eff              : effective separation and traction
    ! - T0                        : effective traction at failure onset
    ! - dm_tmp                    : temporary dm
    real(DP) :: tau_n,   tau_t,   tau_l
    real(DP) :: delta_n, delta_t, delta_l
    real(DP) :: findex
    real(DP) :: Gnc, Gtc, Glc, alpha
    real(DP) :: Gn,  Gt,  Gl
    real(DP) :: lambda_n, lambda_t, lambda_l
    real(DP) :: Gmc, u_eff, T_eff, T0, dm_tmp

    ! initialize intent(out) & local variables
    istat    = STAT_SUCCESS
    emsg     = ''
    tau_n    = ZERO;  tau_t    = ZERO;  tau_l    = ZERO
    delta_n  = ZERO;  delta_t  = ZERO;  delta_l  = ZERO
    findex   = ZERO
    Gnc      = ZERO;  Gtc      = ZERO;  Glc      = ZERO
    alpha    = ZERO
    Gn       = ZERO;  Gt       = ZERO;  Gl       = ZERO
    lambda_n = ZERO;  lambda_t = ZERO;  lambda_l = ZERO
    Gmc      = ZERO
    u_eff    = ZERO;  T_eff    = ZERO;  T0       = ZERO
    dm_tmp   = ZERO


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


    ! extract traction terms
    tau_n   = max(traction(1), ZERO)
    tau_t   = traction(2)
    tau_l   = traction(3)

    ! extract separation terms
    delta_n = max(separation(1), ZERO)
    delta_t = separation(2)
    delta_l = separation(3)

    ! extract toughness parameters
    Gnc    = this_mat%toughness%Gnc
    Gtc    = this_mat%toughness%Gtc
    Glc    = this_mat%toughness%Glc
    alpha  = this_mat%toughness%alpha

    ! check and update fstat and damage variables
    !
    ! possible changes of fstat:
    !
    ! INTACT -> failure criterion -> INTACT : no change
    !                             -> ONSET  : update fstat, u0, uf and dm
    !                             -> FAILED : update fstat, u0, uf and dm
    !
    ! ONSET  -> cohesive law      -> ONSET  : update dm
    !                             -> FAILED : update fstat and dm
    !
    ! so basically, select what to do based in fstat:
    ! - if it is still INTACT, go to failure criterion, calculate
    !   u0, uf and dm if failure criterion is met;
    ! - if it is FAILURE ONSET, go to cohesive law, and update dm
    !     - if dm is close to ONE, update fstat to be FAILED
    ! - if it is other value, flag error stat & msg and return

    fstat_select : select case (fstat)

      ! if fstat is still INTACT, then dm, u0, uf have never been updated;
      ! failure criterion needs to be firstly applied on the traction.
      ! if failure onset is reached, calculate the cohesive law var. u0 and uf,
      ! then update fstat, calculate dm (update fstat again if brittle);
      ! if failure onset is not reached, do nothing, exit
      case (INTACT) fstat_select

        ! apply failure criterion and calculate failure index (for 3D problems)
        call failure_criterion_3d (this_mat, traction, findex)

        ! check to see if the failure onset criterion is reached;
        ! if so, update fstat, calculate u0, uf and dm;
        ! if not, do nothing
        findex_if: if (findex > ONE-SMALLNUM) then
          ! update fstat
          fstat = COH_MAT_ONSET
          ! the rest calculates u0 and uf based on power law
          ! normal tension energy density
          Gn = HALF * tau_n * delta_n
          ! transverse shear energy density
          Gt = HALF * tau_t * delta_t
          ! longitudinal shear energy density
          Gl = HALF * tau_l * delta_l
          ! respective mode ratios
          lambda_n = Gn / (Gn + Gt + Gl)
          lambda_t = Gt / (Gn + Gt + Gl)
          lambda_l = Gl / (Gn + Gt + Gl)
          ! effective separation (could not be ZERO)
          u_eff = sqrt( delta_n**2 + delta_t**2 + delta_l**2 )
          ! effective traction is calculated such that
          ! HALF * u_eff * T_eff = total strain energy density
          T_eff = TWO * (Gn + Gt + Gl) / u_eff
          ! effective separation and traction at failure onset,
          ! adjusted with overshoot
          u0 = u_eff / findex
          T0 = T_eff / findex
          ! mixed mode fracture toughness (power law)
          Gmc = (lambda_n/Gnc)**alpha + &
              & (lambda_t/Gtc)**alpha + (lambda_l/Glc)**alpha
          Gmc = ONE / ( Gmc**(ONE/alpha) )
          ! effective separation at final failure
          uf  = two * Gmc / T0
          ! calculate dm
          if (uf < u0 + SMALLNUM) then
            ! Gmc too small, brittle failure
            ! update dm and fstat
            dm    = d_max
            fstat = COH_MAT_FAILED
          else
            ! dm is calculated such that
            ! the updated tractions just meet the failure criterion
            dm    = ONE - ONE / findex
          end if
        end if findex_if


      ! if fstat is already COH_MAT_ONSET, then cohesive law var. u0 and uf must
      ! have already been defined.
      ! the only thing to do here is to update dm with respect to separation,
      ! according to the already-defined cohesive law.
      ! if dm reaches d_max, update fstat to COH_MAT_FAILED
      case (COH_MAT_ONSET) fstat_select

        ! calculate effective separation
        u_eff = sqrt( delta_n**2 + delta_t**2 + delta_l**2 )
        ! go to the cohesive law ONLY when u_eff is NONZERO
        if (u_eff > SMALLNUM) then
          ! use the defined linear cohesive softening law var. u0 and uf
          ! to calculate dm_tmp
          dm_tmp = (uf / u_eff) * (u_eff-u0) / (uf-u0)
          ! update dm only if dm_tmp is larger
          if (dm_tmp > dm) dm = dm_tmp
        end if
        ! check dm and update fstat if dm reaches the max
        if (dm > d_max-SMALLNUM) then
            dm    = d_max
            fstat = COH_MAT_FAILED
        end if

      case default fstat_select
        ! this case should never be reached; flag error status and msg
        istat = STAT_FAILURE
        emsg  = 'unexpected fstat value in cohesive_law, &
        &cohesive_material_module!'
        return

    end select fstat_select


  end subroutine cohesive_law



  pure subroutine failure_criterion_3d (this_mat, traction, findex)
  ! Purpose:
  ! to calculate the failure index of the coh mat under current traction
  ! for 3D problem with the standard 3 separations
  ! findex >= ONE when failure onset is met

    ! dummy argument list:
    ! - this_mat    : passed-in material object,              pass arg.
    ! - traction    : current traction vector                 passed-in
    ! - findex      : failure index                           to output
    type(cohesive_material),  intent(in)  :: this_mat
    real(DP),                 intent(in)  :: traction(:)
    real(DP),                 intent(out) :: findex

    ! local variables
    ! - tau_n/t/l    : normal and two shear tractions
    ! - tau_nc/tc/lc : normal and two shear strengths
    real(DP) :: tau_n,  tau_t,  tau_l
    real(DP) :: tau_nc, tau_tc, tau_lc

    ! initialize intent(out) & local variables
    findex = ZERO
    tau_n  = ZERO;  tau_t  = ZERO;  tau_l  = ZERO
    tau_nc = ZERO;  tau_tc = ZERO;  tau_lc = ZERO

    ! extract traction terms
    tau_n  = max(traction(1), ZERO)
    tau_t  = traction(2)
    tau_l  = traction(3)

    ! extract strength parameters
    tau_nc = this_mat%strength%tau_nc
    tau_tc = this_mat%strength%tau_tc
    tau_lc = this_mat%strength%tau_lc

    ! --------------------------------------------------------- !
    ! the following assumes quadratic traction criterion
    ! other criteria can also be used in the future
    ! just need to put in a selection criterion and algorithm
    ! --------------------------------------------------------- !
    ! e.g.: select case (strength%FC)
    !           case(0)
    !               quad. traction
    !           case(1)
    !               max. traction
    !           ...

    findex = sqrt( (tau_n/tau_nc)**2 + (tau_t/tau_tc)**2 + (tau_l/tau_lc)**2 )

    ! if findex > ONE, then the current tractions overshoot the strengths
    ! to back-calculate the tractions where the failure criterion is just met,
    ! i.e., where findex = ONE, the current tractions need to be divided by
    ! findex

  end subroutine failure_criterion_3d




! the rest are standard procedures associated with the cohesive_ig_point object



  pure subroutine empty_cohesive_ig_point (ig_point)

    type(cohesive_ig_point), intent(inout) :: ig_point

    ! local varibale; derived type is initialized upon declaration
    type(cohesive_ig_point) :: ig_point_lcl

    ! reset the dummy arg. to initial values
    ig_point = ig_point_lcl

  end subroutine empty_cohesive_ig_point



  pure subroutine update_cohesive_ig_point (ig_point, x, u, &
  & traction, separation, converged_sdv, iterating_sdv)
  ! Purpose:
  ! to update cohesive ig point components
  ! this is an outbound procedure, so its inputs should be checked for validity
  ! x, u, traction, separation can take any value
  ! cohesive sdv objects are only modified within this module
  ! so, checking is spared

    type(cohesive_ig_point),      intent(inout) :: ig_point
    real(DP),           optional, intent(in)    :: x(NDIM), u(NDIM)
    real(DP),           optional, intent(in)    :: traction(NST)
    real(DP),           optional, intent(in)    :: separation(NST)
    type(cohesive_sdv), optional, intent(in)    :: converged_sdv
    type(cohesive_sdv), optional, intent(in)    :: iterating_sdv

    if(present(x)) then
      if (.not. allocated(ig_point%x)) allocate(ig_point%x(NDIM))
      ig_point%x = x
    end if

    if(present(u)) then
      if (.not. allocated(ig_point%u)) allocate(ig_point%u(NDIM))
      ig_point%u = u
    end if

    if(present(traction)) then
      if (.not. allocated(ig_point%traction)) allocate(ig_point%traction(NST))
      ig_point%traction = traction
    end if

    if(present(separation)) then
      if (.not. allocated(ig_point%separation)) allocate(ig_point%separation(NST))
      ig_point%separation = separation
    end if

    if(present(converged_sdv))  ig_point%converged_sdv = converged_sdv
    if(present(iterating_sdv))  ig_point%iterating_sdv = iterating_sdv

  end subroutine update_cohesive_ig_point



  pure subroutine extract_cohesive_ig_point (ig_point, x, u, &
  & traction, separation, converged_sdv, iterating_sdv)

    type(cohesive_ig_point),      intent(in)  :: ig_point
    real(DP),           optional, intent(out) :: x(NDIM), u(NDIM)
    real(DP),           optional, intent(out) :: traction(NST), separation(NST)
    type(cohesive_sdv), optional, intent(out) :: converged_sdv
    type(cohesive_sdv), optional, intent(out) :: iterating_sdv

    ! initialize intent(out) variables
    x = ZERO
    u = ZERO
    traction   = ZERO
    separation = ZERO
    ! derived types are automatically initialized upon declaration

    if(present(x)) then
        if(allocated(ig_point%x)) x = ig_point%x
      end if

      if(present(u)) then
        if(allocated(ig_point%u)) u = ig_point%u
      end if

      if(present(traction)) then
        if(allocated(ig_point%traction)) traction = ig_point%traction
      end if

      if(present(separation)) then
        if(allocated(ig_point%separation)) separation = ig_point%separation
      end if

    if(present(converged_sdv))   converged_sdv = ig_point%converged_sdv
    if(present(iterating_sdv))   iterating_sdv = ig_point%iterating_sdv

  end subroutine extract_cohesive_ig_point



  subroutine display_cohesive_ig_point (this)
  ! Purpose:
  ! to display this cohesive_ig_point's components on cmd window
  ! this is useful for debugging

    type(cohesive_ig_point), intent(in) :: this

    ! local variable to set the output format
    character(len=20) :: display_fmt
    integer           :: i

    ! initialize local variable
    display_fmt = ''
    i = 0

    ! set display format for real
    ! ES for real (scientific notation)
    ! 10 is width, 3 is no. of digits aft decimal point
    ! note that for scientific real, ESw.d, w>=d+7
    display_fmt = '(ES10.3)'

    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the components of the inquired &
                       &cohesive_ig_point object :'
    write(*,'(1X, A)') ''

    if (allocated(this%x)) then
      write(*,'(1X, A)') '- x of this cohesive_ig_point is: '
      do i = 1, NDIM
        write(*,display_fmt,advance="no") this%x(i)
      end do
      write(*,'(1X, A)') ''
      write(*,'(1X, A)') ''
    end if

    if (allocated(this%u)) then
      write(*,'(1X, A)') '- u of this cohesive_ig_point is: '
      do i = 1, NDIM
        write(*,display_fmt,advance="no") this%u(i)
      end do
      write(*,'(1X, A)') ''
      write(*,'(1X, A)') ''
    end if

    if (allocated(this%traction)) then
      write(*,'(1X, A)') '- traction at this cohesive_ig_point is: '
      do i = 1, NST
        write(*,display_fmt,advance="no") this%traction(i)
      end do
      write(*,'(1X, A)') ''
      write(*,'(1X, A)') ''
    end if

    if (allocated(this%separation)) then
      write(*,'(1X, A)') '- separation at this cohesive_ig_point is: '
      do i = 1, NST
        write(*,display_fmt,advance="no") this%separation(i)
      end do
      write(*,'(1X, A)') ''
      write(*,'(1X, A)') ''
    end if

    write(*,'(1X, A)') '- converged cohesive_sdv of this cohesive_ig_point is:'
    call display_cohesive_sdv(this%converged_sdv)
    write(*,'(1X, A)') ''

    write(*,'(1X, A)') '- iterating cohesive_sdv of this cohesive_ig_point is: '
    call display_cohesive_sdv(this%iterating_sdv)
    write(*,'(1X, A)') ''

  end subroutine display_cohesive_ig_point





end module cohesive_material_module
