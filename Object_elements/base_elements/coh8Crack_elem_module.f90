module coh8Crack_elem_module
!
!  Purpose:
!    define a cohCrack element object for cohesive matrix crack in ply
!    with the associated procedures to empty, set, integrate and extract
!    its components
!
!  topological definition of this cohCrack element, local nodal indices:
!
!  8___________________7
!  |\                  |\
!  | \                 | \
!  |  \                |  \
!  |   \               |   \
!  |    \              |    \
!  |     \5____________|_____\6
!  |______|____________|      |
! 4\      |           3\      |
!   \     |             \     |
!    \    |              \    |
!     \   |               \   |
!      \  |                \  |
!       \ |                 \ |
!        \|__________________\|
!         1                   2
!
!  bottom surface nodes (counter-clock vise): 1, 2, 3, 4
!  top    surface nodes (counter-clock vise): 5, 6, 7, 8
!
!  default fibre dir: 1-2
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    28/06/15  B. Y. Chen            Original code
!

use parameter_module, only : NST => NST_COHESIVE, NDIM, DP, ZERO, ONE, HALF, &
                      & ONE_ROOT3, QUARTER, &
                      & MSGLENGTH, STAT_SUCCESS, STAT_FAILURE
! list of external modules used in type definition and other procedures:
! global clock module      : needed in element definition, extract and integrate
! cohesive material module : needed in element definition, extract and integrate
use global_clock_module,      only : program_clock, GLOBAL_CLOCK, clock_in_sync
use cohesive_material_module, only : cohesive_ig_point, cohesive_material, &
                              & cohesive_sdv, ddsdde, update, extract

implicit none
private

! list of private parameters in this module:
! NNODE         : no. of nodes in this element
! NIGPOINT      : no. of integration points in this element
! NDOF          : no. of degree-of-freedom  in this element
!
integer, parameter :: NNODE=8, NIGPOINT=4, NDOF=NDIM*NNODE


type, public :: coh8Crack_elem
  private
  ! list of type components:
  ! fstat         : element failure status
  ! ig_points     : integration points of this element
  ! local_clock   : locally-saved program clock
  ! traction      : traction on the interface, for output
  ! separation    : separation on the interface, for output
  ! dm            : matrix degradation factor for output
  integer  :: fstat         = 0
  type(program_clock)     :: local_clock
  type(cohesive_ig_point) :: ig_points(NIGPOINT)
  real(DP) :: traction(NST)   = ZERO
  real(DP) :: separation(NST) = ZERO
  real(DP) :: dm              = ZERO
end type


interface integrate
  module procedure integrate_coh8Crack_elem
end interface

interface extract
  module procedure extract_coh8Crack_elem
end interface




public :: integrate, extract




contains




pure subroutine extract_coh8Crack_elem (elem, fstat, ig_points, &
& traction, separation, dm)
! Purpose:
! to extract the components of this element
! note that the dummy args connec and ig_points are allocatable arrays
! because their sizes vary with different element types

  type(coh8Crack_elem),                           intent(in)  :: elem
  integer,                              optional, intent(out) :: fstat
  type(cohesive_ig_point), allocatable, optional, intent(out) :: ig_points(:)
  real(DP),                             optional, intent(out) :: traction(NST)
  real(DP),                             optional, intent(out) :: separation(NST)
  real(DP),                             optional, intent(out) :: dm

  if (present(fstat))       fstat = elem%fstat

  if (present(ig_points)) then
    allocate(ig_points(NIGPOINT))
    ig_points = elem%ig_points
  end if

  if (present(traction))    traction    = elem%traction

  if (present(separation))  separation  = elem%separation

  if (present(dm))          dm          = elem%dm

end subroutine extract_coh8Crack_elem



pure subroutine integrate_coh8Crack_elem (elem, nodes, material, &
& K_matrix, F_vector, istat, emsg, nofailure)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element
! note that K and F are allocatable dummy args as their sizes vary with
! different element types

! list of used modules:
! fnode_module                  : fnode derived type and its assoc. procedures
! global_toolkit_module         : global tools for element integration
use fnode_module,                only : fnode, extract
use global_toolkit_module,       only : cross_product3d, normalize_vect, &
                                 & determinant2d
  ! most args are self-explanatory
  ! fdir: fibre direction
  type(coh8Crack_elem),     intent(inout) :: elem
  type(fnode),              intent(in)    :: nodes(NNODE)
  type(cohesive_material),  intent(in)    :: material
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,       optional,  intent(in)    :: nofailure

  ! local copies of intent(inout) dummy arg./its components:
  ! - elfstat         : elem's fstat
  ! - local_clock     : elem's local_clock
  ! - ig_points       : elem's ig_points
  ! - eltraction      : elem's traction
  ! - elseparation    : elem's separation
  ! - eldm            : elem's dm
  integer             :: elfstat
  type(program_clock) :: local_clock
  type(cohesive_ig_point) :: ig_points(NIGPOINT)
  real(DP)            :: eltraction(NST)
  real(DP)            :: elseparation(NST)
  real(DP)            :: eldm

  ! the rest are all local variables

  !** nodal variables:
  ! - xj, uj          : nodal x and u extracted from nodes array
  ! - coords          : element nodal coordinates matrix
  ! - u               : element nodal displacemet vector
  ! - midcoords       : coordinates of the mid-plane
  real(DP), allocatable :: xj(:), uj(:)
  real(DP)            :: coords(NDIM,NNODE), u(NDOF)
  real(DP)            :: midcoords(NDIM,NNODE/2)

  !** local coords and rotational matrix
  ! - normal, tangent1/2 : normal and tangent vectors of the interface,
  !                     obtained from coords matrix
  ! - is_zero_vect    : true of any of the above vectors has ZERO length
  ! - Qmatrix         : rotation matrix from global to local coordinates
  !                     abtained from normal & tangent vectors
  real(DP)            :: normal(NDIM), tangent1(NDIM), tangent2(NDIM)
  logical             :: is_zero_vect
  real(DP)            :: Qmatrix(NDIM,NDIM)

  !** analysis logical control variables:
  ! - last_converged  : true if last iteration has converged
  ! - nofail          : true if no failure is allowed
  logical             :: last_converged, nofail

  !** integration point variables:
  ! - ig_xi           : integration point natural coordinates
  ! - ig_weight       : integration point weights
  ! - ig_sdv_conv/iter: integration point converged and iterating sdvs
  real(DP)            :: ig_xi(NDIM-1, NIGPOINT), ig_weight(NIGPOINT)
  type(cohesive_sdv)  :: ig_sdv_conv, ig_sdv_iter

  !** variables needed for stiffness matrix derivation:
  ! - fn, dn          : shape functions and their derivatives
  ! - jac, det        : element jacobian and its determinant
  ! - Nmatrix         : obtained from fn to compute disp. jump vector ujump
  !                     across the interface: {u}_jump = [N]*{u}
  ! - ujump, delta    : {delta}=[Q]*{u}_jump, {delta} is the separation vector.
  ! - dee             : material stiffness matrix [D]
  ! - tau             : {tau}=[D]*{delta}, traction on the interface
  ! - QN              : [Q]*[N]
  ! - QNtDQN          : [QN]'*[D]*[QN]
  ! - QNtTau          : [QN]'*{tau}
  real(DP)            :: fn(NNODE), dn(NNODE,NDIM-1)
  real(DP)            :: jac(NDIM-1,NDIM-1), det
  real(DP)            :: Nmatrix(NDIM,NDOF)
  real(DP)            :: ujump(NDIM), delta(NDIM)
  real(DP)            :: dee(NDIM,NDIM)
  real(DP)            :: tau(NDIM)
  real(DP)            :: QN(NDIM,NDOF)
  real(DP)            :: QNtDQN(NDOF,NDOF)
  real(DP)            :: QNtTau(NDOF)
  integer             :: ig_fstat
  real(DP)            :: ig_dm

  !** integer counter variables
  integer             :: i, j, k, kig


  ! initialize intent(out) and local variables
  ! except derived types (auto initialized upon declaration) and
  ! allocatable local arrays

  !** intent(out) variables:
  allocate(K_matrix(NDOF,NDOF), F_vector(NDOF))
  K_matrix        = ZERO
  F_vector        = ZERO
  istat           = STAT_SUCCESS
  emsg            = ''
  !** local copies of intent(inout) dummy args./its components:
  elfstat         = 0
  eltraction      = ZERO
  elseparation    = ZERO
  eldm            = ZERO
  !** nodal variables:
  coords          = ZERO
  u               = ZERO
  midcoords       = ZERO
  !** local coords and rotational matrix
  normal          = ZERO
  tangent1        = ZERO
  tangent2        = ZERO
  is_zero_vect    = .false.
  Qmatrix         = ZERO
  !** analysis logical control variables:
  last_converged  = .false.
  nofail          = .false.
  !** integration point variables:
  ig_xi           = ZERO
  ig_weight       = ZERO
  !** variables needed for stiffness matrix derivation:
  fn              = ZERO
  dn              = ZERO
  jac             = ZERO
  det             = ZERO
  Nmatrix         = ZERO
  ujump           = ZERO
  delta           = ZERO
  dee             = ZERO
  tau             = ZERO
  QN              = ZERO
  QNtDQN          = ZERO
  QNtTau          = ZERO
  ig_fstat        = 0
  ig_dm           = ZERO
  !** integer counter variables:
  i=0; j=0; k=0; kig=0


  ! check validity of input/imported variables

  ! assign values to local variables

  !** local copies of intent(inout) dummy arg.:
  ! elem's fstat, traction, separation and dm components are NOT needed
  ! for calculations in this subroutine. here they are only for output,
  ! and they need to be assgined new values during each iteration.
  ! they are therefore NOT passed to their local copies elfstat, eltraction,
  ! elseparation and dm.
  ! elfstat     : no extraction     ! out
  ! eltraction  : no extraction     ! out
  ! elseparation: no extraction     ! out
  ! eldm        : no extraction     ! out
  local_clock   = elem%local_clock  ! inout
  ig_points     = elem%ig_points    ! inout

  !** nodal variables:
  ! extract nodal components and assign to respective local arrays:
  ! nodal x -> coords
  ! nodal u -> u
  do j=1, NNODE
    ! extract x and u from nodes
    call extract(nodes(j), x=xj, u=uj)
    ! assign nodal coordinate values (xj) to coords matrix
    if(allocated(xj)) then
      coords(:,j)=xj(:)
      deallocate(xj)
    else
      istat = STAT_FAILURE
      emsg  = 'x not allocated for node, coh8Crack_elem_module'
      exit
    end if
    ! assign nodal displacement values (uj) to u vector
    if(allocated(uj)) then
      u((j-1)*NDIM+1:j*NDIM)=uj(1:NDIM)
      deallocate(uj)
    else
      istat = STAT_FAILURE
      emsg  = 'u not allocated for node, coh8Crack_elem_module'
      exit
    end if
  end do
  ! if there's any error encountered in the extraction process
  ! clean up and exit the program
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  ! calculate mid-plane coordinates from coords
  do j=1, NNODE/2
      midcoords(:,j) = HALF * (coords(:,j)+coords(:,j+NNODE/2))
  end do

  !** local coords and rotational matrix, and element determinant:
  ! compute tangent1 of the interface: node 2 coords - node 1 coords
  tangent1(:)=midcoords(:,2)-midcoords(:,1)
  ! compute tangent2 of the interface: node 4 coords - node 1 coords
  tangent2(:)=midcoords(:,4)-midcoords(:,1)
  ! compute normal vector of the interface,
  normal = cross_product3d(tangent1,tangent2)
  ! - re-evaluate tangent2 so that it is perpendicular to both
  ! tangent1 and normal
  tangent2 = cross_product3d(normal,tangent1)
  ! - normalize these vectors, exit w error if ZERO vector is encountered
  call normalize_vect(tangent1, is_zero_vect)
  if (is_zero_vect) then
    istat = STAT_FAILURE
    emsg  = 'element tangent1 is ZERO, cohCrack element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  call normalize_vect(tangent2, is_zero_vect)
  if (is_zero_vect) then
    istat = STAT_FAILURE
    emsg  = 'element tangent2 is ZERO, cohCrack element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  call normalize_vect(normal, is_zero_vect)
  if (is_zero_vect) then
    istat = STAT_FAILURE
    emsg  = 'element area is ZERO, cohCrack element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  ! - compute Q matrix
  Qmatrix(1,:) =  normal(:)   ! normal
  Qmatrix(2,:) = -tangent2(:) ! transverse (note minus sign)
  Qmatrix(3,:) =  tangent1(:) ! longitudinal
  ! - rotate midcoords to planar coord sys.
  midcoords = matmul(Qmatrix,midcoords)

  !** analysis logical control variables:
  ! check if last iteration has converged by checking if the global clock has
  ! advanced; if so, last iteration is converged and sync the local clock
  if (.not. clock_in_sync(GLOBAL_CLOCK, local_clock)) then
    last_converged = .true.
    local_clock    = GLOBAL_CLOCK
  end if
  ! nofail:
  if (present(nofailure)) nofail = nofailure

  !** integration point variables:
  ! calculate ig point xi and weight
  call init_ig_point (ig_xi, ig_weight)

  !** other local variables are assigned in the following loop of all ig points



  !**** MAIN CALCULATIONS ****


  ! loop over all ig points,
  ! calculate separation, traction, stiffness matrix, sdv etc. at each ig point

  loop_igpoint: do kig=1, NIGPOINT

      ! get shape function matrix
      call init_shape (fn, dn, ig_xi(:,kig))

      ! calculate jacobian of element at this ig point using planar coordinates
      ! 1st row of midcoords are planar normal-dir. coords of the elem nodes
      ! only the planar tangent dir. coords are used to calculate jacobian w.r.t
      ! the reference rectangle of quadrilateral
      ! only the shape funcs of the first half of nodes are needed
      jac = matmul(midcoords(2:3,:),dn(1:NNODE/2,:))
      det = determinant2d(jac)
      if (det <= ZERO) then
        istat = STAT_FAILURE
        emsg  = 'det of jacob is zero or negative, cohCrack element module'
        return
      end if

      ! Nmatrix: ujump (at each ig point) = Nmatrix * u
      do i = 1, NDIM
        do j = 1, NNODE/2
          Nmatrix( i, i + (j-1)         * NDIM ) = - fn(j)
          Nmatrix( i, i + (j-1+NNODE/2) * NDIM ) =   fn(j)
        end do
      end do

      ! calculate ujump: disp. jump of the two crack surface, in global coords
      ujump = matmul(Nmatrix,u)

      ! calculate separation delta in local coords: delta = Qmatrix * ujump
      delta = matmul(Qmatrix,ujump)

      ! extract sdvs from integration points, ig_sdv_conv/iter
      call extract(ig_points(kig), converged_sdv=ig_sdv_conv, &
      & iterating_sdv=ig_sdv_iter)

      ! update converged sdv with iterating sdv when last iteration is converged
      ! and revert iterating sdv back to the last converged sdv if otherwise
      if(last_converged) then
        ig_sdv_conv = ig_sdv_iter
      else
        ig_sdv_iter = ig_sdv_conv
      end if

      ! use material properties, sdv_iter, and separation
      ! to calculate D and traction, and update sdv_iter
      if(nofail) then
      ! no failure is allowed, use ddsdde_cohesve_intact
        call ddsdde (material, dee=dee, traction=tau, separation=delta)
      else
      ! failure is allowed, use ddsdde_cohesive
        call ddsdde (material, dee=dee, traction=tau, sdv=ig_sdv_iter, &
        & separation=delta, istat=istat, emsg=emsg)
        if (istat == STAT_FAILURE) exit loop_igpoint
      end if

      ! add this ig point contributions to K and F
      QN     =   matmul(Qmatrix, Nmatrix)
      QNtDQN =   matmul(transpose(QN), matmul(dee, QN))
      QNtTau =   matmul(transpose(QN), tau)

      do j=1, NDOF
        do k=1, NDOF
          K_matrix(j,k) = K_matrix(j,k) + QNtDQN(j,k) * det * ig_weight(kig)
        end do
        F_vector(j) = F_vector(j) + QNtTau(j) * ig_weight(kig) * det
      end do

      ! update ig point components, only sdv are stored
      call update(ig_points(kig), converged_sdv=ig_sdv_conv, &
      & iterating_sdv=ig_sdv_iter)

      ! update elem fstat to be the max current ig point fstat
      call extract (ig_sdv_iter, fstat=ig_fstat)
      elfstat  = max(elfstat, ig_fstat)

      ! update elem dm to be the max current ig point dm
      call extract (ig_sdv_iter, dm=ig_dm)
      eldm     = max(eldm,    ig_dm)

      ! update elem stress & strain (avg of ig point stress & strains)
      eltraction   = eltraction   + tau / real(NIGPOINT, DP)
      elseparation = elseparation + delta / real(NIGPOINT, DP)

      ! empty relevant arrays for reuse
      fn      = ZERO
      dn      = ZERO
      jac     = ZERO
      det     = ZERO
      Nmatrix = ZERO
      ujump   = ZERO
      delta   = ZERO
      dee     = ZERO
      tau     = ZERO
      QN      = ZERO
      QNtDQN  = ZERO
      QNtTau  = ZERO

  end do loop_igpoint !-looped over all int points. ig=NIGPOINT

  ! check to see if the loop is exited upon error
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

  !**** END MAIN CALCULATIONS ****


  ! check to see if intent(inout) derived type dummy arg's components are
  ! unintentionally modified

  ! if no unintentinal modification, proceed with final calculations
  ! and updates and exit the program

  ! update intent(inout) dummy arg./its components
  elem%fstat       = elfstat
  elem%dm          = eldm
  elem%traction    = eltraction
  elem%separation  = elseparation
  elem%local_clock = local_clock
  elem%ig_points   = ig_points

  ! deallocate local dynamic arrays
  if(allocated(xj)) deallocate(xj)
  if(allocated(uj)) deallocate(uj)



  contains
  ! internal procedures

    pure subroutine clean_up (K_matrix, F_vector, uj, xj)

      real(DP),              intent(inout) :: K_matrix(:,:), F_vector(:)
      real(DP), allocatable, intent(inout) :: uj(:), xj(:)

      ! ZERO intent(out) variables (do not deallocate)
      K_matrix = ZERO
      F_vector = ZERO
      ! deallocate local alloc. variables
      if (allocated(uj)) deallocate(uj)
      if (allocated(xj)) deallocate(xj)

    end subroutine clean_up

end subroutine integrate_coh8Crack_elem






! the rest are private subroutines






pure subroutine init_ig_point (xi, wt, gauss)
! used parameters:
! DP, ONE, ONE_ROOT3

  real(DP),          intent(inout) :: xi(NDIM-1,NIGPOINT), wt(NIGPOINT)
  logical, optional, intent(in)    :: gauss

  logical :: isgauss

  isgauss = .false.

  if(present(gauss)) isgauss = gauss

  if(isgauss) then
  ! gauss integration points
    xi(1,1) = -ONE_ROOT3
    xi(2,1) = -ONE_ROOT3
    xi(1,2) =  ONE_ROOT3
    xi(2,2) = -ONE_ROOT3
    xi(1,3) = -ONE_ROOT3
    xi(2,3) =  ONE_ROOT3
    xi(1,4) =  ONE_ROOT3
    xi(2,4) =  ONE_ROOT3
  else
  ! newton-cotes integration points
    xi(1,1) = -ONE
    xi(2,1) = -ONE
    xi(1,2) =  ONE
    xi(2,2) = -ONE
    xi(1,3) = -ONE
    xi(2,3) =  ONE
    xi(1,4) =  ONE
    xi(2,4) =  ONE
  end if
  wt = ONE

end subroutine init_ig_point



pure subroutine init_shape (f, df, igxi)
! used parameters:
! DP, ZERO, ONE, QUARTER

  real(DP), intent(inout) :: f(NNODE), df(NNODE,NDIM-1)
  real(DP), intent(in)    :: igxi(NDIM-1)

  ! local variables
  real(DP) :: xi, eta

  xi  = ZERO
  eta = ZERO

  xi  = igxi(1)
  eta = igxi(2)

  f(1)    =  QUARTER*(ONE-xi)*(ONE-eta)
  f(2)    =  QUARTER*(ONE+xi)*(ONE-eta)
  f(3)    =  QUARTER*(ONE+xi)*(ONE+eta)
  f(4)    =  QUARTER*(ONE-xi)*(ONE+eta)
  df(1,1) = -QUARTER*(ONE-eta)
  df(2,1) =  QUARTER*(ONE-eta)
  df(3,1) =  QUARTER*(ONE+eta)
  df(4,1) = -QUARTER*(ONE+eta)
  df(1,2) = -QUARTER*(ONE-xi)
  df(2,2) = -QUARTER*(ONE+xi)
  df(3,2) =  QUARTER*(ONE+xi)
  df(4,2) =  QUARTER*(ONE-xi)

  f(1+4)    =  QUARTER*(ONE-xi)*(ONE-eta)
  f(2+4)    =  QUARTER*(ONE+xi)*(ONE-eta)
  f(3+4)    =  QUARTER*(ONE+xi)*(ONE+eta)
  f(4+4)    =  QUARTER*(ONE-xi)*(ONE+eta)
  df(1+4,1) = -QUARTER*(ONE-eta)
  df(2+4,1) =  QUARTER*(ONE-eta)
  df(3+4,1) =  QUARTER*(ONE+eta)
  df(4+4,1) = -QUARTER*(ONE+eta)
  df(1+4,2) = -QUARTER*(ONE-xi)
  df(2+4,2) = -QUARTER*(ONE+xi)
  df(3+4,2) =  QUARTER*(ONE+xi)
  df(4+4,2) =  QUARTER*(ONE-xi)

end subroutine init_shape



end module coh8Crack_elem_module
