module delam6_element_module
!
!  Purpose:
!    define a 3D 6-Node cohesive element object for delamination
!    with the associated procedures to empty, set, integrate and extract
!    its components
!
!  topological definition of this wedge element, local nodal indices:
!
!                     6
!                    /|\
!                   / | \
!                  /  |  \
!                 /   |   \
!                /    |    \
!               /     |     \
!             4/______|______\5
!             |       |       |
!             |      /3\      |
!             |     /   \     |
!             |    /     \    |
!             |   /       \   |
!             |  /         \  |
!             | /           \ |
!             |/_____________\|
!             1               2
!
!  bottom surface nodes (counter-clock vise): 1, 2, 3
!  top    surface nodes (counter-clock vise): 4, 5, 6
!
!  coord. system: x-y-z, z is normal dir.
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    29/06/15  B. Y. Chen            Original code
!

use parameter_module, only : NST => NST_COHESIVE, NDIM, DP, ZERO, HALF, ONE, &
                      & ONE_SIXTH, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,    &
                      & SMALLNUM, HALFCIRC, PI, INTACT
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
integer, parameter :: NNODE=6, NIGPOINT=3, NDOF=NDIM*NNODE


type, public :: delam6_element
  private
  ! list of type components:
  ! fstat         : element failure status
  ! connec        : indices of element nodes in the global node array
  ! ig_points     : integration points of this element
  ! ig_angles     : delamination longitudial angles at ig points
  ! local_clock   : locally-saved program clock
  ! traction      : traction on the interface, for output
  ! separation    : separation on the interface, for output
  ! dm            : matrix degradation factor for output
  integer  :: fstat         = 0
  integer  :: connec(NNODE) = 0
  type(program_clock)     :: local_clock
  type(cohesive_ig_point) :: ig_points(NIGPOINT)
  real(DP)                :: ig_angles(NIGPOINT) = ZERO
  real(DP) :: traction(NST)   = ZERO
  real(DP) :: separation(NST) = ZERO
  real(DP) :: dm              = ZERO
end type




interface empty
  module procedure empty_delam6_element
end interface

interface set
  module procedure set_delam6_element
end interface

interface integrate
  module procedure integrate_delam6_element
end interface

interface extract
  module procedure extract_delam6_element
end interface




public :: empty, set, integrate, extract




contains




pure subroutine empty_delam6_element (elem)
! Purpose:
! this subroutine is used to format the element for use
! it is used in the initialize_lib_elem procedure in the lib_elem module

  type(delam6_element), intent(inout) :: elem

  ! local variable, derived type var. is initialized upon declaration
  type(delam6_element) :: elem_lcl

  ! reset elem to the initial state
  elem = elem_lcl

end subroutine empty_delam6_element



pure subroutine set_delam6_element (elem, connec, istat, emsg)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
! note that only some of the components need to be set during preproc,
! namely connec, ID_matlist

  type(delam6_element),   intent(inout)   :: elem
  integer,                intent(in)      :: connec(NNODE)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg

  istat = STAT_SUCCESS
  emsg  = ''

  ! check validity of inputs
  if ( any(connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'connec node indices must be >=1, set, &
    &delam6_element_module'
    return
  end if

  elem%connec    = connec

end subroutine set_delam6_element



pure subroutine extract_delam6_element (elem, fstat, connec, ig_points, &
& ig_angles, traction, separation, dm)
! Purpose:
! to extract the components of this element
! note that the dummy args connec and ig_points are allocatable arrays
! because their sizes vary with different element types

  type(delam6_element),                           intent(in)  :: elem
  integer,                              optional, intent(out) :: fstat
  integer,                 allocatable, optional, intent(out) :: connec(:)
  type(cohesive_ig_point), allocatable, optional, intent(out) :: ig_points(:)
  real(DP),                allocatable, optional, intent(out) :: ig_angles(:)
  real(DP),                             optional, intent(out) :: traction(NST)
  real(DP),                             optional, intent(out) :: separation(NST)
  real(DP),                             optional, intent(out) :: dm

  if (present(fstat))       fstat = elem%fstat

  if (present(connec)) then
    allocate(connec(NNODE))
    connec = elem%connec
  end if

  if (present(ig_points)) then
    allocate(ig_points(NIGPOINT))
    ig_points = elem%ig_points
  end if

  if (present(ig_angles)) then
    allocate(ig_angles(NIGPOINT))
    ig_angles = elem%ig_angles
  end if

  if (present(traction))    traction    = elem%traction

  if (present(separation))  separation  = elem%separation

  if (present(dm))          dm          = elem%dm

end subroutine extract_delam6_element



pure subroutine integrate_delam6_element (elem, nodes, material, theta1, theta2,&
& K_matrix, F_vector, istat, emsg, nofailure)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element
! note that K and F are allocatable dummy args as their sizes vary with
! different element types

! list of used modules:
! xnode_module                  : xnode derived type and its assoc. procedures
! global_toolkit_module         : global tools for element integration
use xnode_module,                only : xnode, extract
use global_toolkit_module,       only : cross_product3d, normalize_vect

  ! most args are self-explanatory
  ! fdir: fibre direction
  type(delam6_element),     intent(inout) :: elem
  type(xnode),              intent(in)    :: nodes(NNODE)
  type(cohesive_material),  intent(in)    :: material
  real(DP),                 intent(in)    :: theta1, theta2
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
  real(DP)            :: ig_angles(NIGPOINT)
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
  ! - is_zero_vect    : true of any of the above vectors has zero length
  real(DP)            :: normal(NDIM), tangent1(NDIM), tangent2(NDIM)
  logical             :: is_zero_vect
  !** element determinant
  ! - det             : determinant of jacobian
  real(DP)            :: det

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

  ! - c(s)theta1(2)   : cos and sin of theta1 and theta2
  ! - deltaL1(2)      : longitudinal separation, along theta1(2)
  ! - Qmatrix         : rotation matrix from global to local coordinates
  real(DP)            :: ctheta1, stheta1, ctheta2, stheta2, deltaL1, deltaL2
  real(DP)            :: Qmatrix(NDIM,NDIM)

  !** variables needed for stiffness matrix derivation:
  ! - fn              : shape functions
  ! - Nmatrix         : obtained from fn to compute disp. jump vector ujump
  !                     across the interface: {u}_jump = [N]*{u}
  ! - ujump, delta    : {delta}=[Q]*{u}_jump, {delta} is the separation vector.
  ! - dee             : material stiffness matrix [D]
  ! - tau             : {tau}=[D]*{delta}, traction on the interface
  ! - QN              : [Q]*[N]
  ! - QNtDQN          : [QN]'*[D]*[QN]
  ! - QNtTau          : [QN]'*{tau}
  real(DP)            :: fn(NNODE)
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
  !** normal and tangents
  normal          = ZERO
  tangent1        = ZERO
  tangent2        = ZERO
  is_zero_vect    = .false.
  !** element determinant
  det             = ZERO
  !** analysis logical control variables:
  last_converged  = .false.
  nofail          = .false.
  !** integration point variables:
  ig_xi           = ZERO
  ig_weight       = ZERO
  !** sin, cos of theta1 and theta2 and separ. along these dir.
  ctheta1         = ZERO
  stheta1         = ZERO
  ctheta2         = ZERO
  stheta2         = ZERO
  deltaL1         = ZERO
  deltaL2         = ZERO
  !** rotational matrix
  Qmatrix         = ZERO
  !** variables needed for stiffness matrix derivation:
  fn              = ZERO
  Nmatrix         = ZERO
  ujump           = ZERO
  delta           = ZERO
  dee             = ZERO
  tau             = ZERO
  ig_fstat        = 0
  ig_dm           = ZERO
  QN              = ZERO
  QNtDQN          = ZERO
  QNtTau          = ZERO

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
  ig_angles     = elem%ig_angles    ! inout

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
      emsg  = 'x not allocated for node, delam6_element_module'
      exit
    end if
    ! assign nodal displacement values (uj) to u vector
    if(allocated(uj)) then
      u((j-1)*NDIM+1:j*NDIM)=uj(1:NDIM)
      deallocate(uj)
    else
      istat = STAT_FAILURE
      emsg  = 'u not allocated for node, delam6_element_module'
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

  ! check if all midcoords are on the x-y plane
  normal = [ZERO, ZERO, ONE]
  ! compute tangent1 of the interface: node 2 coords - node 1 coords
  tangent1(:)=midcoords(:,2)-midcoords(:,1)
  ! compute tangent2 of the interface: node 3 coords - node 1 coords
  tangent2(:)=midcoords(:,3)-midcoords(:,1)
  ! check if tangent1 is on x-y plane
  if (dot_product(tangent1,normal) > SMALLNUM) then
    istat = STAT_FAILURE
    emsg  = 'edge 1-2 not on x-y plane, delam6 element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  ! check if tangent2 is on x-y plane
  if (dot_product(tangent2,normal) > SMALLNUM) then
    istat = STAT_FAILURE
    emsg  = 'edge 1-3 not on x-y plane, delam6 element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

  ! compute determinant of jacobian
  ! note that determinant of linear tri elem is constant and is equal to
  ! Atri / Atri_reference = Atri / 0.5 = 2 * Atri
  ! so the length of tangent1 cross tangent2 = determinant of this elem
  normal = cross_product3d(tangent1,tangent2)
  call normalize_vect(normal, is_zero_vect, det)
  if (is_zero_vect) then
    istat = STAT_FAILURE
    emsg  = 'element area is ZERO, delam6 element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

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

      !- get shape function matrix
      call init_shape (fn, ig_xi(:,kig))

      ! Nmatrix: ujump (at each ig point) = Nmatrix * u
      do i = 1, NDIM
        do j = 1, NNODE/2
          Nmatrix( i, i + (j-1)         * NDIM ) = - fn(j)
          Nmatrix( i, i + (j-1+NNODE/2) * NDIM ) =   fn(j)
        end do
      end do

      ! calculate ujump: disp. jump of the two crack surface, in global coords
      ujump = matmul(Nmatrix,u)

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

      ! calculate Qmatrix based on status of ig point
      ! - if ig point already starting to fail, then use the stored theta of this
      !   ig point.
      ! - else if ig point is still intact, then find the preferred theta from
      !   input theta1 and theta2 by comparing the projection of ujump on the
      !   direction of theta1/2. The one with the larger projection is preferred
      call extract (ig_sdv_iter, fstat=ig_fstat)
      if (ig_fstat > INTACT) then
        ctheta1 = cos(ig_angles(kig)/HALFCIRC*PI)
        stheta1 = sin(ig_angles(kig)/HALFCIRC*PI)
        Qmatrix(1,:)=[    ZERO,    ZERO, ONE ]
        Qmatrix(2,:)=[ ctheta1, stheta1, ZERO]
        Qmatrix(3,:)=[-stheta1, ctheta1, ZERO]
      else
        ! find the preferred delam longitudinal direction (theta1 or theta2)
        ctheta1 = cos(theta1/HALFCIRC*PI)
        stheta1 = sin(theta1/HALFCIRC*PI)
        ctheta2 = cos(theta2/HALFCIRC*PI)
        stheta2 = sin(theta2/HALFCIRC*PI)
        deltaL1 = dot_product(ujump(1:2),[ctheta1,stheta1])
        deltaL2 = dot_product(ujump(1:2),[ctheta2,stheta2])
        if(abs(deltaL1) > abs(deltaL2)) then
          ig_angles(kig) = theta1
          Qmatrix(1,:)=[    ZERO,    ZERO, ONE ]
          Qmatrix(2,:)=[ ctheta1, stheta1, ZERO]
          Qmatrix(3,:)=[-stheta1, ctheta1, ZERO]
        else
          ig_angles(kig) = theta2
          Qmatrix(1,:)=[    ZERO,    ZERO, ONE ]
          Qmatrix(2,:)=[ ctheta2, stheta2, ZERO]
          Qmatrix(3,:)=[-stheta2, ctheta2, ZERO]
        end if
      end if

      ! calculate separation delta in local coords: delta = Qmatrix * ujump
      delta = matmul(Qmatrix,ujump)

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
      Nmatrix = ZERO
      ujump   = ZERO
      delta   = ZERO
      dee     = ZERO
      tau     = ZERO
      QN      = ZERO
      QNtDQN  = ZERO
      QNtTau  = ZERO
      ig_fstat= 0
      ig_dm   = ZERO

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
  elem%ig_angles   = ig_angles

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

end subroutine integrate_delam6_element






! the rest are private subroutines






pure subroutine init_ig_point (xi, wt, gauss)
! used parameters:
! DP, HALF, ZERO, ONE, ONE_SIXTH

  real(DP),          intent(inout) :: xi(NDIM-1,NIGPOINT), wt(NIGPOINT)
  logical, optional, intent(in)    :: gauss

  logical :: isgauss

  isgauss = .false.

  if(present(gauss)) isgauss = gauss

  if(isgauss) then
  ! gauss integration points
    xi(1,1) = HALF
    xi(2,1) = HALF
    xi(1,2) = HALF
    xi(2,2) = ZERO
    xi(1,3) = ZERO
    xi(2,3) = HALF
  else
  ! newton-cotes integration points
    xi(1,1) = ZERO
    xi(2,1) = ZERO
    xi(1,2) = ONE
    xi(2,2) = ZERO
    xi(1,3) = ZERO
    xi(2,3) = ONE
  end if
  wt = ONE_SIXTH

end subroutine init_ig_point



pure subroutine init_shape (f, igxi)
! used parameters:
! DP, ZERO, ONE

  real(DP), intent(inout) :: f(NNODE)
  real(DP), intent(in)    :: igxi(NDIM-1)

  real(DP) :: xi, eta ! local variables
  xi  = ZERO
  eta = ZERO

  xi  = igxi(1)
  eta = igxi(2)
  f(1) = ONE-xi-eta
  f(2) = xi
  f(3) = eta
  f(4) = f(1)
  f(5) = f(2)
  f(6) = f(3)

end subroutine init_shape




end module delam6_element_module
