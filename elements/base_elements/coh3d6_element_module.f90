module coh3d6_element_module
!
!  Purpose:
!    define a 3D 6-Node cohesive element object
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
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    08/04/15  B. Y. Chen            Original code
!

use parameter_module, only : NST => NST_COHESIVE, NDIM, DP, ZERO, HALF, ONE, &
                      & ONE_SIXTH, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE
! list of external modules used in type definition and other procedures:
! global clock module    : needed in element definition, extract and integrate
! lamina material module : needed in element definition, extract and integrate
use global_clock_module
use cohesive_material_module

implicit none
private

! list of private parameters in this module:
! NNODE         : no. of nodes in this element
! NIGPOINT      : no. of integration points in this element
! NDOF          : no. of degree-of-freedom  in this element
!
integer, parameter :: NNODE=6, NIGPOINT=3, NDOF=NDIM*NNODE


type, public :: coh3d6_element
  private
  ! list of type components:
  ! fstat         : element failure status
  ! connec        : indices of element nodes in the global node array
  ! ID_matlist    : index of element material in its material list
  ! ig_points     : integration points of this element
  ! local_clock   : locally-saved program clock
  ! traction      : traction on the interface, for output
  ! separation    : separation on the interface, for output
  ! dm            : matrix degradation factor for output
  integer  :: fstat         = 0
  integer  :: connec(NNODE) = 0
  integer  :: ID_matlist    = 0
  type(program_clock)     :: local_clock
  type(cohesive_ig_point) :: ig_points(NIGPOINT)
  real(DP) :: traction(NST)   = ZERO
  real(DP) :: separation(NST) = ZERO
  real(DP) :: dm              = ZERO
end type




interface empty
  module procedure empty_coh3d6_element
end interface

interface set
  module procedure set_coh3d6_element
end interface

interface integrate
  module procedure integrate_coh3d6_element
end interface

interface extract
  module procedure extract_coh3d6_element
end interface




public :: empty, set, integrate, extract




contains




pure subroutine empty_coh3d6_element (elem)
! Purpose:
! this subroutine is used to format the element for use
! it is used in the initialize_lib_elem procedure in the lib_elem module

  type(coh3d6_element), intent(inout) :: elem

  ! local variable, derived type var. is initialized upon declaration
  type(coh3d6_element) :: elem_lcl

  ! reset elem to the initial state
  elem = elem_lcl

end subroutine empty_coh3d6_element



pure subroutine set_coh3d6_element (elem, connec, ID_matlist, istat, emsg)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
! note that only some of the components need to be set during preproc,
! namely connec, ID_matlist

  type(coh3d6_element),   intent(inout)   :: elem
  integer,                intent(in)      :: connec(NNODE)
  integer,                intent(in)      :: ID_matlist
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  
  istat = STAT_SUCCESS
  emsg  = ''
  
  ! check validity of inputs
  if ( any(connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'connec node indices must be >=1, set, &
    &coh3d6_element_module'
    return
  end if
  
  if ( ID_matlist < 1 ) then
    istat = STAT_FAILURE
    emsg  = 'ID_matlist must be >=1, set, &
    &coh3d6_element_module'
    return
  end if
  
  elem%connec    = connec
  elem%ID_matlist = ID_matlist

end subroutine set_coh3d6_element



pure subroutine extract_coh3d6_element (elem, fstat, connec, ID_matlist, &
& local_clock, ig_points, traction, separation, dm)
! Purpose:
! to extract the components of this element
! note that the dummy args connec and ig_points are allocatable arrays
! because their sizes vary with different element types

  type(coh3d6_element),                           intent(in)  :: elem
  integer,                              optional, intent(out) :: fstat
  integer,                 allocatable, optional, intent(out) :: connec(:)
  integer,                              optional, intent(out) :: ID_matlist
  type(program_clock),                  optional, intent(out) :: local_clock
  type(cohesive_ig_point), allocatable, optional, intent(out) :: ig_points(:)
  real(DP),                             optional, intent(out) :: traction(NST)
  real(DP),                             optional, intent(out) :: separation(NST)
  real(DP),                             optional, intent(out) :: dm

  if (present(fstat))       fstat = elem%fstat

  if (present(connec)) then
    allocate(connec(NNODE))
    connec = elem%connec
  end if

  if (present(ID_matlist))  ID_matlist  = elem%ID_matlist

  if (present(local_clock)) local_clock = elem%local_clock

  if (present(ig_points)) then
    allocate(ig_points(NIGPOINT))
    ig_points = elem%ig_points
  end if

  if (present(traction))    traction    = elem%traction

  if (present(separation))  separation  = elem%separation

  if (present(dm))          dm          = elem%dm

end subroutine extract_coh3d6_element



pure subroutine integrate_coh3d6_element (elem, K_matrix, F_vector, istat, &
& emsg, nofailure, mnodes)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element
! note that K and F are allocatable dummy args as their sizes vary with
! different element types

! list of used modules:
! xnode_module                  : xnode derived type and its assoc. procedures
! global_node_list_module       : global node list
! global_material_list_module   : global material list
! global_toolkit_module         : global tools for element integration
use xnode_module
use global_node_list_module,     only : global_node_list
use global_material_list_module, only : global_cohesive_list
use global_toolkit_module,       only : cross_product3d, normalize_vect

  ! mnodes: material (interpolated) nodes
  type(coh3d6_element),     intent(inout) :: elem
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,       optional,  intent(in)    :: nofailure
  type(xnode),   optional,  intent(in)    :: mnodes(NNODE)

  ! local copies of intent(inout) dummy arg./its components:
  ! - elfstat         : elem's fstat
  ! - connec          : elem's connectivity matrix
  ! - ID_matlist      : elem's ID_matlist
  ! - local_clock     : elem's local_clock
  ! - ig_points       : elem's ig_points
  ! - eltraction      : elem's traction
  ! - elseparation    : elem's separation
  ! - eldm            : elem's dm
  integer             :: elfstat
  integer             :: connec(NNODE)
  integer             :: ID_matlist
  type(program_clock) :: local_clock
  type(cohesive_ig_point) :: ig_points(NIGPOINT)
  real(DP)            :: eltraction(NST)
  real(DP)            :: elseparation(NST)
  real(DP)            :: eldm

  ! the rest are all local variables

  !** nodal variables:
  ! - nodes           : array of element nodes
  ! - xj, uj          : nodal x and u extracted from nodes array
  ! - coords          : element nodal coordinates matrix
  ! - u               : element nodal displacemet vector
  ! - midcoords       : coordinates of the mid-plane
  type(xnode)         :: nodes(NNODE)
  real(DP), allocatable :: xj(:), uj(:)
  real(DP)            :: coords(NDIM,NNODE), u(NDOF)
  real(DP)            :: midcoords(NDIM,NNODE/2)

  !** local coords and rotational matrix
  ! - normal, tangent1/2 : normal and tangent vectors of the interface,
  !                     obtained from coords matrix
  ! - is_zero_vect    : true of any of the above vectors has zero length
  ! - Qmatrix         : rotation matrix from global to local coordinates
  !                     abtained from normal & tangent vectors
  real(DP)            :: normal(NDIM), tangent1(NDIM), tangent2(NDIM)
  logical             :: is_zero_vect
  real(DP)            :: Qmatrix(NDIM,NDIM)

  !** element determinant
  ! - det             : determinant of jacobian
  real(DP)            :: det

  !** material definition:
  ! - this_mat        : the lamina material definition of this element
  type(cohesive_material) :: this_mat

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
  connec          = 0
  ID_matlist      = 0
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
  !** element determinant
  det             = ZERO
  !** analysis logical control variables:
  last_converged  = .false.
  nofail          = .false.
  !** integration point variables:
  ig_xi           = ZERO
  ig_weight       = ZERO
  !** variables needed for stiffness matrix derivation:
  fn              = ZERO
  Nmatrix         = ZERO
  ujump           = ZERO
  delta           = ZERO
  dee             = ZERO
  tau             = ZERO
  QN              = ZERO
  QNtDQN          = ZERO
  QNtTau          = ZERO
  !** integer counter variables:
  i=0; j=0; k=0; kig=0


  ! check validity of input/imported variables
  ! here, check if global_node_list and global_lamina_list are allocated
  if (.not. allocated(global_node_list)) then
    istat = STAT_FAILURE
    emsg  = 'global_node_list not allocated, coh3d6_element_module'
  else if (.not. allocated(global_cohesive_list)) then
    istat = STAT_FAILURE
    emsg  = 'global_cohesive_list not allocated, coh3d6_element_module'
  end if
  ! if there's any error encountered above
  ! clean up and exit the program
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

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
  connec        = elem%connec       ! in
  ID_matlist    = elem%ID_matlist   ! in
  local_clock   = elem%local_clock  ! inout
  ig_points     = elem%ig_points    ! inout

  !** nodal variables:
  nodes = global_node_list(connec)
  if (present(mnodes)) then
  ! - extract nodes from passed-in node array
    nodes = mnodes
  else
  ! - extract nodes from global node list
    nodes = global_node_list(connec)
  end if
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
      emsg  = 'x not allocated for node, coh3d6_element_module'
      exit
    end if
    ! assign nodal displacement values (uj) to u vector
    if(allocated(uj)) then
      u((j-1)*NDIM+1:j*NDIM)=uj(1:NDIM)
      deallocate(uj)
    else
      istat = STAT_FAILURE
      emsg  = 'u not allocated for node, coh3d6_element_module'
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
  tangent1(1)=midcoords(1,2)-midcoords(1,1)
  tangent1(2)=midcoords(2,2)-midcoords(2,1)
  tangent1(3)=midcoords(3,2)-midcoords(3,1)
  ! compute tangent2 of the interface: node 3 coords - node 1 coords
  tangent2(1)=midcoords(1,3)-midcoords(1,1)
  tangent2(2)=midcoords(2,3)-midcoords(2,1)
  tangent2(3)=midcoords(3,3)-midcoords(3,1)
  ! compute normal vector of the interface,
  normal = cross_product3d(tangent1,tangent2)
  ! - re-evaluate tangent1 so that it is perpendicular to both
  ! tangent2 and normal
  tangent1 = cross_product3d(tangent2,normal)
  ! - normalize these vectors, exit w error if ZERO vector is encountered
  ! the magnitude of normal is twice the area of the triangle, 2 * Atri
  ! note that determinant of linear tri elem is constant and is equal to
  ! Atri / Atri_reference = Atri / 0.5 = 2 * Atri
  ! so the length of vector normal = determinant of this elem
  call normalize_vect(normal, is_zero_vect, det)
  if (is_zero_vect) then
    istat = STAT_FAILURE
    emsg  = 'element area is ZERO, coh3d6 element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  call normalize_vect(tangent1, is_zero_vect)
  if (is_zero_vect) then
    istat = STAT_FAILURE
    emsg  = 'element area is ZERO, coh3d6 element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  call normalize_vect(tangent2, is_zero_vect)
  if (is_zero_vect) then
    istat = STAT_FAILURE
    emsg  = 'element area is ZERO, coh3d6 element module'
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  ! - compute Q matrix
  do j=1,3
    Qmatrix(1,j)=normal(j)
    Qmatrix(2,j)=tangent1(j)
    Qmatrix(3,j)=tangent2(j)
  end do

  !** material definition:
  ! extract material definition from global material list
  this_mat = global_cohesive_list(ID_matlist)

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
        call ddsdde (this_mat, dee=dee, traction=tau, separation=delta)
      else
      ! failure is allowed, use ddsdde_cohesive
        call ddsdde (this_mat, dee=dee, traction=tau, sdv=ig_sdv_iter, &
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
      elfstat  = max(elfstat, ig_sdv_iter%fstat)

      ! update elem dm to be the max current ig point dm
      eldm     = max(eldm,    ig_sdv_iter%dm)

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

  end do loop_igpoint !-looped over all int points. ig=NIGPOINT

  ! check to see if the loop is exited upon error
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

  !**** END MAIN CALCULATIONS ****


  ! check to see if intent(inout) derived type dummy arg's components are
  ! unintentionally modified
  if ( any (connec /= elem%connec) ) then
    istat = STAT_FAILURE
    emsg  = 'elem%connec is unintentionally modified in coh3d6_element module'
  else if (ID_matlist /= elem%ID_matlist) then
    istat = STAT_FAILURE
    emsg  = 'elem%ID_matlist is unintentionally modified in coh3d6_element module'
  end if
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

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

end subroutine integrate_coh3d6_element






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




end module coh3d6_element_module
