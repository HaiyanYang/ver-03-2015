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

use parameter_module, only : NST => NST_COHESIVE, NDIM, DP, ZERO,            &
                      & MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,               &
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



pure subroutine set_coh3d6_element (elem, connec, ID_matlist)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
! note that only some of the components need to be set during preproc,
! namely connec, ID_matlist

  type(coh3d6_element),   intent(inout)   :: elem
  integer,                intent(in)      :: connec(NNODE)
  integer,                intent(in)      :: ID_matlist
  
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
  real(DP)                              optional, intent(out) :: traction(NST)
  real(DP)                              optional, intent(out) :: separation(NST)
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
use global_toolkit_module,       only : 

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
  real(DP)            :: normal(NDIM), tangent1(NDIM), tangent2(NDIM)
  real(DP)            :: Qmatrix(NDIM,NDIM)
  
  !** element determinant
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
  ! - ig_x, ig_u      : temporary x and u vectors for an ig point
  ! - ig_strain, ig_stress : temporary strain and stress vectors for an ig point
  real(DP)            :: ig_xi(NDIM-1, NIGPOINT), ig_weight(NIGPOINT)
  type(cohesive_sdv)  :: ig_sdv_conv, ig_sdv_iter
  real(DP)            :: ig_x(NDIM), ig_u(NDIM)
  real(DP)            :: ig_traction(NST), ig_separation(NST)

  !** variables needed for stiffness matrix derivation:
  ! - fn              : shape functions
  ! - Nmatrix         : obtained from fn to compute disp. jump vector ujump 
  !                     across the interface: {u}_jump = [N]*{u}
  ! - normal, tangent1/2 : normal and tangent vectors of the interface, 
  !                     obtained from coords matrix
  ! - det             : determinant of jacobian
  ! - Qmatrix         : rotation matrix from global to local coordinates 
  !                     abtained from normal & tangent vectors
  ! - ujump, delta    : {delta}=[Q]*{u}_jump, {delta} is the separation vector.
  ! - Dee             : material stiffness matrix [D]
  ! - Tau             : {Tau}=[D]*{delta}, traction on the interface 
  ! - QN, DQN         : [Q]*[N], [D]*[Q]*[N]
  ! - NtQt, NtQtDQN   : [N']*[Q'], [N']*[Q']*[D]*[Q]*[N] 
  ! - NtQtTau         : [N']*[Q']*{Tau}
  real(DP)            :: fn(NNODE)
  real(DP)            :: Nmatrix(NDIM,NDOF)
  real(DP)            :: ujump(NDIM), delta(NDIM)
  real(DP)            :: Dee(NDIM,NDIM)
  real(DP)            :: Tau(NDIM)
  real(DP)            :: QN(NDIM,NDOF), DQN(NDIM,NDOF)
  real(DP)            :: NtQt(NDOF,NDIM), NtQtDQN(NDOF,NDOF)
  real(DP)            :: NtQtTau(NDOF)

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
  ply_angle       = ZERO
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
  Qmatrix         = ZERO
  !** element determinant
  det             = ZERO
  !** analysis logical control variables:
  last_converged  = .false.
  nofail          = .false.
  !** integration point variables:
  ig_xi           = ZERO
  ig_weight       = ZERO
  ig_x            = ZERO
  ig_u            = ZERO
  ig_traction     = ZERO
  ig_separation   = ZERO
  !** variables needed for stiffness matrix derivation:
  fn              = ZERO
  Nmatrix         = ZERO
  ujump           = ZERO
  delta           = ZERO
  Dee             = ZERO
  Tau             = ZERO
  QN              = ZERO
  DQN             = ZERO
  NtQt            = ZERO
  NtQtDQN         = ZERO
  NtQtTau         = ZERO
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
    call clean_up()
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
    call clean_up()
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
  ! its magnitude is twice the area of the triangle, 2 * A_actual
  ! determinant of linear triangular elem is constant and is equal to
  ! A_actual / A_reference = A_actual / 0.5 = 2 * A_actual
  ! so the length of vector normal = determinant of this elem 
  normal = cross_product3d(tangent1,tangent2)
  ! - re-evaluate tangent1 so that it is perpendicular to both 
  ! tangent2 and normal
  tangent1 = cross_product3d(tangent2,normal)
  ! - normalize these vectors, exit w error if zero vector is encountered
  call normalize(normal, is_zero_vect, det)
  if (is_zero_vect) then 
    istat = STAT_FAILURE
    emsg  = 'element area is zero, coh3d6 element module'
    call clean_up()
    return
  end if
  call normalize(tangent1, is_zero_vect)
  if (is_zero_vect) then 
    istat = STAT_FAILURE
    emsg  = 'element area is zero, coh3d6 element module'
    call clean_up()
    return
  end if
  call normalize(tangent2, is_zero_vect)
  if (is_zero_vect) then 
    istat = STAT_FAILURE
    emsg  = 'element area is zero, coh3d6 element module'
    call clean_up()
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
  ! calculate ig point xi and weight; other ig point variables are updated
  ! later in the loop
  call init_ig_point (ig_xi, ig_weight)
  
  !** other local variables are assigned in the following loop of all ig points



  !**** MAIN CALCULATIONS ****


  ! loop over all ig points,
  ! calculate strain, stress, stiffness matrix, sdv etc. at each ig point

  loop_igpoint: do kig=1, NIGPOINT 
      
      !- get shape function matrix
      call init_shape(fn, igxi(:,kig)) 
      
      ! Nmatrix: ujump (at each int pnt) = Nmatrix*u
      do i = 1, NDIM
        do j = 1, NNODE/2
          Nmatrix( i, i + (j-1)         * NDIM ) = - fn(j)
          Nmatrix( i, i + (j-1+NNODE/2) * NDIM ) =   fn(j)
        end do
      end do
      
      ! calculate ujump: disp. jump of the two crack surface, in global coords
      ujump=matmul(Nmatrix,u)
      
      ! calculate separation delta in local coords: delta=Qmatrix*ujump
      delta=matmul(Qmatrix,ujump)
      
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
      
      ! get D matrix dee accord. to material properties, and update intg point variables
      if(nofail) then
      ! do not pass in sdv, then no cohesive law can be done, only linear elasticity stiffness and stress calculated
          call ddsdde(this_mat, Dee, jump=delta, stress=Tau)
      else
      ! calculate D matrix, update stress, and iterating sdv
          call ddsdde(this_mat, Dee, jump=delta, stress=Tau, sdv=ig_sdv(2), dfail=dfail) 
      end if
      
      !------------------------------------------------!
      !  add this ig point contributions to K and F
      !------------------------------------------------!
      QN      =   matmul(Qmatrix, Nmatrix)
      NtQt    =   transpose(QN)
      DQN     =   matmul(Dee, QN)
      NtQtDQN =   matmul(NtQt, DQN)
      NtQtTau =   matmul(NtQt, Tau)
      do j=1, NDOF
        do k=1, NDOF
          K_matrix(j,k) = K_matrix(j,k) + NtQtDQN(j,k) * det * ig_weight(kig)
        end do
        F_vector(j) = F_vector(j) + NtQtTau(j) * ig_weight(kig) * det
      end do
         
      !------------------------------------------------!
      !   update ig point kig
      !------------------------------------------------!
      
      !- calculate integration point physical coordinates (initial)
      ig_x = matmul(coords,fn)
      !- calculate integration point displacement
      do j=1, NDIM
        do i=1, NNODE
          ig_u(j) = ig_u(j) + fn(i) * u((i-1)*NDIM+j)
        end do
      end do
      
      ! update element ig point arrays
      !call update(elem%ig_points(kig),x=tmpx,u=tmpu,strain=delta,stress=Tau,sdv=ig_sdv)
      call update(elem%ig_points(kig),sdv=ig_sdv)
      
      ! update elem curr status variable
      igstat=0
      if(allocated(ig_sdv(2)%i)) igstat=ig_sdv(2)%i(1)
      elem%fstat=max(elem%fstat,igstat)
      
      !~! update elem sdv arrays
      !~if(.not.allocated(elem%sdv)) then
      !~    allocate(elem%sdv(1))
      !~    elem%sdv(1)=ig_sdv(2)
      !~else
      !~    if(ig_sdv(2)%i(1) > elem%sdv(1)%i(1))
      !~    ! if fstat value is higher, update elem sdv to this ig point sdv values
      !~        elem%sdv(1)=ig_sdv(2)
      !~    else if(ig_sdv(2)%i(1)==elem%sdv(1)%i(1) .and. ig_sdv(2)%r(1) > elem%sdv(1)%r(1))
      !~    ! if fstat value is the same but dm value is higher, also update to this ig point sdv
      !~        elem%sdv(1)=ig_sdv(2)
      !~    else
      !~    ! no update
      !~        continue
      !~    end if
      !~end if
      
      ! - empty relevant arrays for reuse
      fn=zero; Nmatrix=zero
      ujump=zero; delta=zero; dee=zero; Tau=zero
      QN=zero; NtQt=zero; DQN=zero; NtQtDQN=zero; NtQtTau=zero
      tmpx=zero; tmpu=zero
      
      
  end do loop_igpoint !-looped over all int points. ig=NIGPOINT

  if(allocated(xj)) deallocate(xj) 
  if(allocated(uj)) deallocate(uj)
  
  
  
  contains
  ! internal procedures
  
    subroutine clean_up()
      ! zero intent(out) variables (do not deallocate)
      K_matrix = ZERO
      F_vector = ZERO
      ! deallocate local alloc. variables
      if (allocated(uj)) deallocate(uj)
      if (allocated(xj)) deallocate(xj)
    end subroutine clean_up

end subroutine integrate_coh3d6_element












! the rest are private subroutines









pure subroutine init_ig (xi, wt, gauss)

  real(DP), intent(inout)  :: xi(NDIM-1,NIGPOINT), wt(NIGPOINT)
  logical, optional, intent(in) :: gauss
  

    

    if (NIGPOINT .eq. 3) then   

        if(present(gauss).and.gauss) then
            xi(1,1)=half
            xi(2,1)=half
            
            xi(1,2)=half
            xi(2,2)=zero
            
            xi(1,3)=zero
            xi(2,3)=half
        else
            xi(1,1)=zero
            xi(2,1)=zero
            
            xi(1,2)=one
            xi(2,2)=zero
            
            xi(1,3)=zero
            xi(2,3)=one
        end if
        wt=one_sixth
        
    else
        write(msg_file,*) 'no. of integration points incorrect for coh3d6_ig!'
        call exit_function
    end if

end subroutine init_ig



pure subroutine init_shape (igxi, f)
  
    real(DP),intent(inout) :: f(NNODE)
    real(DP),intent(in) :: igxi(NDIM-1)
    
    real(DP) :: xi, eta ! local variables
    xi=zero
    eta=zero

    xi=igxi(1)
    eta=igxi(2)
    f(1)=one-xi-eta
    f(2)=xi
    f(3)=eta
    f(4)=f(1)
    f(5)=f(2)
    f(6)=f(3)

end subroutine init_shape




end module coh3d6_element_module
