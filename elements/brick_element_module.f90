module brick_element_module
!
!  Purpose:
!    define a brick element object with lamina material definition
!    with the associated procedures to empty, set, integrate and extract
!    its components
!
!  topological definition of this brick element, local nodal indices:
!
!  8___________________7
!  |\                  |\
!  | \                 | \
!  |  \                |  \
!  |   \               |   \
!  |    \              |    \
!  |     \5____________|_____\6
!  |______|____________|      |
! 4\      |   E3      3\      |
!   \     |             \     |
!    \    |              \    |
!   E4\   |             E2\   | 
!      \  |                \  |
!       \ |                 \ |
!        \|__________________\|
!         1         E1        2
!
!  bottom surface nodes (counter-clock vise): 1, 2, 3, 4
!  top    surface nodes (counter-clock vise): 5, 6, 7, 8
!  bottom surface edges (counter-clock vise): E1, E2, E3, E4
!  end nodes of E1 : 1, 2
!  end nodes of E2 : 2, 3
!  end nodes of E3 : 3, 4
!  end nodes of E4 : 4, 1
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    08/04/15  B. Y. Chen            Original code
!
use parameter_module, only : NDIM, NST => NST_STANDARD, DP, SDV_ARRAY, ZERO,
use global_clock_module
use lamina_material_module

implicit none
private

! list of private parameters in this module:
! NNODE         : no. of nodes in this element
! NEDGE_BOTTOM  : no. of edges in this element on the bottom surface
! NIGPOINT      : no. of integration points in this element
! NDOF          : no. of degree-of-freedom  in this element
!
! NODE_ON_BOTTOM_EDGE(i,j) : this is a matrix to describe the nodal connectivity
!                            of the bottom edges in this element. Each edge
!                            topologically composed of its two end nodes.
!                            This term stores the local nodal index of  
!                            the i_th node on edge j
!
! ** Note: only the bottom edge nodal connec is actually needed; top edge nodal
!          connec can be derived as:
!          NODE_ON_TOP_EDGE(i,j) = NODE_ON_BOTTOM_EDGE(i,j) + NNODE/2
!
integer, parameter :: NIGPOINT=8, NNODE=8, NEDGE_BOTTOM=4, NDOF=NDIM*NNODE
integer, parameter :: NODE_ON_BOTTOM_EDGE(2, NEDGE_BOTTOM) = &
                    & reshape([ 1,2, 2,3, 3,4, 4,1 ], [2, NEDGE_BOTTOM])


type, public :: brick_element 
  private
  ! list of type components:
  ! curr_status   : element current status
  ! connec        : indices of element nodes in the global node array
  ! ID_matlist    : index of element material in its material array
  ! ply_angle     : ply angle for composite lamina (rotation around z axis)
  ! stress        : stresses for output
  ! strain        : strains for output
  ! ig_point      : integration point array of this element
  ! local_clock   : locally-saved program clock
  ! converged_sdv : element lamina_sdv of last converged increment
  ! iterating_sdv : element lamina_sdv of current iteration
  integer  :: curr_status   = 0
  integer  :: connec(NNODE) = 0
  integer  :: ID_matlist    = 0
  real(DP) :: ply_angle     = ZERO
  real(DP) :: stress(NST)   = ZERO
  real(DP) :: strain(NST)   = ZERO
  type(program_clock)   :: local_clock
  type(lamina_ig_point) :: ig_points(NIGPOINT)
  type(lamina_sdv)      :: converged_sdv
  type(lamina_sdv)      :: iterating_sdv  
end type




interface empty
  module procedure empty_brick_element
end interface

interface set
  module procedure set_brick_element
end interface

interface integrate
  module procedure integrate_brick_element
end interface

interface extract
  module procedure extract_brick_element
end interface




public :: empty, set, integrate, extract




contains




pure subroutine empty_brick_element (elem)
! Purpose:
! this subroutine is used to format the element for use
! it is used in the initialize_lib_elem procedure in the lib_elem module

  type(brick_element), intent(inout) :: elem
  
  ! local variable, derived type var. is initialized upon declaration
  type(brick_element) :: elem_lcl
  
  ! reset elem to the initial state
  elem = elem_lcl

end subroutine empty_brick_element



pure subroutine set_brick_element (elem, connec, ID_matlist, ply_angle)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
! note that only some of the components need to be set during preproc,
! namely connec, ID_matlist, ply_angle and local_clock

  type(brick_element),    intent(inout)   :: elem
  integer,                intent(in)      :: connec(NNODE)
  integer,                intent(in)      :: ID_matlist
  real(DP),               intent(in)      :: ply_angle
  
  elem%connec    = connec
  elem%ID_matlist = ID_matlist
  elem%ply_angle = ply_angle
  elem%local_clock = GLOBAL_CLOCK

end subroutine set_brick_element



pure subroutine extract_brick_element (elem, curr_status, connec, &
& ID_matlist, ply_angle, stress, strain, local_clock, ig_points,    &
& converged_sdv, iterating_sdv)
! Purpose:
! to extract the components of this element
! note that the dummy args connec and ig_point are allocatable arrays
! because their sizes vary with different element types

  type(brick_element),                          intent(in)  :: elem  
  integer,                            optional, intent(out) :: curr_status
  integer,               allocatable, optional, intent(out) :: connec(:)
  integer,                            optional, intent(out) :: ID_matlist
  real(DP),                           optional, intent(out) :: ply_angle
  real(DP),                           optional, intent(out) :: stress(NST)
  real(DP),                           optional, intent(out) :: strain(NST)
  type(program_clock),                optional, intent(out) :: local_clock
  type(lamina_ig_point), allocatable, optional, intent(out) :: ig_points(:)
  type(lamina_sdv),                   optional, intent(out) :: converged_sdv
  type(lamina_sdv),                   optional, intent(out) :: iterating_sdv
  
  if (present(curr_status)) curr_status = elem%curr_status
  
  if (present(connec)) then
    allocate(connec(NNODE))
    connec = elem%connec
  end if
  
  if (present(ID_matlist))   ID_matlist   = elem%ID_matlist
  
  if (present(ply_angle))   ply_angle   = elem%ply_angle
  
  if (present(stress)) stress = elem%stress

  if (present(strain)) strain = elem%strain

  if (present(local_clock)) local_clock = elem%local_clock
  
  if (present(ig_points)) then
    allocate(ig_points(NIGPOINT))
    ig_points = elem%ig_points
  end if
  
  if (present(converged_sdv))  converged_sdv = elem%converged_sdv
  
  if (present(iterating_sdv))  iterating_sdv = elem%iterating_sdv

end subroutine extract_brick_element



pure subroutine integrate_brick_element (elem, K_matrix, F_vector, istat, emsg,&
& nofailure)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element
! note that K and F are allocatable dummy args as their sizes vary with
! different element types

! list of used modules:
! global tools for element integration
! global material list
! global node list
! global analysis progress (curr. step, inc, time, dtime)
use toolkit_module,      only :     
use list_mat_module,     only : list_lamina_mat  
use list_node_module,    only : list_node   
use global_clock_module, only : program_clock, GLOBAL_CLOCK, clock_in_sync   

  type(brick_element),      intent(inout) :: elem 
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure
  
  ! local copies of intent(inout) dummy arg./its components:
  ! - elstatus        : elem's curr_status
  ! - connec          : elem's connectivity matrix
  ! - ID_matlist      : elem's ID_matlist
  ! - ply_angle       : elem's ply_angle
  ! - elstress        : elem's stress
  ! - elstrain        : elem's strain
  ! - local_clock     : elem's local_clock
  ! - ig_points       : elem's ig_points
  
  ! the rest are all local variables
  
  !** nodal variables:
  ! - nodes           : array of element nodes  
  ! - xj, uj          : nodal x and u extracted from nodes array
  ! - coords          : element nodal coordinates matrix
  ! - u               : element nodal displacemet vector
  
  !** material definition:
  ! - this_mat        : the lamina material definition of this element
  
  !** element characteristic length variables:
  ! - clength, ctip   : variables used for calculating clength
  
  !** analysis logical control variables:
  ! - last_converged  : true if last iteration has converged
  ! - nofail          : true if no failure is allowed
  
  !** integration point variables:
  ! - ig_xi           : integration point natural coordinates
  ! - ig_weight       : integration point weights
  ! - ig_sdv_conv/iter: integration point converged and iterating sdvs
  ! - ig_fstat        : failure status of integration point (from its sdv)
  ! - ig_x, ig_u      : temporary x and u vectors for an ig point
  ! - ig_strain, ig_stress : temporary strain and stress vectors for an ig point
  
  !** variables needed for stiffness matrix derivation:
  ! - fn, dn          : shape functions & their derivatives in natural space
  ! - gn              : gradients of shape functions in physical space
  ! - jac, detj       : jacobian matrix and its determinant
  ! - bee             : B matrix
  ! - dee             : D matrix
  ! - btdb            : B' * D * B
  
  !** integer counter variables
  ! - i, j, kig

  ! local copies of intent(inout) dummy arg./its components:
  integer  :: elstatus
  integer  :: connec(NNODE)
  integer  :: ID_matlist
  real(DP) :: ply_angle
  real(DP) :: elstress(NST)
  real(DP) :: elstrain(NST)
  type(program_clock)   :: local_clock
  type(lamina_ig_point) :: ig_points(NIGPOINT)
  !** nodal variables:
  type(xnode)           :: nodes(NNODE)
  real(DP), allocatable :: xj(:), uj(:)
  real(DP)              :: coords(NDIM,NNODE), u(NDOF)
  !** material definition:
  type(lamina_material) :: this_mat
  !** element characteristic length variables: 
  real(DP) :: clength, ctip(2,2)
  !** analysis logical control variables:
  logical  :: last_converged, nofail
  !** integration point variables:
  real(DP)         :: ig_xi(NDIM, NIGPOINT), ig_weight(NIGPOINT)
  type(lamina_sdv) :: ig_sdv_conv, ig_sdv_iter
  integer          :: ig_fstat
  real(DP)         :: ig_x(NDIM), ig_u(NDIM), ig_strain(NST), ig_stress(NST)
  !** variables needed for stiffness matrix derivation:
  real(DP) :: fn(NNODE), dn(NNODE,NDIM), gn(NNODE,NDIM)
  real(DP) :: jac(NDIM,NDIM), detj
  real(DP) :: bee(NST,NDOF), dee(NST,NST), btdb(NDOF,NDOF)
  !** integer counter variables:
  integer  :: i, j, kig
 
  
  ! initialize intent(out) and local variables
  ! except derived types (auto initialized upon declaration) and 
  ! allocatable local arrays
  
  !** intent(out) variables:
  allocate(K_matrix(NDOF,NDOF), F_vector(NDOF))
  K_matrix = ZERO
  F_vector = ZERO
  istat = 0
  emsg  = ''
  !** local copies of intent(inout) dummy args./its components:
  elstatus   = 0
  connec     = 0
  ID_matlist = 0
  ply_angle  = ZERO
  elstress   = ZERO
  elstrain   = ZERO
  !** nodal variables:
  call empty(nodes)
  coords    = ZERO
  u         = ZERO
  !** element characteristic length variables: 
  clength = ZERO
  ctip    = ZERO
  !** analysis logical control variables:
  last_converged = .false.
  nofail         = .false.
  !** integration point variables:
  ig_xi     = ZERO
  ig_weight = ZERO
  ig_fstat  = 0
  ig_x      = ZERO
  ig_u      = ZERO
  ig_strain = ZERO
  ig_stress = ZERO
  !** variables needed for stiffness matrix derivation:
  fn        = ZERO
  dn        = ZERO
  gn        = ZERO
  jac       = ZERO
  detj      = ZERO
  bee       = ZERO
  dee       = ZERO
  btdb      = ZERO
  !** integer counter variables:
  i=0; j=0; kig=0
  
  
  ! assign values to local variables
  
  !** local copies of intent(inout) dummy arg.:
  ! elem's curr_status, stress and strain components are NOT needed 
  ! for calculations in this subroutine. they are only used for output, 
  ! and they need to be assgined new values during each iteration.
  ! they are therefore NOT passed to their local copies elstatus, elstress 
  ! and elstrain.
  ! elstatus  : no extraction     ! out
  ! elstress  : no extraction     ! out
  ! elstrain  : no extraction     ! out
  connec      = elem%connec       ! in
  ID_matlist  = elem%ID_matlist   ! in
  ply_angle   = elem%ply_angle    ! in
  local_clock = elem%local_clock  ! inout
  ig_points   = elem%ig_points    ! inout
  
  !** nodal variables:
  ! copy nodes from global nodes array to local nodes array 
  ! using element node indices stored in connectivity array
  nodes = list_node(connec)
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
      emsg  = 'x not allocated for node, brick_element_module'
      exit
    end if
    ! assign nodal displacement values (uj) to u vector
    if(allocated(uj)) then 
      u((j-1)*NDIM+1:j*NDIM)=uj(1:NDIM)
      deallocate(uj)
    else
      istat = STAT_FAILURE
      emsg  = 'u not allocated for node, brick_element_module'
      exit
    end if 
  end do
  if (istat == STAT_FAILURE) then
    ! zero intent(out) variables (do not deallocate)
    K_matrix = ZERO
    F_vector = ZERO
    ! deallocate local alloc. variables
    if (allocated(uj)) deallocate(uj)
    if (allocated(xj)) deallocate(xj)
    ! exit program
    return
  end if
  
  !** material variables:
  ! extract material definition from global material list
  this_mat = list_lamina_mat(ID_matlist)
  
  !** element characteristic length variables: 
  call elem_ctips_origin('brick', theta, coords, &
  & NODE_ON_BOTTOM_EDGE, NEDGE_BOTTOM, ctip)
  clength=sqrt((ctip(1,2)-ctip(1,1))**2+(ctip(2,2)-ctip(2,1))**2)
  
  !** analysis logical control variables:
  ! check if last iteration has converged by checking if the global clock has 
  ! advanced; if so, last iteration is converged and sync the local clock
  if (.not. clock_in_sync(GLOBAL_CLOCK, local_clock)) then
    last_converged = .true.
    local_clock    = GLOBAL_CLOCK
  end if
  ! nofail:
  if(present(nofailure)) nofail = nofailure

  !** integration point variables:
  ! update ig point xi and weight
  call init_ig(ig_xi,ig_weight)
  
  !** other local variables are assigned in the following loop of all ig points
  
  
  ! loop over all ig points, 
  ! calculate strain, stress, stiffness matrix, sdv etc. at each ig point
  loop_igpoint: do kig = 1, NIGPOINT 
    
    ! get shape matrix and derivatives
    call init_shape(ig_xi(:,kig), fn, dn) 
    
    ! calculate integration point physical coordinates (initial)
    ig_x = matmul(coords,fn)
    
    ! calculate integration point displacement
    do j=1, NDIM
      do i=1, NNODE
        ig_u(j) = ig_u(j) + fn(i) * u((i-1)*NDIM+j)
      end do
    end do
    
    ! get jacobian
    jac = matmul(coords,dn)
    
    ! get determinant of jacobian
    detj = determinant(jac)
    
    ! invert jac onto itself
    jac = inverse(jac,detj)
    
    ! calculate gradient of shape function matrix
    gn = matmul(dn,jac)
    
    ! obtain b matrix (NST*NDOF) from rearranging terms of gn
    bee = beemat(gn)
    
    ! calculate global strains at this ig point
    ig_strain = matmul(bee,u)
    
    ! transfer strain to local coordinates when ply_angle is non-zero
    if(abs(ply_angle) > SMALLNUM) ig_strain = lcl_strain(ig_strain,ply_angle)
    
    ! extract sdvs from integration points
    call extract(ig_point(kig), converged_sdv=ig_sdv_conv, &
    & iterating_sdv=ig_sdv_iter)
    
    ! update converged sdv with iterating sdv when last iteration is converged
    ! and update iterating sdv back to the last converged sdv if otherwise
    if(last_converged) then
      ig_sdv_conv = ig_sdv_iter
    else               
      ig_sdv_iter = ig_sdv_conv
    end if
    
    ! get D matrix dee accord. to material properties, sdv, strain and clength
    ! and update ig point stresses, sdv
    if(nofail) then
      ! ddsdde_lamina_intact
      call ddsdde (this_mat, dee=dee, stress=ig_stress, strain=ig_strain)
    else
      ! ddsdde_lamina
      call ddsdde (this_mat, dee=dee, stress=ig_stress, sdv=ig_sdv_iter, &
      & strain=ig_strain, clength=clength, istat=istat, emsg=emsg)
    end if

    ! get D matrix in global coordinates
    if(abs(ply_angle) > SMALLNUM) dee = glb_dee(dee,ply_angle)
    
    ! calculate B' D B
    btdb = matmul( matmul(transpose(bee),dee), bee )

    ! integrate and update K matrix
    do i=1, NDOF
      do j=1, NDOF
        K_matrix(i,j) = K_matrix(i,j) + btdb(i,j) * detj * ig_weight(kig)
      end do
    end do	
    
    ! update ig point arrays
    call update(ig_point(kig), x=ig_x, u=ig_u,                        &
    &           stress=ig_stress, strain=ig_strain,                   &
    &           converged_sdv=ig_sdv_conv, iterating_sdv=ig_sdv_iter)
    
    ! update elem curr status variable
    ig_fstat = ig_sdv_iter%fstat
    elstatus = max(elstatus,ig_fstat)

    ! update elem stress & strain (avg of ig point stress & strains)
    elstress = elstress + ig_stress/NIGPOINT
    elstrain = elstrain + ig_strain/NIGPOINT
    
    ! empty relevant arrays for reuse
    fn  = ZERO; dn   = ZERO; gn   = ZERO
    jac = ZERO; detj = ZERO
    bee = ZERO; dee  = ZERO; btdb = ZERO
    ig_x = ZERO; ig_u = ZERO
    ig_strain=ZERO; ig_stress=ZERO
      
  end do loop_igpoint !-looped over all int points. ig=NIGPOINT

  ! check to see if intent(inout) derived type dummy arg's components are 
  ! unintentionally modified
  if (any(connec /= elem%connec)) then
    istat = STAT_FAILURE
    emsg  = 'elem%connec is unintentionally modified in brick_element module'
  end if
  if (ID_matlist /= elem%ID_matlist) then
    istat = STAT_FAILURE
    emsg  = 'elem%ID_matlist is unintentionally modified in brick_element module'
  end if
  if (abs(ply_angle - elem%ply_angle) > SMALLNUM) then
    istat = STAT_FAILURE
    emsg  = 'elem%ply_angle is unintentionally modified in brick_element module'
  end if
  if (istat == STAT_FAILURE) then
    ! zero intent(out) dummy args
    K_matrix = ZERO
    F_vector = ZERO
    ! deallocate local variables
    if(allocated(xj)) deallocate(xj) 
    if(allocated(uj)) deallocate(uj) 
    ! exit program
    return
  end if
  
  ! if no unintentinal modification, proceed with final calculations
  ! and updates and exit the program
  
  ! calculate F_vector
  F_vector=matmul(K_matrix,u) 
  
  ! update intent(inout) dummy arg./its components
  elem%curr_status = elstatus
  elem%stress      = elstress
  elem%strain      = elstrain
  elem%local_clock = local_clock
  elem%ig_points   = ig_points
  
  ! deallocate local dynamic arrays
  if(allocated(xj)) deallocate(xj) 
  if(allocated(uj)) deallocate(uj) 
    
    
end subroutine integrate_brick_element












! the rest are private subroutines









pure subroutine init_ig (xi, wt)

  real(DP),intent(inout) :: xi(NDIM,NIGPOINT),wt(NIGPOINT)

    if (NIGPOINT .eq. 1) then
      xi(1,1)= ZERO
      xi(2,1)= ZERO
      xi(3,1)= ZERO
      wt = eight
    else if (NIGPOINT .eq. 8) then
      xi(1,1)= -root3
      xi(2,1)= -root3
      xi(3,1)= -root3
      xi(1,2)=  root3
      xi(2,2)= -root3
      xi(3,2)= -root3
      xi(1,3)= -root3
      xi(2,3)=  root3
      xi(3,3)= -root3
      xi(1,4)=  root3
      xi(2,4)=  root3
      xi(3,4)= -root3
      xi(1,5)= -root3
      xi(2,5)= -root3
      xi(3,5)=  root3
      xi(1,6)=  root3
      xi(2,6)= -root3
      xi(3,6)=  root3
      xi(1,7)= -root3
      xi(2,7)=  root3
      xi(3,7)=  root3
      xi(1,8)=  root3
      xi(2,8)=  root3
      xi(3,8)=  root3          
      wt = one         
    else
        write(msg_file,*) 'no. of integration points incorrect for brick_ig!'
        call exit_function
    end if

end subroutine init_ig



pure subroutine init_shape (ig_xi, f, df)
  
    real(DP),intent(inout) :: f(NNODE),df(NNODE,NDIM)
    real(DP),intent(in) :: ig_xi(NDIM)
    
    real(DP) :: xi,eta,zeta ! local variables
    xi=ZERO
    eta=ZERO
    zeta=ZERO
    

    xi=ig_xi(1)
    eta=ig_xi(2)
    zeta=ig_xi(3)
    
    f(1)=one_eighth*(one-xi)*(one-eta)*(one-zeta)
    f(2)=one_eighth*(one+xi)*(one-eta)*(one-zeta)
    f(3)=one_eighth*(one+xi)*(one+eta)*(one-zeta)
    f(4)=one_eighth*(one-xi)*(one+eta)*(one-zeta)
    f(5)=one_eighth*(one-xi)*(one-eta)*(one+zeta)
    f(6)=one_eighth*(one+xi)*(one-eta)*(one+zeta)
    f(7)=one_eighth*(one+xi)*(one+eta)*(one+zeta)
    f(8)=one_eighth*(one-xi)*(one+eta)*(one+zeta)
    
    df(1,1) = -one_eighth*(one-eta)*(one-zeta)
    df(2,1) =  one_eighth*(one-eta)*(one-zeta)
    df(3,1) =  one_eighth*(one+eta)*(one-zeta)
    df(4,1) = -one_eighth*(one+eta)*(one-zeta)
    df(5,1) = -one_eighth*(one-eta)*(one+zeta)
    df(6,1) =  one_eighth*(one-eta)*(one+zeta)
    df(7,1) =  one_eighth*(one+eta)*(one+zeta)
    df(8,1) = -one_eighth*(one+eta)*(one+zeta)
    
    df(1,2) = -one_eighth*(one-xi)*(one-zeta)
    df(2,2) = -one_eighth*(one+xi)*(one-zeta)
    df(3,2) =  one_eighth*(one+xi)*(one-zeta)
    df(4,2) =  one_eighth*(one-xi)*(one-zeta)
    df(5,2) = -one_eighth*(one-xi)*(one+zeta)
    df(6,2) = -one_eighth*(one+xi)*(one+zeta)
    df(7,2) =  one_eighth*(one+xi)*(one+zeta)
    df(8,2) =  one_eighth*(one-xi)*(one+zeta)
    
    df(1,3) = -one_eighth*(one-xi)*(one-eta)
    df(2,3) = -one_eighth*(one+xi)*(one-eta)
    df(3,3) = -one_eighth*(one+xi)*(one+eta)
    df(4,3) = -one_eighth*(one-xi)*(one+eta)
    df(5,3) =  one_eighth*(one-xi)*(one-eta)
    df(6,3) =  one_eighth*(one+xi)*(one-eta)
    df(7,3) =  one_eighth*(one+xi)*(one+eta)
    df(8,3) =  one_eighth*(one-xi)*(one+eta)
    

end subroutine init_shape




end module brick_element_module
