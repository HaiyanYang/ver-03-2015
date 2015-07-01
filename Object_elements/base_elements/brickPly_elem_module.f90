module brickPly_elem_module
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
!    21/04/15  B. Y. Chen            removed access of glb node/mat lists, 
!                                    remove ID_matlist, pass mat and nodes arg.
!                                    to integrate procedure
!

use parameter_module, only : NST => NST_STANDARD, NDIM, DP,                  &
                      & MSGLENGTH, STAT_SUCCESS, STAT_FAILURE,               &
                      & ZERO, ONE, EIGHT, ONE_ROOT3, ONE_EIGHTH, SMALLNUM
! list of external modules used in type definition and other procedures:
! global clock module    : needed in element definition, extract and integrate
! lamina material module : needed in element definition, extract and integrate
use global_clock_module,    only : program_clock, GLOBAL_CLOCK, clock_in_sync
use lamina_material_module, only : lamina_ig_point, lamina_material, lamina_sdv,&
                            & ddsdde, update, extract

implicit none
private

! list of private parameters in this module:
! NNODE         : no. of nodes in this element
! NEDGE_BOTTOM  : no. of edges in this element on the bottom surface
! NIGPOINT      : no. of integration points in this element
! NDOF          : no. of degree-of-freedom  in this element
!
! NODES_ON_BOTTOM_EDGES(i,j) : this is a matrix to describe the nodal connectivity
!                            of the bottom edges in this element. Each edge
!                            topologically composed of its two end nodes.
!                            This term stores the local nodal index of
!                            the i_th node on edge j
!
! ** Note: only the bottom edge nodal connec is actually needed; top edge nodal
!          connec can be derived as:
!          NODE_ON_TOP_EDGE(i,j) = NODES_ON_BOTTOM_EDGES(i,j) + NNODE/2
!
integer, parameter :: NIGPOINT=8, NNODE=8, NEDGE_BOTTOM=4, NDOF=NDIM*NNODE
integer, parameter :: NODES_ON_BOTTOM_EDGES(2, NEDGE_BOTTOM) = &
                    & reshape([ 1,2, 2,3, 3,4, 4,1 ], [2, NEDGE_BOTTOM])


type, public :: brickPly_elem
  private
  ! list of type components:
  ! fstat         : element failure status
  ! connec        : indices of elem nodes in the glb node list, need for output
  ! ply_angle     : ply angle for composite lamina (rotation around z axis)
  ! ig_points     : integration points of this element
  ! local_clock   : locally-saved program clock
  ! stress        : stresses for output
  ! strain        : strains for output
  ! df            : fibre degradation factor for output
  integer  :: fstat         = 0
  integer  :: connec(NNODE) = 0
  real(DP) :: ply_angle     = ZERO
  type(program_clock)   :: local_clock
  type(lamina_ig_point) :: ig_points(NIGPOINT)
  real(DP) :: stress(NST)   = ZERO
  real(DP) :: strain(NST)   = ZERO
  real(DP) :: df            = ZERO
end type




interface empty
  module procedure empty_brickPly_elem
end interface

interface set
  module procedure set_brickPly_elem
end interface

interface integrate
  module procedure integrate_brickPly_elem
end interface

interface extract
  module procedure extract_brickPly_elem
end interface




public :: empty, set, integrate, extract




contains




pure subroutine empty_brickPly_elem (elem)
! Purpose:
! this subroutine is used to format the element for use
! it is used in the initialize_lib_elem procedure in the lib_elem module

  type(brickPly_elem), intent(inout) :: elem

  ! local variable, derived type var. is initialized upon declaration
  type(brickPly_elem) :: elem_lcl

  ! reset elem to the initial state
  elem = elem_lcl

end subroutine empty_brickPly_elem



pure subroutine set_brickPly_elem (elem, connec, ply_angle, istat, emsg)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
! note that only some of the components need to be set during preproc,
! namely connec, ID_matlist, ply_angle

  type(brickPly_elem),      intent(inout) :: elem
  integer,                  intent(in)    :: connec(NNODE)
  real(DP),                 intent(in)    :: ply_angle
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  
  istat = STAT_SUCCESS
  emsg  = ''

  ! check validity of inputs
  if ( any(connec < 1) ) then
    istat = STAT_FAILURE
    emsg  = 'connec node indices must be >=1, set, &
    &brickPly_elem_module'
    return
  end if

  elem%connec     = connec
  elem%ply_angle  = ply_angle

end subroutine set_brickPly_elem



pure subroutine extract_brickPly_elem (elem, fstat, connec, ply_angle, &
& ig_points, stress, strain, df)
! Purpose:
! to extract the components of this element
! note that the dummy args connec and ig_points are allocatable arrays
! because their sizes vary with different element types

  type(brickPly_elem),                          intent(in)  :: elem
  integer,                            optional, intent(out) :: fstat
  integer,               allocatable, optional, intent(out) :: connec(:)
  real(DP),                           optional, intent(out) :: ply_angle
  type(lamina_ig_point), allocatable, optional, intent(out) :: ig_points(:)
  real(DP),                           optional, intent(out) :: stress(NST)
  real(DP),                           optional, intent(out) :: strain(NST)
  real(DP),                           optional, intent(out) :: df

  if (present(fstat))       fstat = elem%fstat

  if (present(connec)) then
    allocate(connec(NNODE))
    connec = elem%connec
  end if

  if (present(ply_angle))   ply_angle   = elem%ply_angle

  if (present(ig_points)) then
    allocate(ig_points(NIGPOINT))
    ig_points = elem%ig_points
  end if

  if (present(stress))      stress = elem%stress

  if (present(strain))      strain = elem%strain

  if (present(df))          df     = elem%df

end subroutine extract_brickPly_elem



pure subroutine integrate_brickPly_elem (elem, nodes, material, K_matrix, &
& F_vector, istat, emsg, nofailure)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element
! note that K and F are allocatable dummy args as their sizes vary with
! different element types

! used parameters:
! ZERO, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, SMALLNUM

! list of additionally used modules:
! fnode_module                  : fnode derived type and its assoc. procedures
! global_toolkit_module         : global tools for element integration
use fnode_module,                only : fnode, extract
use global_toolkit_module,       only : crack_elem_centroid2d, determinant3d, &
                                      & invert_self3d, beemat3d, lcl_strain3d,&
                                      & glb_dee3d, distance

  type(brickPly_elem),      intent(inout) :: elem
  type(fnode),              intent(in)    :: nodes(NNODE)
  type(lamina_material),    intent(in)    :: material
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure

  ! local copies of intent(inout) dummy arg./its components:
  ! - elfstat         : elem's fstat
  ! - ply_angle       : elem's ply_angle
  ! - local_clock     : elem's local_clock
  ! - ig_points       : elem's ig_points
  ! - elstress        : elem's stress
  ! - elstrain        : elem's strain
  ! - eldf            : elem's df
  integer             :: elfstat
  real(DP)            :: ply_angle
  type(program_clock) :: local_clock
  type(lamina_ig_point) :: ig_points(NIGPOINT)
  real(DP)            :: elstress(NST)
  real(DP)            :: elstrain(NST)
  real(DP)            :: eldf

  ! the rest are all local variables

  !** nodal variables:
  ! - xj, uj          : nodal x and u extracted from nodes array
  ! - coords          : element nodal coordinates matrix
  ! - u               : element nodal displacemet vector
  real(DP), allocatable :: xj(:), uj(:)
  real(DP)            :: coords(NDIM,NNODE), u(NDOF)

  !** element characteristic length variables:
  ! - clength         : characteristic length of this element
  ! - crosspoints     : coordinates of cross points btw characteristic line and
  !                     two element edges
  real(DP)            :: clength, crosspoints(2,2)

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
  ! - ig_fstat        : extracted fstat from ig point
  ! - ig_df           : extracted fstat from ig point
  real(DP)            :: ig_xi(NDIM, NIGPOINT), ig_weight(NIGPOINT)
  type(lamina_sdv)    :: ig_sdv_conv, ig_sdv_iter
  real(DP)            :: ig_x(NDIM), ig_u(NDIM), ig_strain(NST), ig_stress(NST)
  integer             :: ig_fstat
  real(DP)            :: ig_df

  !** variables needed for stiffness matrix derivation:
  ! - fn, dn          : shape functions & their derivatives in natural space
  ! - gn              : gradients of shape functions in physical space
  ! - jac, detj       : jacobian matrix and its determinant
  ! - bee             : B matrix
  ! - dee             : D matrix
  ! - btdb            : B' * D * B
  real(DP)            :: fn(NNODE), dn(NNODE,NDIM), gn(NNODE,NDIM)
  real(DP)            :: jac(NDIM,NDIM), detj
  real(DP)            :: bee(NST,NDOF), dee(NST,NST), btdb(NDOF,NDOF)

  !** integer counter variables
  ! - i, j, kig
  integer             :: i, j, kig
  

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
  ply_angle       = ZERO
  elstress        = ZERO
  elstrain        = ZERO
  eldf            = ZERO
  !** nodal variables:
  coords          = ZERO
  u               = ZERO
  !** element characteristic length variables:
  clength         = ZERO
  crosspoints     = ZERO
  !** analysis logical control variables:
  last_converged  = .false.
  nofail          = .false.
  !** integration point variables:
  ig_xi           = ZERO
  ig_weight       = ZERO
  ig_x            = ZERO
  ig_u            = ZERO
  ig_strain       = ZERO
  ig_stress       = ZERO
  ig_fstat        = 0
  ig_df           = ZERO
  !** variables needed for stiffness matrix derivation:
  fn              = ZERO
  dn              = ZERO
  gn              = ZERO
  jac             = ZERO
  detj            = ZERO
  bee             = ZERO
  dee             = ZERO
  btdb            = ZERO
  !** integer counter variables:
  i=0; j=0; kig=0


  ! check validity of input/imported variables

  ! assign values to local variables

  !** local copies of intent(inout) dummy arg.:
  ! elem's fstat, stress, strain and df components are NOT needed
  ! for calculations in this subroutine. here they are only for output,
  ! and they need to be assgined new values during each iteration.
  ! they are therefore NOT passed to their local copies elfstat, elstress,
  ! elstrain and df.
  ! elfstat   : no extraction     ! out
  ! elstress  : no extraction     ! out
  ! elstrain  : no extraction     ! out
  ! eldf      : no extraction     ! out
  ply_angle   = elem%ply_angle    ! in
  local_clock = elem%local_clock  ! inout
  ig_points   = elem%ig_points    ! inout

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
      emsg  = 'x not allocated for node, brickPly_elem_module'
      exit
    end if
    ! assign nodal displacement values (uj) to u vector
    if(allocated(uj)) then
      u((j-1)*NDIM+1:j*NDIM)=uj(1:NDIM)
      deallocate(uj)
    else
      istat = STAT_FAILURE
      emsg  = 'u not allocated for node, brickPly_elem_module'
      exit
    end if
  end do
  ! if there's any error encountered in the extraction process
  ! clean up and exit the program
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

  !** element characteristic length:
  ! this is the length of the line segment orthogonal to the fibre angle,
  ! and passes the element centroid and crosses two edges of the element.
  ! its length is the distance between the cross points on the two edges.
  ! the subroutine below was intended to calculate the cross points between
  ! the matrix crack (passing element centroid) and two element edges. it is
  ! used here to calculate the element characteristic length, by imagining
  ! a crack line orthogonal to the fibre angle, thus: crack_angle=ply_angle+90
  call crack_elem_centroid2d (nedge=NEDGE_BOTTOM, crack_angle=ply_angle+90._DP,&
  & coords=coords(1:2,:), nodes_on_edges=NODES_ON_BOTTOM_EDGES, istat=istat,   &
  & emsg=emsg, edge_crack_points=crosspoints)
  ! if there's any error encountered in the above process
  ! clean up and exit the program
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if
  ! if no error, proceed with the clength calculation
  ! the distance btw two cross points is the characteristic length
  clength = distance (crosspoints(:,2), crosspoints(:,1), 2)

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
  ! update ig point xi and weight
  call init_ig_point (ig_xi, ig_weight)
  ! ig_x, u, strain, stress are updated later in the loop

  !** other local variables are assigned in the following loop of all ig points



  !**** MAIN CALCULATIONS ****


  ! loop over all ig points,
  ! calculate strain, stress, stiffness matrix, sdv etc. at each ig point

  loop_igpoint: do kig = 1, NIGPOINT

      ! get values of shape functions and derivatives, fn and dn, at this ig point
      call init_shape (fn, dn, ig_xi(:,kig))

      ! get jacobian, jac
      jac = matmul(coords,dn)

      ! get determinant of jacobian, detj
      detj = determinant3d (jac)

      ! invert jac onto itself
      call invert_self3d (jac, istat, emsg, detj)
      if (istat == STAT_FAILURE) exit loop_igpoint

      ! calculate gradient of shape function matrix, gn
      gn = matmul(dn,jac)

      ! obtain b matrix (NST*NDOF) from rearranging terms of gn
      bee = beemat3d (gn, NNODE)

      ! calculate global strains at this ig point, ig_strain
      ig_strain = matmul(bee,u)

      ! transfer strain to local coordinates when ply_angle is non-zero
      if(abs(ply_angle) > SMALLNUM) &
      & ig_strain = lcl_strain3d (ig_strain, ply_angle)

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

      ! use material properties, sdv_iter, strain and clength
      ! to calculate D and stress, and update sdv_iter
      if(nofail) then
      ! no failure is allowed, use ddsdde_lamina_intact
        call ddsdde (material, dee=dee, stress=ig_stress, strain=ig_strain)
      else
      ! failure is allowed, use ddsdde_lamina
        call ddsdde (material, dee=dee, stress=ig_stress, sdv=ig_sdv_iter, &
        & strain=ig_strain, clength=clength, istat=istat, emsg=emsg)
        if (istat == STAT_FAILURE) exit loop_igpoint
      end if

      ! get D matrix in global coordinates
      if(abs(ply_angle) > SMALLNUM) dee = glb_dee3d (dee, ply_angle)

      ! calculate B' D B
      btdb = matmul( matmul(transpose(bee),dee), bee )

      ! integrate and update K matrix
      do i=1, NDOF
        do j=1, NDOF
          K_matrix(i,j) = K_matrix(i,j) + btdb(i,j) * detj * ig_weight(kig)
        end do
      end do

      ! update ig point components
      call update(ig_points(kig), converged_sdv=ig_sdv_conv, &
      & iterating_sdv=ig_sdv_iter)

      ! update elem fstat to be the max current ig point fstat
      call extract (ig_sdv_iter, fstat=ig_fstat)
      elfstat  = max(elfstat, ig_fstat)

      ! update elem df to be the max current ig point df
      call extract (ig_sdv_iter, df=ig_df)
      eldf     = max(eldf,    ig_df)

      ! update elem stress & strain (avg of ig point stress & strains)
      elstress = elstress + ig_stress/real(NIGPOINT, DP)
      elstrain = elstrain + ig_strain/real(NIGPOINT, DP)

      ! empty relevant arrays for reuse
      fn  = ZERO
      dn  = ZERO
      gn  = ZERO
      jac = ZERO
      detj = ZERO
      bee = ZERO
      dee  = ZERO
      btdb = ZERO
      ig_x = ZERO
      ig_u = ZERO
      ig_strain=ZERO
      ig_stress=ZERO
      ig_fstat = 0
      ig_df    =ZERO

  end do loop_igpoint !-looped over all int points. ig=NIGPOINT

  ! check to see if the loop is exited upon error
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

  !**** END MAIN CALCULATIONS ****


  ! check to see if intent(inout) derived type dummy arg's components are
  ! unintentionally modified
  if (abs(ply_angle - elem%ply_angle) > SMALLNUM) then
    istat = STAT_FAILURE
    emsg  = 'elem%ply_angle is unintentionally modified in brickPly_elem module'
  end if
  if (istat == STAT_FAILURE) then
    call clean_up (K_matrix, F_vector, uj, xj)
    return
  end if

  ! if no unintentinal modification, proceed with final calculations
  ! and updates and exit the program

  ! calculate F_vector
  F_vector = matmul(K_matrix,u)

  ! update intent(inout) dummy arg./its components
  elem%fstat       = elfstat
  elem%df          = eldf
  elem%stress      = elstress
  elem%strain      = elstrain
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

end subroutine integrate_brickPly_elem












! the rest are private subroutines









pure subroutine init_ig_point (xi, wt)
! used parameters:
! ZERO, EIGHT, ONE_ROOT3, ONE

  real(DP), intent(inout) :: xi(NDIM,NIGPOINT), wt(NIGPOINT)

    if (NIGPOINT .eq. 1) then
      xi(1,1)= ZERO
      xi(2,1)= ZERO
      xi(3,1)= ZERO
      wt = EIGHT
    else if (NIGPOINT .eq. 8) then
      xi(1,1)= -ONE_ROOT3
      xi(2,1)= -ONE_ROOT3
      xi(3,1)= -ONE_ROOT3
      xi(1,2)=  ONE_ROOT3
      xi(2,2)= -ONE_ROOT3
      xi(3,2)= -ONE_ROOT3
      xi(1,3)= -ONE_ROOT3
      xi(2,3)=  ONE_ROOT3
      xi(3,3)= -ONE_ROOT3
      xi(1,4)=  ONE_ROOT3
      xi(2,4)=  ONE_ROOT3
      xi(3,4)= -ONE_ROOT3
      xi(1,5)= -ONE_ROOT3
      xi(2,5)= -ONE_ROOT3
      xi(3,5)=  ONE_ROOT3
      xi(1,6)=  ONE_ROOT3
      xi(2,6)= -ONE_ROOT3
      xi(3,6)=  ONE_ROOT3
      xi(1,7)= -ONE_ROOT3
      xi(2,7)=  ONE_ROOT3
      xi(3,7)=  ONE_ROOT3
      xi(1,8)=  ONE_ROOT3
      xi(2,8)=  ONE_ROOT3
      xi(3,8)=  ONE_ROOT3
      wt = ONE
    end if

end subroutine init_ig_point



pure subroutine init_shape (f, df, ig_xi)
! used parameters:
! ZERO, ONE_EIGHTH, ONE

    real(DP), intent(inout)  :: f(NNODE),df(NNODE,NDIM)
    real(DP), intent(in)     :: ig_xi(NDIM)

    ! local variables
    real(DP) :: xi, eta, zeta 
    
    xi   = ZERO
    eta  = ZERO
    zeta = ZERO

    xi   = ig_xi(1)
    eta  = ig_xi(2)
    zeta = ig_xi(3)

    f(1)=ONE_EIGHTH*(ONE-xi)*(ONE-eta)*(ONE-zeta)
    f(2)=ONE_EIGHTH*(ONE+xi)*(ONE-eta)*(ONE-zeta)
    f(3)=ONE_EIGHTH*(ONE+xi)*(ONE+eta)*(ONE-zeta)
    f(4)=ONE_EIGHTH*(ONE-xi)*(ONE+eta)*(ONE-zeta)
    f(5)=ONE_EIGHTH*(ONE-xi)*(ONE-eta)*(ONE+zeta)
    f(6)=ONE_EIGHTH*(ONE+xi)*(ONE-eta)*(ONE+zeta)
    f(7)=ONE_EIGHTH*(ONE+xi)*(ONE+eta)*(ONE+zeta)
    f(8)=ONE_EIGHTH*(ONE-xi)*(ONE+eta)*(ONE+zeta)

    df(1,1) = -ONE_EIGHTH*(ONE-eta)*(ONE-zeta)
    df(2,1) =  ONE_EIGHTH*(ONE-eta)*(ONE-zeta)
    df(3,1) =  ONE_EIGHTH*(ONE+eta)*(ONE-zeta)
    df(4,1) = -ONE_EIGHTH*(ONE+eta)*(ONE-zeta)
    df(5,1) = -ONE_EIGHTH*(ONE-eta)*(ONE+zeta)
    df(6,1) =  ONE_EIGHTH*(ONE-eta)*(ONE+zeta)
    df(7,1) =  ONE_EIGHTH*(ONE+eta)*(ONE+zeta)
    df(8,1) = -ONE_EIGHTH*(ONE+eta)*(ONE+zeta)

    df(1,2) = -ONE_EIGHTH*(ONE-xi)*(ONE-zeta)
    df(2,2) = -ONE_EIGHTH*(ONE+xi)*(ONE-zeta)
    df(3,2) =  ONE_EIGHTH*(ONE+xi)*(ONE-zeta)
    df(4,2) =  ONE_EIGHTH*(ONE-xi)*(ONE-zeta)
    df(5,2) = -ONE_EIGHTH*(ONE-xi)*(ONE+zeta)
    df(6,2) = -ONE_EIGHTH*(ONE+xi)*(ONE+zeta)
    df(7,2) =  ONE_EIGHTH*(ONE+xi)*(ONE+zeta)
    df(8,2) =  ONE_EIGHTH*(ONE-xi)*(ONE+zeta)

    df(1,3) = -ONE_EIGHTH*(ONE-xi)*(ONE-eta)
    df(2,3) = -ONE_EIGHTH*(ONE+xi)*(ONE-eta)
    df(3,3) = -ONE_EIGHTH*(ONE+xi)*(ONE+eta)
    df(4,3) = -ONE_EIGHTH*(ONE-xi)*(ONE+eta)
    df(5,3) =  ONE_EIGHTH*(ONE-xi)*(ONE-eta)
    df(6,3) =  ONE_EIGHTH*(ONE+xi)*(ONE-eta)
    df(7,3) =  ONE_EIGHTH*(ONE+xi)*(ONE+eta)
    df(8,3) =  ONE_EIGHTH*(ONE-xi)*(ONE+eta)


end subroutine init_shape




end module brickPly_elem_module
