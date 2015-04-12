module brick_element_module
!
!  Purpose:
!    define a brick element object
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
  ! curr_status : element current status
  ! connec    : indices of element nodes in the global node array
  ! ID_matkey : index of element material in the global matkey array
  ! ply_angle : ply angle for composite lamina (rotation around z axis)
  ! stress    : stresses for output
  ! strain    : strains for output
  ! ig_point  : x, xi, weight, stress, strain, sdv; initialize in set procedure
  ! local_clock     : locally-saved program clock
  ! equilibrium_sdv : element lamina_sdv at each incremental equilibrium 
  ! iterating_sdv   : element lamina_sdv during iterations of increments
  integer  :: curr_status   = 0
  integer  :: connec(NNODE) = 0
  integer  :: ID_matkey     = 0
  real(DP) :: ply_angle     = ZERO
  real(DP) :: stress(NST)   = ZERO
  real(DP) :: strain(NST)   = ZERO
  type(program_clock)   :: local_clock
  type(lamina_ig_point) :: ig_point(NIGPOINT)
  type(lamina_sdv)      :: equilibrium_sdv
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



pure subroutine set_brick_element (elem, connec, ID_matkey, ply_angle)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
! note that only some of the components need to be set during preproc,
! namely connec, ID_matkey, ply_angle and local_clock

  type(brick_element),    intent(inout)   :: elem
  integer,                intent(in)      :: connec(NNODE)
  integer,                intent(in)      :: ID_matkey
  real(DP),               intent(in)      :: ply_angle
  
  elem%connec    = connec
  elem%ID_matkey = ID_matkey
  elem%ply_angle = ply_angle
  elem%local_clock = GLOBAL_CLOCK

end subroutine set_brick_element



pure subroutine extract_brick_element (elem, curr_status, connec, &
& ID_matkey, ply_angle, stress, strain, local_clock, ig_point,    &
& equilibrium_sdv, iterating_sdv)
! Purpose:
! to extract the components of this element
! note that the dummy args connec and ig_point are allocatable arrays
! because their sizes vary with different element types

  type(brick_element),                          intent(in)  :: elem  
  integer,                            optional, intent(out) :: curr_status
  integer,               allocatable, optional, intent(out) :: connec(:)
  integer,                            optional, intent(out) :: ID_matkey
  real(DP),                           optional, intent(out) :: ply_angle
  real(DP),                           optional, intent(out) :: stress(NST)
  real(DP),                           optional, intent(out) :: strain(NST)
  type(program_clock),                optional, intent(out) :: local_clock
  type(lamina_ig_point), allocatable, optional, intent(out) :: ig_point(:)
  type(lamina_sdv),                   optional, intent(out) :: equilibrium_sdv
  type(lamina_sdv),                   optional, intent(out) :: iterating_sdv
  
  if (present(curr_status)) curr_status = elem%curr_status
  
  if (present(connec)) then
    allocate(connec(NNODE))
    connec = elem%connec
  end if
  
  if (present(ID_matkey))   ID_matkey   = elem%ID_matkey
  
  if (present(ply_angle))   ply_angle   = elem%ply_angle
  
  if (present(stress)) stress = elem%stress

  if (present(strain)) strain = elem%strain

  if (present(local_clock)) local_clock = elem%local_clock
  
  if (present(ig_point)) then
    allocate(ig_point(NIGPOINT))
    ig_point = elem%ig_point
  end if
  
  if (present(equilibrium_sdv))  equilibrium_sdv = elem%equilibrium_sdv
  
  if (present(iterating_sdv))    iterating_sdv   = elem%iterating_sdv

end subroutine extract_brick_element



pure subroutine integrate_brick_element (elem, K_matrix, F_vector, istat, emsg,&
& nofailure)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element
! note that K and F are allocatable dummy args as their sizes vary with
! different element types

use toolkit_module     ! global tools for element integration
use lib_mat_module     ! global material library
use lib_node_module    ! global node library
use glb_clock_module   ! global analysis progress (curr. step, inc, time, dtime)

  type(brick_element),      intent(inout) :: elem 
  real(DP),    allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  integer,                  intent(out)   :: istat
  character(len=MSGLENGTH), intent(out)   :: emsg
  logical,        optional, intent(in)    :: nofailure
  
  ! the rest are all local variables
  
  ! variables to be extracted from global arrays :
  ! - nodes   : array of element nodes containing x, u, du, v, extra dof/ddof
  ! - matkey  : global_matkey of this element containing matname, mattype and 
  !             type_array_index
  ! - matname : material name
  ! - mattype : material type name
  ! - type_list_index  : the index in the corresponding material type list
  ! - last_converged   : true if last iteration has converged
  type(xnode)                   :: nodes(NNODE)
  type(global_matkey)           :: matkey
  character(len=matnamelength)  :: matname
  character(len=mattypelength)  :: mattype
  integer                       :: type_list_index
  logical                       :: last_converged
  
  ! variables defined locally
  ! - coords        : element nodal coordinates matrix
  ! - u             : element nodal displacemet vector
  ! - theta         : element local orientation
  ! - ig_xi         : integration point natural coordinates
  ! - ig_weight     : integration point weights
  ! - ig_sdv_eq/it  : integration point equilibrium and iterating sdvs

  ! - fn, dn        : shape functions & their derivatives in natural space
  ! - gn            : gradients of shape functions in physical space
  ! - jac, detj     : jacobian matrix and its determinant
  ! - bee, beet     : b matrix and its transpose
  ! - dee           : D matrix
  ! - btd, btdb     : b'*d & b'*d*b
 
  ! - xj, uj        : nodal vars extracted from glb lib_node array
  ! - clength, ctip : variables used for calculating clength
  ! - tmpx, tmpu    : temporary x and u vectors
  ! - tmpstrain, tmpstress : temporary strain and stress vectors
  
  real(DP) :: coords(NDIM,NNODE), u(NDOF), theta
  real(DP) :: ig_xi(NDIM, NIGPOINT), ig_weight(NIGPOINT)
  type(lamina_sdv) :: ig_sdv_eq, ig_sdv_it
  
  real(DP) :: fn(NNODE), dn(NNODE,NDIM), gn(NNODE,NDIM)
  real(DP) :: jac(NDIM,NDIM), detj
  real(DP) :: bee(NST,NDOF), beet(NDOF,NST)
  real(DP) :: dee(NST,NST)
  real(DP) :: btd(NDOF,NST), btdb(NDOF,NDOF)
  
  real(DP), allocatable :: xj(:), uj(:)
  real(DP) :: clength, ctip(2,2)
  real(DP) :: tmpx(NDIM), tmpu(NDIM), tmpstrain(NST), tmpstress(NST)
  
  logical  :: nofail
  integer  :: i, j, kig, igstat
 
  
  ! initialize intent(out) and local variables
  allocate(K_matrix(NDOF,NDOF), F_vector(NDOF))
  K_matrix = ZERO
  F_vector = ZERO
  istat = 0
  emsg  = ''
  
  do i=1, NNODE
      call empty(nodes(i))
  end do
  
  matname = ''
  mattype = ''
  type_list_index = 0
  last_converged = .false.
  
  coords    = ZERO
  u         = ZERO
  theta     = ZERO
  ig_xi     = ZERO
  ig_weight = ZERO

  i=0; j=0; kig=0
  
  nofail=.false.
  if(present(nofailure)) nofail=nofailure
  
  ! copy nodes from global nodes array 
  nodes(:) = lib_node(elem%connec(:))
  
  ! assign values to local arrays from nodal values
  do j=1, NNODE
  
    ! extract useful values from nodes
    call extract(nodes(j), x=xj, u=uj)
    
    ! assign values to coords matrix and u vector
    if(allocated(xj)) then
      coords(:,j)=xj(:)
      deallocate(xj)
    else
      write(msg_file,*)'WARNING: x not allocated for node:',elem%connec(j)
    end if
    
    if(allocated(uj)) then 
      u((j-1)*NDIM+1:j*NDIM)=uj(1:NDIM)
      deallocate(uj)
    else
      write(msg_file,*)'WARNING: u not allocated for node:',elem%connec(j)
    end if
      
  end do
  
  ! extract material values from global material array
  matkey = list_matkey(elem%ID_matkey)
  call extract(matkey, matname, mattype, type_list_index)
  
  ! extract ply angle from element definition
  theta=elem%ply_angle
  
  ! - check if last iteration has converged by checking if the global clock has 
  !   advanced; if so, last iteration is converged and sync the local clock
  if(.not. clock_in_sync(GLOBAL_CLOCK, elem%local_clock)) then
    last_converged=.true.
    elem%local_clock = GLOBAL_CLOCK
  end if
  
  
  
  !-----------------------------------------------------------!
  !           calculate approximate clength
  !-----------------------------------------------------------!
  ! initialize relevant variables
  clength=ZERO
  call elem_ctips_origin('brick',theta,coords,NODE_ON_BOTTOM_EDGE,NEDGE_BOTTOM,ctip)
  clength=sqrt((ctip(1,2)-ctip(1,1))**2+(ctip(2,2)-ctip(2,1))**2)
  
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!
  
  
  ! update ig point xi and weight
  call init_ig(ig_xi,ig_weight)
  
  ! ZERO elem curr status for update
  elem%curr_status=ZERO
    
  ! ZERO element stress and strain (used for output only) for update
  elem%stress=ZERO
  elem%strain=ZERO
  
  !-calculate strain,stress,stiffness,sdv etc. at each int point
  do kig=1, NIGPOINT 
  
    ! empty relevant arrays for reuse
    fn=ZERO; dn=ZERO
    jac=ZERO; gn=ZERO; detj=ZERO
    bee=ZERO; beet=ZERO; dee=ZERO
    btd=ZERO; btdb=ZERO
    tmpx=ZERO; tmpu=ZERO; tmpstrain=ZERO; tmpstress=ZERO
    
    !- get shape matrix and derivatives
    call init_shape(ig_xi(:,kig),fn,dn) 
    
    !- calculate integration point physical coordinates (initial)
    tmpx=matmul(coords,fn)
    
    !- calculate integration point displacement
    do j=1,NDIM
      do i=1,NNODE
        tmpu(j)=tmpu(j)+fn(i)*u((i-1)*NDIM+j)
      end do
    end do
    
    ! get jacobian
    jac=matmul(coords,dn)
    
    !-get determinant of jacobian
    detj=determinant(jac)
    
    ! invert jac onto itself
    jac=inverse(jac,detj)
    
    ! calculate gradient of shape function matrix
    gn=matmul(dn,jac)
    
    !-obtain b matrix (NST*NDOF) from rearranging terms of gn
    bee=beemat(gn)
    
    ! calculate global strains
    tmpstrain=matmul(bee,u)
    
    ! transfer strain to local coordinates
    if(theta/=ZERO) tmpstrain=lcl_strain(tmpstrain,theta)
    
    
    ! - extract sdvs from integration points; ig_sdv automatically deallocated when passed in
    call extract(elem%ig_point(kig),sdv=ig_sdv)
    
    ! allocate ig_sdv arrays for 1st iteration of analysis
    if(.not.allocated(ig_sdv)) then
    ! allocate 2 sets of sdv arrays, 1 for converged sdvs and 1 for iterating sdvs
      allocate(ig_sdv(2))
    end if
    
    ! update converged sdvs (sdv1) with iterating sdvs (sdv2) when last iteration has converged
    ! and revalue iterating sdvs (sdv2) to the last converged sdvs (sdv1) if otherwise
    if(last_converged) then
      ig_sdv(1)=ig_sdv(2)
    else               
      ig_sdv(2)=ig_sdv(1)
    end if
    
    
    ! get D matrix dee accord. to material properties, and update intg point stresses
    select case (mattype)
            
      case ('lamina')
        
        if(nofail) then
        
        else
        
        end if
        
        if(failure) then
            ! calculate D matrix, update tmpstress and sdv
            call ddsdde(lib_lamina(typekey),clength,dee,strain=tmpstrain,stress=tmpstress,sdv=ig_sdv(2),dfail=dfail)
        else
            ! calculate D matrix, update tmpstress
            call ddsdde(lib_lamina(typekey),dee=dee,strain=tmpstrain,stress=tmpstress)
        end if
          
      case default
        write(msg_file,*) 'material type not supported for tri element!'
        call exit_function
            
    end select
    
    ! get D matrix in global coordinates deeg
    if(theta/=ZERO) dee=glb_dee(dee,theta)
    
    ! calculate B' D B
    beet=transpose(bee)
    btd=matmul(beet,dee)
    btdb=matmul(btd,bee)

    ! integrate and update K matrix
    do i=1,NDOF
      do j=1,NDOF
        K_matrix(i,j) = K_matrix(i,j)+btdb(i,j)*detj*ig_weight(kig) !-gauss integration
      end do
    end do	
    
    
    ! update ig point arrays
    !call update(elem%ig_point(kig),x=tmpx,u=tmpu,strain=tmpstrain,stress=tmpstress,sdv=ig_sdv)
    call update(elem%ig_point(kig),sdv=ig_sdv)
    
    ! update elem curr status variable
    igstat=0
    if(allocated(ig_sdv(2)%i)) igstat=ig_sdv(2)%i(1)
    elem%curr_status=max(elem%curr_status,igstat)

    ! update elem stress & strain (avg of ig point stress & strains)
    elem%stress=elem%stress+tmpstress/NIGPOINT
    elem%strain=elem%strain+tmpstrain/NIGPOINT
    
    deallocate(ig_sdv)
      
  end do !-looped over all int points. ig=NIGPOINT
  
  F_vector=matmul(K_matrix,u) 

  
  ! deallocate local dynamic arrays

  if(allocated(xj)) deallocate(xj) 
  if(allocated(uj)) deallocate(uj) 
  if(allocated(ig_sdv)) deallocate(ig_sdv)
    
    

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
