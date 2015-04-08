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
use parameter_module, only : DP, SDV_ARRAY, ZERO,
use integration_point_module

implicit none
private

! list of private parameters in this module:
!
! NDIM          : dimension of this element
! NSTRAIN       : no. of strains
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
integer, parameter :: NDIM=3, NSTRAIN=6, NIGPOINT=8, NNODE=8, NEDGE_BOTTOM=4, &
                    & NDOF=NDIM*NNODE
integer, parameter :: NODE_ON_BOTTOM_EDGE(2, NEDGE_BOTTOM) = &
                    & reshape([ 1,2, 2,3, 3,4, 4,1 ], [2, NEDGE_BOTTOM])


type, public :: brick_element 
  private
  ! list of type components:
  ! curr_status : element current status
  ! ID_elem   : index of this element in the global element array
  ! connec    : indices of element nodes in the global node array
  ! ID_matkey : index of element material in the global matkey array
  ! ply_angle : ply angle for composite lamina (rotation around z axis)
  ! stress    : stresses for output
  ! strain    : strains for output
  ! ig_point  : x, xi, weight, stress, strain, sdv; initialize in set procedure
  ! sdv       : element solution dependent variables
  integer  :: curr_status     = 0
  integer  :: ID_elem         = 0
  integer  :: connec(NNODE)   = 0
  integer  :: ID_matkey       = 0
  real(DP) :: ply_angle       = ZERO
  real(DP) :: stress(NSTRAIN) = ZERO
  real(DP) :: strain(NSTRAIN) = ZERO
  
  type(integration_point) :: ig_point(NIGPOINT)
  
  type(SDV_ARRAY), allocatable :: sdv(:)
    
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
  
  integer :: i
  i = 0
  
  elem%curr_status = 0
  elem%ID_elem     = 0
  elem%connec      = 0
  elem%ID_matkey   = 0
  elem%ply_angle   = ZERO        
  elem%stress      = ZERO
  elem%strain      = ZERO
 
  do i=1, NIGPOINT
    call empty (elem%ig_point(i))
  end do      

  if(allocated(elem%sdv)) deallocate(elem%sdv) 

end subroutine empty_brick_element



pure subroutine set_brick_element (elem, ID_elem, connec, ID_matkey, ply_angle)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module

  type(brick_element),    intent(inout)   :: elem
  integer,                intent(in)      :: connec(NNODE)
  integer,                intent(in)      :: ID_elem,ID_matkey
  real(DP),               intent(in)      :: ply_angle

  ! local variables
  real(DP) :: x(NDIM), u(NDIM), stress(NSTRAIN), strain(NSTRAIN)
  integer  :: i
  
  ! initialize local variables
  x      = ZERO
  u      = ZERO
  stress = ZERO
  strain = ZERO
  i      = 0
  
  elem%ID_elem   = ID_elem
  elem%connec    = connec
  elem%ID_matkey = ID_matkey
  elem%ply_angle = ply_angle
  
  do i=1, NIGPOINT
    call update (elem%ig_point(i), x=x, u=u, stress=stress, strain=strain)
  end do

end subroutine set_brick_element



pure subroutine extract_brick_element (elem, curr_status, ID_elem, connec, &
& ID_matkey, ply_angle, stress, strain, ig_point, sdv)
! Purpose:
! to extract the components of this element

  type(brick_element),                            intent(in)  :: elem
  integer,                              optional, intent(out) :: curr_status
  integer,                              optional, intent(out) :: ID_elem
  integer,                 allocatable, optional, intent(out) :: connec(:)
  integer,                              optional, intent(out) :: ID_matkey
  real(DP),                             optional, intent(out) :: ply_angle
  real(DP),                allocatable, optional, intent(out) :: stress(:)
  real(DP),                allocatable, optional, intent(out) :: strain(:)
  type(integration_point), allocatable, optional, intent(out) :: ig_point(:)
  type(SDV_ARRAY),         allocatable, optional, intent(out) :: sdv(:)
  
  if (present(curr_status)) curr_status = elem%curr_status
  
  if (present(ID_elem))     ID_elem     = elem%ID_elem
  
  if (present(ID_matkey))   ID_matkey   = elem%ID_matkey
  
  if (present(ply_angle))   ply_angle   = elem%ply_angle
  
  if (present(stress)) then
    allocate(stress(NSTRAIN))
    stress = elem%stress
  end if

  if (present(strain)) then
    allocate(strain(NSTRAIN))
    strain = elem%strain
  end if
  
  if (present(connec)) then
    allocate(connec(NNODE))
    connec = elem%connec
  end if
  
  if (present(ig_point)) then
    allocate(ig_point(NIGPOINT))
    ig_point = elem%ig_point
  end if
  
  if (present(sdv)) then        
    if (allocated(elem%sdv)) then
      allocate(sdv(size(elem%sdv)))
      sdv = elem%sdv
    end if
  end if

end subroutine extract_brick_element



pure subroutine integrate_brick_element (elem, K_matrix, F_vector, nofailure)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element

use toolkit_module     ! global tools for element integration
use lib_mat_module     ! global material library
use lib_node_module    ! global node library
use glb_clock_module   ! global analysis progress (curr. step, inc, time, dtime)

  type(brick_element),   intent(inout) :: elem 
  real(DP), allocatable, intent(out)   :: K_matrix(:,:), F_vector(:)
  logical,  optional,    intent(in)    :: nofailure
  
  ! the rest are all local variables
  
  ! variables to be extracted from global arrays
  type(xnode)                   :: node(NNODE) ! x, u, du, v, extra dof ddof etc
  type(global_matkey)           :: matkey ! matname, mattype and ID_matkey to glb mattype array
  character(len=matnamelength)  :: matname
  character(len=mattypelength)  :: mattype
  integer                       :: type_array_index
  ! - glb clock step and increment no. extracted from glb clock module
  integer   :: curr_step, curr_inc
  ! - variables extracted from element isdv
  integer   :: nstep, ninc                  ! step and increment no. of the last iteration, stored in the element
  logical   :: last_converged               ! true if last iteration has converged: a new increment/step has started
  ! - variables extracted from intg point sdvs
  type(SDV_ARRAY), allocatable  :: ig_sdv(:)
  
  ! variables defined locally
  real(DP) :: igxi(NDIM,NIGPOINT),igwt(NIGPOINT) ! ig point natural coords and weights
  real(DP) :: coords(NDIM,NNODE) ! coordinates of the element nodes
  real(DP) :: theta ! orientation of element local coordinates
  real(DP) :: u(NDOF) ! nodal disp. vector
  logical  :: failure ! true for failure analysis
  real(DP) :: fn(NNODE),dn(NNODE,NDIM) ! shape functions & their deriv. physical space
  real(DP) :: jac(NDIM,NDIM),gn(NNODE,NDIM),detj ! jacobian & shape func. deriv. natural space
  real(DP) :: bee(NSTRAIN,NDOF),beet(NDOF,NSTRAIN) ! b matrix and its transpose
  real(DP) :: dee(NSTRAIN,NSTRAIN) ! d matrix
  real(DP) :: btd(NDOF,NSTRAIN),btdb(NDOF,NDOF) ! b'*d & b'*d*b
  real(DP) :: tmpx(NDIM),tmpu(NDIM),tmpstrain(NSTRAIN),tmpstress(NSTRAIN) ! temp. x, strain & stress arrays for intg pnts      
  real(DP), allocatable :: xj(:),uj(:)! nodal vars extracted from glb lib_node array
  
  integer :: i,j,kig,igstat
  
  ! variables used for calculating clength
  real(DP)   :: clength,ctip(2,2)
  
  logical :: nofail
  
  ! initialize variables
  allocate(K_matrix(NDOF,NDOF),F_vector(NDOF))
  K_matrix=ZERO; F_vector=ZERO
  
  i=0; j=0; kig=0
  do i=1,NNODE
      call empty(node(i))
  end do
  call empty(mat)
  curr_step=0; curr_inc=0
  
  igxi=ZERO; igwt=ZERO
  coords=ZERO; theta=ZERO; u=ZERO
  failure=.false.
  
  nofail=.false.
  if(present(nofailure)) nofail=nofailure
  
  ! copy nodes from global node array 
  node(:)=lib_node(elem%connec(:))
  
  ! assign values to local arrays from nodal values
  do j=1,NNODE
  
    ! extract useful values from nodes
    call extract(node(j),x=xj,u=uj)
    
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
  mat=lib_mat(elem%ID_matkey)
  call extract(mat,matname,mattype,typekey)
  
  ! extract ply angle from element definition
  theta=elem%ply_angle
  
  ! - extract curr step and inc values from glb clock module
  call extract_glb_clock(kstep=curr_step,kinc=curr_inc)
  
  ! - check if last iteration has converged, and update the current step & increment no.
  if(elem%nstep.ne.curr_step .or. elem%ninc.ne.curr_inc) then
    last_converged=.true.
    elem%nstep = curr_step
    elem%ninc = curr_inc
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
  call init_ig(igxi,igwt)
  
  ! ZERO elem curr status for update
  elem%curr_status=ZERO
    
  ! ZERO element stress and strain (used for output only) for update
  elem%stress=ZERO
  elem%strain=ZERO
  
  !-calculate strain,stress,stiffness,sdv etc. at each int point
  do kig=1,NIGPOINT 
  
    ! empty relevant arrays for reuse
    fn=ZERO; dn=ZERO
    jac=ZERO; gn=ZERO; detj=ZERO
    bee=ZERO; beet=ZERO; dee=ZERO
    btd=ZERO; btdb=ZERO
    tmpx=ZERO; tmpu=ZERO; tmpstrain=ZERO; tmpstress=ZERO
    
    !- get shape matrix and derivatives
    call init_shape(igxi(:,kig),fn,dn) 
    
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
    
    !-obtain b matrix (NSTRAIN*NDOF) from rearranging terms of gn
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
    
        ! check if failure analysis is needed (check if strength parameters are present)
        call extract(lib_lamina(typekey),strength_active=failure)
        
        if(nofail) failure=.false.
        
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
        K_matrix(i,j) = K_matrix(i,j)+btdb(i,j)*detj*igwt(kig) !-gauss integration
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



pure subroutine init_shape (igxi, f, df)
  
    real(DP),intent(inout) :: f(NNODE),df(NNODE,NDIM)
    real(DP),intent(in) :: igxi(NDIM)
    
    real(DP) :: xi,eta,zeta ! local variables
    xi=ZERO
    eta=ZERO
    zeta=ZERO
    

    xi=igxi(1)
    eta=igxi(2)
    zeta=igxi(3)
    
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
