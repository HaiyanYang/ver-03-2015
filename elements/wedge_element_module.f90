module wedge_element_module
!
!  Purpose:
!    define a wedge element object
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
!             | E3/     E2\   | 
!             |  /         \  |
!             | /           \ |
!             |/_____________\|
!             1      E1       2
!
!  bottom surface nodes (counter-clock vise): 1, 2, 3
!  top    surface nodes (counter-clock vise): 4, 5, 6
!  bottom surface edges (counter-clock vise): E1, E2, E3
!  end nodes of E1 : 1, 2
!  end nodes of E2 : 2, 3
!  end nodes of E3 : 3, 1
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    08/04/15  B. Y. Chen            Original code
!

use parameter_module, only : NDIM, NST => NST_STANDARD, DP, ZERO,
! list of external modules used in type definition and other procedures:
! global clock module    : needed in element definition, extract and integrate
! lamina material module : needed in element definition, extract and integrate 
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
integer, parameter :: NIGPOINT=6, NNODE=6, NEDGE_BOTTOM=3, NDOF=NDIM*NNODE
integer, parameter :: NODE_ON_BOTTOM_EDGE(2, NEDGE_BOTTOM) = &
                    & reshape([ 1,2, 2,3, 3,1 ], [2, NEDGE_BOTTOM])

type, public :: wedge_element 
  private
  ! list of type components:
  ! fstat         : element failure status
  ! connec        : indices of element nodes in the global node list
  ! ID_matlist    : index of element material in its material list
  ! ply_angle     : ply angle for composite lamina (rotation around z axis)
  ! ig_points     : integration points of this element
  ! local_clock   : locally-saved program clock
  ! stress        : stresses for output
  ! strain        : strains for output  
  ! df            : fibre degradation factor for output
  integer  :: fstat         = 0
  integer  :: connec(NNODE) = 0
  integer  :: ID_matlist    = 0
  real(DP) :: ply_angle     = ZERO
  type(program_clock)   :: local_clock
  type(lamina_ig_point) :: ig_points(NIGPOINT)
  real(DP) :: stress(NST)   = ZERO
  real(DP) :: strain(NST)   = ZERO
  real(DP) :: df            = ZERO
end type




interface empty
  module procedure empty_wedge_element
end interface

interface set
  module procedure set_wedge_element
end interface

interface integrate
  module procedure integrate_wedge_element
end interface

interface extract
  module procedure extract_wedge_element
end interface




public :: empty, set, integrate, extract




contains




pure subroutine empty_wedge_element (elem)
! Purpose:
! this subroutine is used to format the element for use
! it is used in the initialize_lib_elem procedure in the lib_elem module

  type(wedge_element), intent(inout) :: elem
  
  ! local variable, derived type var. is initialized upon declaration
  type(wedge_element) :: elem_lcl
  
  ! reset elem to the initial state
  elem = elem_lcl

end subroutine empty_wedge_element



pure subroutine set_wedge_element (elem, connec, ID_matlist, ply_angle)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
! note that only some of the components need to be set during preproc,
! namely connec, ID_matlist, ply_angle and local_clock

  type(wedge_element),    intent(inout)   :: elem
  integer,                intent(in)      :: connec(NNODE)
  integer,                intent(in)      :: ID_matlist
  real(DP),               intent(in)      :: ply_angle
  
  elem%connec    = connec
  elem%ID_matlist = ID_matlist
  elem%ply_angle = ply_angle

end subroutine set_wedge_element



pure subroutine extract_wedge_element (elem, fstat, connec, ID_matlist, &
& ply_angle, local_clock, ig_points, stress, strain, df)
! Purpose:
! to extract the components of this element
! note that the dummy args connec and ig_points are allocatable arrays
! because their sizes vary with different element types

  type(wedge_element),                          intent(in)  :: elem  
  integer,                            optional, intent(out) :: fstat
  integer,               allocatable, optional, intent(out) :: connec(:)
  integer,                            optional, intent(out) :: ID_matlist
  real(DP),                           optional, intent(out) :: ply_angle
  type(program_clock),                optional, intent(out) :: local_clock
  type(lamina_ig_point), allocatable, optional, intent(out) :: ig_points(:)
  real(DP),                           optional, intent(out) :: stress(NST)
  real(DP),                           optional, intent(out) :: strain(NST)
  real(DP),                           optional, intent(out) :: df
  
  if (present(fstat))       fstat = elem%fstat
  
  if (present(connec)) then
    allocate(connec(NNODE))
    connec = elem%connec
  end if
  
  if (present(ID_matlist))  ID_matlist  = elem%ID_matlist
  
  if (present(ply_angle))   ply_angle   = elem%ply_angle

  if (present(local_clock)) local_clock = elem%local_clock
  
  if (present(ig_points)) then
    allocate(ig_points(NIGPOINT))
    ig_points = elem%ig_points
  end if
  
  if (present(stress))      stress = elem%stress

  if (present(strain))      strain = elem%strain
  
  if (present(df))          df     = elem%df

end subroutine extract_wedge_element



pure subroutine integrate_wedge_element (elem, K_matrix, F_vector, nofailure)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element

use toolkit_module                  ! global tools for element integration
use lib_mat_module                  ! global material library
use lib_node_module                 ! global node library

  type(wedge_element),intent(inout)       :: elem 
  real(kind=DP),allocatable,intent(out)   :: K_matrix(:,:), F_vector(:)
  logical,  optional,  intent(in)         :: nofailure
  
  ! the rest are all local variables
  
  ! variables to be extracted from global arrays
  type(xnode) :: node(NNODE) ! x, u, du, v, extra dof ddof etc
  type(material) :: mat ! matname, mattype and ID_matlist to glb mattype array
  character(len=matnamelength) :: matname
  character(len=mattypelength) :: mattype
  integer :: typekey
  ! - glb clock step and increment no. extracted from glb clock module
  integer         :: curr_step, curr_inc
  
  ! - variables extracted from element isdv
  integer         :: nstep, ninc                  ! step and increment no. of the last iteration, stored in the element
  logical         :: last_converged               ! true if last iteration has converged: a new increment/step has started
  
  ! - variables extracted from intg point sdvs
  type(sdv_array),  allocatable   :: ig_sdv(:)
  
  ! variables defined locally
  real(kind=DP)   :: igxi(NDIM,NIGPOINT),igwt(NIGPOINT) ! ig point natural coords and weights
  real(kind=DP)   :: coords(NDIM,NNODE) ! coordinates of the element nodes
  real(kind=DP)   :: theta ! orientation of element local coordinates
  real(kind=DP)   :: u(NDOF) ! nodal disp. vector
  logical         :: failure ! true for failure analysis
  real(kind=DP)   :: fn(NNODE),dn(NNODE,NDIM) ! shape functions & their deriv. physical space
  real(kind=DP)   :: jac(NDIM,NDIM),gn(NNODE,NDIM),detj ! jacobian & shape func. deriv. natural space
  real(kind=DP)   :: bee(NST,NDOF),beet(NDOF,NST) ! b matrix and its transpose
  real(kind=DP)   :: dee(NST,NST) ! d matrix
  real(kind=DP)   :: btd(NDOF,NST),btdb(NDOF,NDOF) ! b'*d & b'*d*b
  real(kind=DP)   :: tmpx(NDIM),tmpu(NDIM),tmpstrain(NST),tmpstress(NST) ! temp. x, strain & stress arrays for intg pnts      
  real(kind=DP),allocatable :: xj(:),uj(:)! nodal vars extracted from glb lib_node array
  
  integer :: i,j,kig,igstat
  
  ! variables used for calculating clength
  real(kind=DP)   :: clength,ctip(2,2)

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
  mat=lib_mat(elem%ID_matlist)
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
  call elem_ctips_origin('wedge',theta,coords,NODE_ON_BOTTOM_EDGE,NEDGE_BOTTOM,ctip)
  clength=sqrt((ctip(1,2)-ctip(1,1))**2+(ctip(2,2)-ctip(2,1))**2)
  
  !-----------------------------------------------------------!
  !-----------------------------------------------------------! 
  
  
  ! update ig point xi and weight
  call init_ig(igxi,igwt)
  
  ! ZERO elem curr status for update
  elem%fstat=ZERO
    
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
      
      !-obtain b matrix (NST*NDOF) from rearranging terms of gn
      bee=beemat(gn)
      
      ! calculate global strains
      tmpstrain=matmul(bee,u)
      
      ! transfer strain to local coordinates
      if(theta/=ZERO) tmpstrain=lcl_strain(tmpstrain,theta)
      
      
      ! - extract sdvs from integration points; ig_sdv automatically deallocated when passed in
      call extract(elem%ig_points(kig),sdv=ig_sdv)
      
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
          case ('isotropic')
              
              ! check if failure analysis is needed (check if strength parameters are present)
              call extract(lib_iso(typekey),strength_active=failure)
              
              if(nofail) failure=.false.
              
              if(failure) write(msg_file,*) "WARNING: failure analysis is not yet supported for &
              & wedge_element type isotropic material; only linear elastic stiffness matrix is integrated."
              
              ! calculate D matrix, update tmpstress
              call ddsdde(lib_iso(typekey),dee,strain=tmpstrain,stress=tmpstress) 
              
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
      !call update(elem%ig_points(kig),x=tmpx,u=tmpu,strain=tmpstrain,stress=tmpstress,sdv=ig_sdv)
      call update(elem%ig_points(kig),sdv=ig_sdv)
      
      ! update elem curr status variable
      igstat=0
      if(allocated(ig_sdv(2)%i)) igstat=ig_sdv(2)%i(1)
      elem%fstat=max(elem%fstat,igstat)

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
  
    

end subroutine integrate_wedge_element












! the rest are private subroutines









pure subroutine init_ig (xi, wt)

  real(kind=DP),intent(inout) :: xi(NDIM,NIGPOINT),wt(NIGPOINT)

    if (NIGPOINT .eq. 2) then
        xi(1,1)= one_third
        xi(2,1)= one_third
        xi(3,1)= -root3
        xi(1,2)= one_third
        xi(2,2)= one_third
        xi(3,2)= root3
        wt = half
    else if(NIGPOINT == 6) then
        xi(1,1)=one/six
        xi(2,1)=one/six
        xi(3,1)=-root3

        xi(1,2)=four/six
        xi(2,2)=one/six
        xi(3,2)=-root3

        xi(1,3)=one/six
        xi(2,3)=four/six
        xi(3,3)=-root3

        xi(1,4)=one/six
        xi(2,4)=one/six
        xi(3,4)=root3

        xi(1,5)=four/six
        xi(2,5)=one/six
        xi(3,5)=root3

        xi(1,6)=one/six
        xi(2,6)=four/six
        xi(3,6)=root3

        wt=one/six

!~        else if(NIGPOINT == 9) then
!~            xi(1,1)=one/six
!~            xi(2,1)=one/six
!~            xi(3,1)=-sqrt(three/five)
!~
!~            xi(1,2)=four/six
!~            xi(2,2)=one/six
!~            xi(3,2)=-sqrt(three/five)
!~
!~            xi(1,3)=one/six
!~            xi(2,3)=four/six
!~            xi(3,3)=-sqrt(three/five)
!~
!~            xi(1,4)=one/six
!~            xi(2,4)=one/six
!~            xi(3,4)=sqrt(three/five)
!~
!~            xi(1,5)=four/six
!~            xi(2,5)=one/six
!~            xi(3,5)=sqrt(three/five)
!~
!~            xi(1,6)=one/six
!~            xi(2,6)=four/six
!~            xi(3,6)=sqrt(three/five)
!~
!~            xi(1,7)=one/six
!~            xi(2,7)=one/six
!~            xi(3,7)=ZERO
!~
!~            xi(1,8)=four/six
!~            xi(2,8)=one/six
!~            xi(3,8)=ZERO
!~
!~            xi(1,9)=one/six
!~            xi(2,9)=four/six
!~            xi(3,9)=ZERO
!~
!~            wt(1:6)=five/54._dp
!~            wt(7:9)=eight/54._dp

    else
        write(msg_file,*) 'no. of integration points incorrect for wedge_ig!'
        call exit_function
    end if

end subroutine init_ig



pure subroutine init_shape (igxi, f, df)
  
    real(kind=DP),intent(inout) :: f(NNODE),df(NNODE,NDIM)
    real(kind=DP),intent(in) :: igxi(NDIM)
    
    real(kind=DP) :: xi,eta,zeta ! local variables
    xi=ZERO
    eta=ZERO
    zeta=ZERO

    xi=igxi(1)
    eta=igxi(2)
    zeta=igxi(3)
    
    f(1)=half*(one-xi-eta)*(one-zeta)
    f(2)=half*xi*(one-zeta)
    f(3)=half*eta*(one-zeta)
    f(4)=half*(one-xi-eta)*(one+zeta)
    f(5)=half*xi*(one+zeta)
    f(6)=half*eta*(one+zeta)
    
    
    df(1,3) = -half*(one-xi-eta)
    df(2,3) = -half*xi
    df(3,3) = -half*eta
    df(4,3) =  half*(one-xi-eta)
    df(5,3) =  half*xi
    df(6,3) =  half*eta
    
    
    df(1,1) = -half*(one-zeta)
    df(2,1) =  half*(one-zeta)
    df(3,1) =  ZERO
    df(4,1) = -half*(one+zeta)
    df(5,1) =  half*(one+zeta)
    df(6,1) =  ZERO
    
    df(1,2) = -half*(one-zeta)
    df(2,2) =  ZERO
    df(3,2) =  half*(one-zeta)
    df(4,2) = -half*(one+zeta)
    df(5,2) =  ZERO
    df(6,2) =  half*(one+zeta)
    
    

end subroutine init_shape




end module wedge_element_module
