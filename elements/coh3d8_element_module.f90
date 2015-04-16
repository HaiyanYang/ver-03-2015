module coh3d8_element_module
!
!  Purpose:
!    define a coh3d8 element object
!    with the associated procedures to empty, set, integrate and extract
!    its components
!
!  topological definition of this coh3d8 element, local nodal indices:
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
integer, parameter :: NNODE=8, NIGPOINT=4, NDOF=NDIM*NNODE

type, public :: coh3d8_element 
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
  module procedure empty_coh3d8_element
end interface

interface set
  module procedure set_coh3d8_element
end interface

interface integrate
  module procedure integrate_coh3d8_element
end interface

interface extract
  module procedure extract_coh3d8_element
end interface




public :: empty, set, integrate, extract




contains




pure subroutine empty_coh3d8_element (elem)
! Purpose:
! this subroutine is used to format the element for use
! it is used in the initialize_lib_elem procedure in the lib_elem module

  type(coh3d8_element), intent(inout) :: elem
  
  ! local variable, derived type var. is initialized upon declaration
  type(coh3d8_element) :: elem_lcl
  
  ! reset elem to the initial state
  elem = elem_lcl

end subroutine empty_coh3d8_element



pure subroutine set_coh3d8_element (elem, connec, ID_matlist)
! Purpose:
! this subroutine is used to set the components of the element
! it is used in the initialize_lib_elem procedure in the lib_elem module
! note that only some of the components need to be set during preproc,
! namely connec, ID_matlist

  type(coh3d8_element),   intent(inout)   :: elem
  integer,                intent(in)      :: connec(NNODE)
  integer,                intent(in)      :: ID_matlist
  
  elem%connec    = connec
  elem%ID_matlist = ID_matlist

end subroutine set_coh3d8_element



pure subroutine extract_coh3d8_element (elem, fstat, connec, ID_matlist, &
& local_clock, ig_points, traction, separation, dm)
! Purpose:
! to extract the components of this element
! note that the dummy args connec and ig_points are allocatable arrays
! because their sizes vary with different element types

  type(coh3d8_element),                           intent(in)  :: elem
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
   
end subroutine extract_coh3d8_element



pure subroutine integrate_coh3d8_element (elem, K_matrix, F_vector, nofailure, &
& gauss, mnode)
! Purpose:
! updates K matrix, F vector, integration point stress and strain,
! and the solution dependent variables (sdvs) of ig points and element

  use toolkit_module                  ! global tools for element integration
  use lib_mat_module                  ! global material library
  use lib_node_module                 ! global node library

  type(coh3d8_element),       intent(inout)   :: elem 
  real(kind=dp), allocatable, intent(out)     :: K_matrix(:,:), F_vector(:)
  logical,  optional,  intent(in)         :: nofailure
  logical,        optional,   intent(in)      :: gauss
  type(xnode),    optional,   intent(in)      :: mnode(:)  ! material (interpolated) nodes

  !-----------------------------------------!
  ! local variables, extracted from glb libs
  !-----------------------------------------!
  
  ! - element nodes extracted from global node library
  type(xnode)     :: node(NNODE)                  ! x, u, du, v, extra dof ddof etc
  
  ! - element material extracted from global material library
  type(material)  :: mat                          ! matname, mattype and ID_matlist to glb mattype array   

  ! - glb clock step and increment no. extracted from glb clock module
  integer         :: curr_step, curr_inc
  
  logical :: nofail
  
  
  !-----------------------------------------!
  ! - the rest are all pure local variables
  !-----------------------------------------!
  
  ! - variables extracted from element nodes
  real(kind=dp),allocatable :: xj(:),uj(:)        ! nodal vars extracted from glb lib_node array
  real(kind=dp)   :: coords(NDIM,NNODE)           ! coordinates of the element nodes, formed from xj of each node
  real(kind=dp)   :: midcoords(NDIM,NNODE/2)      ! coordinates of the mid-plane
  real(kind=dp)   :: u(NDOF)                      ! element nodal disp. vector, formed from uj of each node
  
  ! - variables extracted from element material
  character(len=matnamelength)    :: matname      ! name of the material assigned to this element
  character(len=mattypelength)    :: mattype      ! type of the material
  integer                         :: typekey       ! index of the material in the material library of its type 
  
  ! - variables extracted from element isdv
  integer         :: nstep, ninc                  ! step and increment no. of the last iteration, stored in the element
  logical         :: last_converged               ! true if last iteration has converged: a new increment/step has started
  
  ! - variables extracted from intg point sdvs
  type(sdv_array),  allocatable   :: ig_sdv(:)
  
  
  ! - variables defined locally
  
  real(kind=dp)   :: igxi(NDIM-1,NIGPOINT),igwt(NIGPOINT)   ! ig point natural coords and weights
  real(kind=dp)   :: tmpx(NDIM), tmpu(NDIM)       ! temporary arrays to hold x and u values for intg points
         
  real(kind=dp)   :: fn(NNODE),dn(NNODE,NDIM-1)   ! shape functions and derivatives
  real(kind=dp)   :: Nmatrix(NDIM,NDOF)           ! obtained from fn, to compute disp. jump acrss intfc: {u}_jump = [N]*{u}
  real(kind=dp)   :: normal(NDIM),tangent1(NDIM),tangent2(NDIM)   ! normal and tangent vectors of the interface, obtained from coords
  real(kind=dp)   :: jac(NDIM-1,NDIM-1), det      ! jacobian and determinant of jacobian
  real(kind=dp)   :: Qmatrix(NDIM,NDIM)           ! rotation matrix from global to local coordinates (from normal & tangent)
  real(kind=dp)   :: ujump(NDIM),delta(NDIM)      ! {delta}=[Q]*{u}_jump, {delta} is the jump vector in lcl coords.
  real(kind=dp)   :: Dee(NDIM,NDIM)               ! material stiffness matrix [D]
  real(kind=dp)   :: Tau(NDIM)                    ! {Tau}=[D]*{delta}, traction on the interface
  
  real(kind=dp)   :: QN(NDIM,NDOF),DQN(NDIM,NDOF)         ! [Q]*[N], [D]*[Q]*[N]
  real(kind=dp)   :: NtQt(NDOF,NDIM),NtQtDQN(NDOF,NDOF)   ! [N']*[Q'], [N']*[Q']*[D]*[Q]*[N] 
  real(kind=dp)   :: NtQtTau(NDOF)                        ! [N']*[Q']*{Tau}
  integer         :: i,j,k,kig,igstat

  
  
  
  !------------------------------------------------!
  !           initialize all variables 
  !------------------------------------------------!
  
  ! intent(out) variables, automatically deallocated when passed in
  allocate(K_matrix(NDOF,NDOF),F_vector(NDOF)); K_matrix=zero; F_vector=zero 
  
  ! integer counters
  i=0; j=0; k=0; kig=0
  
  ! local variables, extracted from glb libs
  do i=1,NNODE
      call empty(node(i))
  end do 
  call empty(mat)
  curr_step=0; curr_inc=0
  
  nofail=.false.
  if(present(nofailure)) nofail=nofailure
  
  ! pure local variables
  coords=zero; midcoords=zero; u=zero 
  matname=''; mattype=''; typekey=0 
  nstep=0; ninc=0; last_converged=.false.
  igxi=zero; igwt=zero
  tmpx=zero; tmpu=zero
  fn=zero; dn=zero; Nmatrix=zero
  tangent1=zero; tangent2=zero; normal=zero 
  jac=zero; det=zero; Qmatrix=zero
  ujump=zero; delta=zero
  Dee=zero; Tau=zero
  QN=zero; DQN=zero; NtQt=zero; NtQtDQN=zero; NtQtTau=zero
  
  
  
  
  
  !------------------------------------------------!
  !   extract variables from global arrays 
  !------------------------------------------------!
  
  if(present(mnode)) then
      ! - extract nodes from passed-in node array
      node(:)=mnode(:)
  else
      ! - extract nodes from global node array 
      node(:)=lib_node(elem%connec(:))
  end if
  
  ! - extract material values from global material array
  mat=lib_mat(elem%ID_matlist)
  
  ! - extract curr step and inc values from glb clock module
  call extract_glb_clock(kstep=curr_step,kinc=curr_inc)
  
  ! - check if last iteration has converged, and update the current step & increment no.
  if(elem%nstep.ne.curr_step .or. elem%ninc.ne.curr_inc) then
      last_converged=.true.
      elem%nstep = curr_step
      elem%ninc = curr_inc
  end if
  
  
  
  !------------------------------------------------!
  !   assign values to material, coords and u 
  !------------------------------------------------!      
  
  ! - extract x and u values from nodes and assign to local arrays 
  do j=1,NNODE
      ! extract x and u values from nodes
      call extract(node(j),x=xj,u=uj)     
      ! assign x to coords matrix
      if(allocated(xj)) then
          coords(:,j)=xj(:)
          deallocate(xj)
      else
          write(msg_file,*)'WARNING: x not allocated for node:',elem%connec(j)
      end if
      ! and u to u vector
      if(allocated(uj)) then 
          u((j-1)*NDIM+1:j*NDIM)=uj(1:NDIM)
          deallocate(uj)
      else
          write(msg_file,*)'WARNING: u not allocated for node:',elem%connec(j)
      end if    
  end do
  
  ! calculate mid-plane coordinates
  do j=1,NNODE/2
      midcoords(:,j)=half*(coords(:,j)+coords(:,j+NNODE/2))
  end do
  
  
  ! - extract values from mat (material type) and assign to local vars (matname, mattype & ID_matlist)
  call extract(mat,matname,mattype,typekey) 
  
  
  
  
  
  !------------------------------------------------!
  !   compute Q matrix (rotation) and determinat 
  !------------------------------------------------!
  
  ! - compute tangent1 of the interface: node 2 coords - node 1 coords
  tangent1(1)=midcoords(1,2)-midcoords(1,1)
  tangent1(2)=midcoords(2,2)-midcoords(2,1)
  tangent1(3)=midcoords(3,2)-midcoords(3,1)
  
  ! - compute tangent2 of the interface: node 4 coords - node 1 coords
  tangent2(1)=midcoords(1,4)-midcoords(1,1)
  tangent2(2)=midcoords(2,4)-midcoords(2,1)
  tangent2(3)=midcoords(3,4)-midcoords(3,1)
  
  ! - compute normal vector of the interface, its magnitude is the area of the triangle-> det
  normal=CrossProduct3D(tangent1,tangent2)
  
  ! - re-evaluate tangent1 so that it is perpendicular to both tangent2 and normal
  tangent1=CrossProduct3D(tangent2,normal)
  
  ! - normalize these vectors
  call normalize(normal)
  call normalize(tangent1)
  call normalize(tangent2)
  
  ! - compute Q matrix
  do j=1,3
      Qmatrix(1,j)=normal(j)
      Qmatrix(2,j)=tangent1(j)
      Qmatrix(3,j)=tangent2(j)
  end do      

  ! - transform midcoords matrix into local isoplanar coordinates
  midcoords=matmul(Qmatrix,midcoords)
  
  
  
  !------------------------------------------------!
  !   perform integration at ig points
  !------------------------------------------------!
  
  ! - calculate ig point xi and weight
  if(present(gauss).and.gauss)  then
      call init_ig(igxi,igwt,gauss)
  else                          
      call init_ig(igxi,igwt)
  end if

  ! zero elem curr status for update
  elem%fstat=zero
   
  !-calculate strain,stress,stiffness,sdv etc. at each int point
  do kig=1,NIGPOINT 
  
      ! - empty relevant arrays for reuse
      fn=zero; dn=zero; Nmatrix=zero
      jac=zero; det=zero
      ujump=zero; delta=zero; dee=zero; Tau=zero
      QN=zero; NtQt=zero; DQN=zero; NtQtDQN=zero; NtQtTau=zero
      tmpx=zero; tmpu=zero
      
      !- get shape function matrix and their derivatives (used to calculate jacobian)
      call init_shape(igxi(:,kig),fn,dn)

      
      ! - calculate jacobian; midcoords 1st row are now zero (out-of-plane coordinates of the elem nodes)
      jac=matmul(midcoords(2:3,:),dn(1:NNODE/2,:)) ! only the first half of nodes' shape func are used
      det=determinant(jac)
      
  ! Nmatrix: ujump (at each int pnt) = Nmatrix*u
      do i = 1,NDIM
          do j = 1,NNODE/2
          Nmatrix(i,i+(j-1)*NDIM)= -fn(j)
          Nmatrix(i,i+(j-1+NNODE/2)*NDIM)=fn(j)
          end do
      end do
      
    ! calculate ujump: disp. jump of the two crack surface, in global coords
      ujump=matmul(Nmatrix,u)
      
      
      ! calculate separation delta in local coords: delta=Qmatrix*ujump
      delta=matmul(Qmatrix,ujump)


      
      
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
      
      ! get D matrix dee accord. to material properties, and update intg point variables
      select case (mattype)
          case ('interface')
          
              if(nofail) then
              ! do not pass in sdv, then no cohesive law can be done, only linear elasticity stiffness and stress calculated
                  call ddsdde(lib_interface(typekey), Dee, jump=delta, stress=Tau)
              else
              ! calculate D matrix, update stress, and iterating sdv
                  call ddsdde(lib_interface(typekey), Dee, jump=delta, stress=Tau, sdv=ig_sdv(2), dfail=dfail) 
              end if
              
          case default
              write(msg_file,*) 'material type not supported for cohesive element!'
              call exit_function
      end select

     


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
              K_matrix(j,k) = K_matrix(j,k) + NtQtDQN(j,k) * det * igwt(kig)
          end do
          F_vector(j) = F_vector(j) + NtQtTau(j) * igwt(kig) * det
      end do
         
      
      
      !------------------------------------------------!
      !   update ig point kig
      !------------------------------------------------!

      !- calculate integration point physical coordinates (initial)
      tmpx    =   matmul(coords,fn)
      
      !- calculate integration point displacement
      do j=1, NDIM
          do i=1, NNODE
              tmpu(j) = tmpu(j) + fn(i) * u((i-1)*NDIM+j)
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
      

      
      deallocate(ig_sdv)
      
  end do !-looped over all int points. ig=NIGPOINT


  if(allocated(xj)) deallocate(xj) 
  if(allocated(uj)) deallocate(uj) 
  if(allocated(ig_sdv)) deallocate(ig_sdv)          

end subroutine integrate_coh3d8_element












! the rest are private subroutines









pure subroutine init_ig (xi, wt, gauss)

  real(kind=dp), intent(inout)  :: xi(NDIM-1,NIGPOINT), wt(NIGPOINT)
  logical, optional, intent(in) :: gauss
  

    

    if (NIGPOINT .eq. 4) then   

        if(present(gauss).and.gauss) then
            xi(1,1)= -root3
            xi(2,1)= -root3
            
            xi(1,2)= root3
            xi(2,2)= -root3
            
            xi(1,3)= -root3
            xi(2,3)= root3
            
            xi(1,4)= root3
            xi(2,4)= root3
        else
            xi(1,1)=-one
            xi(2,1)=-one
            
            xi(1,2)=one
            xi(2,2)=-one
            
            xi(1,3)=-one
            xi(2,3)=one
            
            xi(1,4)=one
            xi(2,4)=one
        end if
        wt=one
        
    else
        write(msg_file,*) 'no. of integration points incorrect for coh3d8_ig!'
        call exit_function
    end if

end subroutine init_ig



pure subroutine init_shape (igxi, f, df)
  
    real(kind=dp),intent(inout) :: f(NNODE),df(NNODE,NDIM-1)
    real(kind=dp),intent(in) :: igxi(NDIM-1)
    
    real(kind=dp) :: xi, eta ! local variables
    xi=zero
    eta=zero

    xi=igxi(1)
    eta=igxi(2)

    f(1)=quarter*(one-xi)*(one-eta)
    f(2)=quarter*(one+xi)*(one-eta)
    f(3)=quarter*(one+xi)*(one+eta)
    f(4)=quarter*(one-xi)*(one+eta)
    df(1,1) = -quarter*(one-eta)
    df(2,1) =  quarter*(one-eta)
    df(3,1) =  quarter*(one+eta)
    df(4,1) = -quarter*(one+eta)
    df(1,2) = -quarter*(one-xi)
    df(2,2) = -quarter*(one+xi)
    df(3,2) =  quarter*(one+xi)
    df(4,2) =  quarter*(one-xi)
    
    f(1+4)=quarter*(one-xi)*(one-eta)
    f(2+4)=quarter*(one+xi)*(one-eta)
    f(3+4)=quarter*(one+xi)*(one+eta)
    f(4+4)=quarter*(one-xi)*(one+eta)
    df(1+4,1) = -quarter*(one-eta)
    df(2+4,1) =  quarter*(one-eta)
    df(3+4,1) =  quarter*(one+eta)
    df(4+4,1) = -quarter*(one+eta)
    df(1+4,2) = -quarter*(one-xi)
    df(2+4,2) = -quarter*(one+xi)
    df(3+4,2) =  quarter*(one+xi)
    df(4+4,2) =  quarter*(one-xi)

end subroutine init_shape




end module coh3d8_element_module
