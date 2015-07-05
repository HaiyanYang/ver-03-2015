!***************************************!
!   Abaqus UEL interface of the FNM     !
!                                       !
!                                       !
!***************************************!


!------ include FNM modules --------------------------
include 'globals/parameter_module.f90'
include 'globals/global_clock_module.f90'
include 'globals/global_toolkit_module.f90'
include 'object_materials/lamina_material_module.f90'
include 'object_materials/cohesive_material_module.f90'
include 'object_node/fnode_module.f90'
include 'object_elements/base_elements/brickPly_elem_module.f90'
include 'object_elements/base_elements/wedgePly_elem_module.f90'
include 'object_elements/base_elements/coh8Crack_elem_module.f90'
include 'object_elements/base_elements/coh6Delam_elem_module.f90'
include 'object_elements/base_elements/coh8Delam_elem_module.f90'
include 'object_elements/base_elements/abstPly_elem_module.f90'
include 'object_elements/base_elements/abstDelam_elem_module.f90'
include 'object_elements/fBrickPly_elem_module.f90'
include 'object_elements/fCoh8Delam_subelem_module.f90'
include 'object_elements/fCoh8Delam_elem_module.f90'
include 'object_elements/fBrickLam_elem_module.f90'
include 'datalists/material_list_module.f90'
include 'datalists/node_list_module.f90'
include 'datalists/edge_list_module.f90'
include 'datalists/elem_list_module.f90'
include 'inputs/input_module.f90'
include 'outputs/output_module.f90'
!------------------------------------------------------

!write(msg_file,*) 'reach here'

!---------------------------------------------------------!
!   Abaqus user subroutine for I/O to external files
!---------------------------------------------------------!
subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
use parameter_module,    only: DP, DIRLENGTH, MSG_FILE, EXIT_FUNCTION
use global_clock_module, only: GLOBAL_CLOCK, set_program_clock
use input_module
use output_module

  implicit none

  real(kind=DP), intent(in)   :: time(2)
  real(kind=DP), intent(in)   :: dtime
  integer,       intent(in)   :: lop,lrestart,kstep,kinc

  ! local variables
  character(len=DIRLENGTH)    :: workdir
  integer                     :: lenworkdir
  integer                     :: freq

  workdir    = ''
  lenworkdir = 0
  freq       = 0

  ! output frequency x: once every x increments
  freq = 1
  
  select case (lop)
  
  ! start of the analysis
  case (0)
      ! initialize global clock and list of nodes, edges, elems and materials
      call set_program_clock(GLOBAL_CLOCK, curr_step=kstep, curr_inc=kinc)
      call set_fnm_nodes
      call set_fnm_edges
      call set_fnm_elems
      call set_fnm_materials
      ! get output directory (global variable defined in output module)
      outdir = ''
      call getoutdir(workdir, lenworkdir)
      if (  DIRLENGTH < lenworkdir+len('/outputs/')  ) then
        write(MSG_FILE,*)'increase DIRLENGTH parameter to:',lenworkdir+len('/outputs/')
        call EXIT_FUNCTION
      end if
      outdir = trim(workdir)//'/outputs/'
  
  ! start of increment
  case (1)
      call set_program_clock(GLOBAL_CLOCK, curr_step=kstep, curr_inc=kinc)
    
  ! end of increment
  case (2)
      ! print element outputs after certain increment         
      if ( kinc==1 .or. mod(kinc,freq)==0 ) call output(kstep,kinc,outdir)

  ! end of analysis
  case (3)
      call cleanup_all
  
  end select


end subroutine uexternaldb





!---------------------------------------------------------!
!   Abaqus user element subroutine
!---------------------------------------------------------!
subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars, &
 &       props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime, &
 &       kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
 &       lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)
! load FNM modules
use parameter_module,       only: DP, ZERO, MSG_FILE, MSGLENGTH, STAT_SUCCESS, &
                                & STAT_FAILURE, EXIT_FUNCTION
use fnode_module,           only: fnode
use global_node_list_module,only: global_node_list
use global_edge_list_module,only: global_edge_list
use global_elem_list_module,only: global_elem_list, elem_node_connec, elem_edge_connec
use global_material_module, only: UDSinglePly_material, matrixCrack_material, &
                                & interface_material

  ! use Abaqus default implict type declaration for passed-in variables only
  include 'aba_param.inc'

  dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*)
  dimension svars(*),energy(8),coords(mcrd,nnode),u(ndofel)
  dimension du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*)
  dimension jdltyp(mdload,*),adlmag(mdload,*),ddlmag(mdload,*)
  dimension predef(2,npredf,nnode),lflags(*),jprops(*)


  ! explicitly define all local variables
  integer                 :: istat
  character(len=MSGLENGTH):: emsg
  character(len=MSGLENGTH):: msgloc
  real(DP),allocatable    :: Kmat(:,:), Fvec(:)
  real(DP)                :: uj(mcrd,nnode)
  type(fBrickLam_element) :: elem
  type(fnode)             :: nodes(nnode)
  integer, allocatable    :: edge_status(:)
  integer, allocatable    :: node_cnc(:), edge_cnc(:)
  integer                 :: j

  amatrx    = ZERO
  rhs(:,1)  = ZERO
  istat     = STAT_SUCCESS
  emsg      = ''
  msgloc    = ', abaqusUEL_fBrickLam.f'
  Kmat      = ZERO
  Fvec      = ZERO
  uj        = ZERO
  j         = 0
  

  ! extract this element from global elem list, using the element key jelem
  elem        = global_elem_list(jelem)
  node_cnc(:) = elem_node_connec(:,jelem)
  edge_cnc(:) = elem_edge_connec(:,jelem)

  ! extract passed-in nodal solutions obtained by Abaqus nonlinear Solver
  do j=1, nnode
    uj(1:mcrd,j) = u( (j-1)*mcrd+1 : j*mcrd )
  end do
  
  ! update the nodal solutions to global_node_list
  do j=1, nnode
    call update(global_node_list(node_cnc(j)),u=uj(:,j))
  end do

  ! extract nodes and edge status from global node and edge lists
  nodes       = global_node_list(node_cnc)
  edge_status = global_edge_list(edge_cnc)

  ! integrate this element
  call integrate (elem, nodes, edge_status, UDSinglePly_material, &
  &  matrixCrack_material, interface_material, Kmat, Fvec, istat, emsg)
  if (istat == STAT_FAILURE) then
    emsg = emsg//trim(msgloc)
    write(MSG_FILE,*) emsg
    call cleanup (Kmat, Fvec, edge_status, node_cnc, edge_cnc)
    call cleanup_all
    call EXIT_FUNCTION
  end if

  ! update to global lists
  global_elem_list(jelem)    = elem
  global_node_list(node_cnc) = nodes
  global_edge_list(edge_cnc) = edge_status

  ! in the end, pass Kmat and Fvec to Abaqus UEL amatrx and rhs
  amatrx   =  Kmat
  rhs(:,1) = -Fvec(:)

  ! clean up memory used in local dynamic arrays
  call cleanup (Kmat, Fvec, edge_status, node_cnc, edge_cnc)
  return
  
  contains
  
  pure subroutine cleanup (Kmat, Fvec, edge_status, node_cnc, edge_cnc)
    real(DP),allocatable, intent(inout)    :: Kmat(:,:), Fvec(:)
    integer, allocatable, intent(inout)    :: edge_status(:)
    integer, allocatable, intent(inout)    :: node_cnc(:), edge_cnc(:)
    if(allocated(Kmat))         deallocate(Kmat)
    if(allocated(Fvec))         deallocate(Fvec)
    if(allocated(edge_status))  deallocate(edge_status)
    if(allocated(node_cnc))     deallocate(node_cnc)
    if(allocated(edge_cnc))     deallocate(edge_cnc)
  end subroutine cleanup

end subroutine uel





!---------------------------------------------------------!
!   subroutine to clean up all the datalist
!---------------------------------------------------------!
subroutine cleanup_all()
use material_list_module, only: empty_material_list
use node_list_module,     only: empty_node_list
use edge_list_module,     only: empty_edge_list
use elem_list_module,     only: empty_elem_list

  call empty_material_list
  call empty_node_list
  call empty_edge_list
  call empty_elem_list

end subroutine cleanup_all