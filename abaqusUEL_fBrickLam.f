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
use global_clock_module, only: GLOBAL_CLOCK, set
use input_module,        only: set_fnm_nodes, set_fnm_edges, set_fnm_elems, set_fnm_materials
use output_module,       only: outdir, output

  implicit none

  real(kind=DP), intent(in)   :: time(2)
  real(kind=DP), intent(in)   :: dtime
  integer,       intent(in)   :: lop,lrestart,kstep,kinc

  ! local variables
  character(len=DIRLENGTH)    :: workdir
  integer                     :: lenworkdir
  integer                     :: freq
  !~integer                     :: istat
  !~character(len=MSGLENGTH)    :: emsg

  workdir    = ''
  lenworkdir = 0
  freq       = 0
  !~istat      = STAT_SUCCESS
  !~emsg       = ''

  ! output frequency x: once every x increments
  freq = 1
  
  select case (lop)
  
  ! start of the analysis
  case (0)
      ! get output directory (global variable defined in output module)
      outdir = ''
      call getoutdir(workdir, lenworkdir)
      if (  DIRLENGTH < lenworkdir+len('/outputs/')  ) then
        write(MSG_FILE,*)'increase DIRLENGTH parameter to:',lenworkdir+len('/outputs/')
        call EXIT_FUNCTION
      end if
      outdir = trim(workdir)//'/outputs/'


      
      ! initialize global clock and list of nodes, edges, elems and materials
      call set(GLOBAL_CLOCK, curr_step=0, curr_inc=0)
      
      call set_fnm_nodes
      
      call set_fnm_edges
      
      call set_fnm_elems
      
      call set_fnm_materials
      
      ! open a file 
      open(110, file=trim(outdir)//'record.dat', status="replace", action="write")
      write(110,'(1X, a)')'reach mark 1'
      close(110)
      
  
  ! start of increment
  case (1)
      call set(GLOBAL_CLOCK, curr_step=kstep, curr_inc=kinc)
    
  ! end of increment
  case (2)
      ! print element outputs after certain increment         
      if ( kinc==1 .or. mod(kinc,freq)==0 ) then
        call output(kstep,kinc,outdir)
        !~if(istat == STAT_FAILURE)
        !~  write(MSG_FILE,*) emsg
        !~  call EXIT_FUNCTION
        !~end if
      end if

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
use parameter_module,     only: NDIM, DP, ZERO, MSG_FILE, MSGLENGTH, STAT_SUCCESS, &
                                & STAT_FAILURE, EXIT_FUNCTION
use fnode_module,         only: fnode, update
use node_list_module,     only: node_list
use edge_list_module,     only: edge_list
use elem_list_module,     only: elem_list, elem_node_connec, elem_edge_connec
use material_list_module, only: UDSinglePly_material, matrixCrack_material, &
                                & interface_material
use fBrickLam_elem_module,only: fBrickLam_elem, integrate
use output_module,        only: outdir

  ! use Abaqus default implict type declaration for passed-in variables only
  include 'aba_param.inc'

  dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*)
  dimension svars(*),energy(8),coords(mcrd,nnode),u(ndofel)
  dimension du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*)
  dimension jdltyp(mdload,*),adlmag(mdload,*),ddlmag(mdload,*)
  dimension predef(2,npredf,nnode),lflags(*),jprops(*)

  ! all local variables
  integer                 :: istat
  character(len=MSGLENGTH):: emsg
  character(len=MSGLENGTH):: msgloc
  real(DP),allocatable    :: Kmat(:,:), Fvec(:)
  real(DP)                :: uj(NDIM,nnode)
  !~type(fBrickLam_elem)    :: elem
  type(fnode)             :: nodes(nnode)
  integer                 :: node_cnc(nnode)
  integer                 :: nedge
  integer, allocatable    :: edge_status(:)
  integer, allocatable    :: edge_cnc(:)
  integer                 :: j

  amatrx    = ZERO
  rhs(:,1)  = ZERO
  istat     = STAT_SUCCESS
  emsg      = ''
  msgloc    = ', abaqusUEL_fBrickLam.f'
  uj        = ZERO
  node_cnc  = 0
  nedge     = 0
  j         = 0
  
!~  ! check input validity during first run
!~  if (kstep == 1 .and. kinc == 1) then
!~    ! check lflag(1)
!~    ! for static analysis, lflag(1)=1, 2
!~    if (.not. (lflags(1)==1 .or. lflags(1)==2)) then
!~      istat = STAT_FAILURE
!~      emsg  = 'this uel is only for static analysis'//trim(msgloc)
!~      goto 10
!~    end if
!~    ! check lflag(3)
!~    if (lflags(3) /= 1) then
!~      istat = STAT_FAILURE
!~      emsg  = 'this uel is only for normal increment'//trim(msgloc)
!~      goto 10
!~    end if
!~    ! check NRHS
!~    ! for static analysis, nrhs = 1
!~    ! for modified riks static analysis, nrhs = 2
!~    if (nrhs /= 1) then
!~      istat = STAT_FAILURE
!~      emsg  = 'this uel is NOT for riks method'//trim(msgloc)
!~      goto 10
!~    end if
!~    ! check mcrd
!~    if (mcrd /= NDIM) then
!~      istat = STAT_FAILURE
!~      emsg  = 'uel mcrd does not match param NDIM'//trim(msgloc)
!~      goto 10
!~    end if
!~    ! check nnode
!~    if (nnode /= size(elem_node_connec(:,jelem))) then
!~      istat = STAT_FAILURE
!~      emsg  = 'uel nnode does not match elem_node_connec size'//trim(msgloc)
!~      goto 10
!~    end if
!~    ! check ndofel
!~    if (ndofel /= NDIM*nnode) then
!~      istat = STAT_FAILURE
!~      emsg  = 'uel ndofel does not match NDIM*NNODE'//trim(msgloc)
!~      goto 10
!~    end if
!~
!~10  if (istat == STAT_FAILURE) then
!~      write(MSG_FILE,*) emsg
!~      call cleanup_all
!~      call EXIT_FUNCTION
!~    end if
!~  end if

  
  ! extract this element from global elem list, using the element key jelem
  associate ( elem  => elem_list(jelem))
  
    ! extract the node and edge connec of this elem
    node_cnc(:) = elem_node_connec(:,jelem)
    nedge       = size(elem_edge_connec(:,jelem))
    allocate(edge_status(nedge),edge_cnc(nedge))
    edge_cnc(:) = elem_edge_connec(:,jelem)

    ! extract passed-in nodal solutions obtained by Abaqus Solver
    do j=1, nnode
      uj(1:NDIM,j) = u( (j-1)*NDIM+1 : j*NDIM )
    end do
    
    ! update the nodal solutions to global_node_list
    do j=1, nnode
      call update(node_list(node_cnc(j)),u=uj(:,j))
    end do

    ! extract nodes and edge status from global node and edge lists
    nodes       = node_list(node_cnc)
    edge_status = edge_list(edge_cnc)

    ! open a file 
    open(110, file=trim(outdir)//'record.dat', status="replace", action="write")
    write(110,'(1X, a)')'reach mark 2'
    close(110)

    ! integrate this element
    call integrate (elem, nodes, edge_status, UDSinglePly_material, &
    &  matrixCrack_material, interface_material, Kmat, Fvec, istat, emsg)
    if (istat == STAT_FAILURE) then
      emsg = emsg//trim(msgloc)
      write(MSG_FILE,*) emsg
      call cleanup (Kmat, Fvec, edge_status, edge_cnc)
      call cleanup_all
      call EXIT_FUNCTION
    end if
  
    ! open a file 
    open(110, file=trim(outdir)//'record.dat', status="replace", action="write")
    write(110,'(1X, a)')'reach mark 3'
    close(110)

    ! update to global lists
    !~elem_list(jelem)    = elem
    node_list(node_cnc) = nodes
    edge_list(edge_cnc) = edge_status
  
  end associate
  
  ! open a file 
  open(110, file=trim(outdir)//'record.dat', status="replace", action="write")
  write(110,'(1X, a)')'reach mark 4'
  close(110)

  ! in the end, pass Kmat and Fvec to Abaqus UEL amatrx and rhs
  amatrx   =  Kmat
  rhs(:,1) = -Fvec(:)

  ! clean up memory used in local dynamic arrays
  call cleanup (Kmat, Fvec, edge_status, edge_cnc)
  return
  
  contains
  
  pure subroutine cleanup (Kmat, Fvec, edge_status, edge_cnc)
    real(DP),allocatable, intent(inout)    :: Kmat(:,:), Fvec(:)
    integer, allocatable, intent(inout)    :: edge_status(:)
    integer, allocatable, intent(inout)    :: edge_cnc(:)
    if(allocated(Kmat))         deallocate(Kmat)
    if(allocated(Fvec))         deallocate(Fvec)
    if(allocated(edge_status))  deallocate(edge_status)
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