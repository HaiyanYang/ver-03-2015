!***************************************!
!   Abaqus UEL interface of the FNM     !
!                                       !
!                                       !
!***************************************!


!------ include useful modules--------------------------
include 'globals/parameter_module.f90'
include 'globals/global_clock_module.f90'
include 'globals/global_toolkit_module.f90'
include 'libraries/lib_node_module.f90'
include 'libraries/lib_edge_module.f90'
include 'libraries/lib_mat_module.f90'
include 'libraries/lib_bcd_module.f90'
include 'libraries/lib_elem_module.f90'
include 'libraries/initialize_lib_module.f90'
include 'outputs/output_module.f90'
!------------------------------------------------------

!write(msg_file,*) 'reach here'

!---------------------------------------------------------!
!   Abaqus user subroutine for I/O to external files
!---------------------------------------------------------!
subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
use parameter_module
use glb_clock_module
use initialize_lib_module
use output_module

implicit none

real(kind=dp), intent(in)   :: time(2)
real(kind=dp), intent(in)   :: dtime
integer,       intent(in)   :: lop,lrestart,kstep,kinc

! local variables
character(len=dirlength)    :: workdir
integer                     :: lenworkdir
integer                     :: freq
real(dp)                    :: prvdtime,crrpnewdt

    workdir=''
    lenworkdir=0
    freq=0
    prvdtime=zero
    crrpnewdt=zero

    ! output frequency x: once every x increments
    freq=1

    if (lop .eq. 0) then
!       start of the analysis
        

        ! initialize global clock and libraries
        call initialize_glb_clock
        call initialize_lib_node
        call initialize_lib_edge
        call initialize_lib_elem
        call initialize_lib_mat 
        call initialize_lib_bcd
         
        ! get output directory (global variable defined in output module)
        
        outdir=''
        
        call getoutdir(workdir, lenworkdir)
        
        if(dirlength<lenworkdir+len('/outputs/')) then
            write(msg_file,*)'increase dirlength parameter to:',lenworkdir+len('/outputs/')
            call exit_function
        end if
        
        outdir=trim(workdir)//'/outputs/'

    else if (lop .eq. 1) then
!       start of the current increment
        !~call extract_glb_clock(dtime=prvdtime)            
        !~
        !~if(abs(dtime-prvdtime)>tiny(one)) then
        !~! a change in dtime, update pnewdt
        !~    if(prvdtime>tiny(one)) then
        !~        crrpnewdt=dtime/prvdtime
        !~        call update_glb_clock(pnewdt=crrpnewdt)
        !~    end if
        !~end if
        
        call update_glb_clock(kstep,kinc,time(2),dtime)
        
        
    
    else if (lop .eq. 2) then
!	    end of the increment 

        ! print element outputs after certain increment         
       if(kinc==1 .or. mod(kinc,freq)==0) call output(kstep,kinc,outdir)
        
    
    else if (lop .eq. 3) then
!	    end of the analysis
        call empty_lib_mat
        call empty_lib_node
        call empty_lib_edge
        call empty_lib_elem
        call empty_lib_bcd
        
  
    end if


return
end subroutine uexternaldb





!---------------------------------------------------------!
!   Abaqus user element subroutine
!---------------------------------------------------------!
subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars, &
 &       props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime, &
 &       kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
 &       lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)
! load FNM modules
use parameter_module,       only:
use xnode_module,           only: xnode
use global_nodelist_module, only: global_node_list
use global_edgelist_module, only: global_edge_list
use global_elemlist_module, only: global_fBrickLam_list
use global_material_module, only: plylam_mat, plycoh_mat, interf_mat

! use Abaqus default implict type declaration for passed-in variables only
include 'aba_param.inc'

dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*)
dimension svars(*),energy(8),coords(mcrd,nnode),u(ndofel)
dimension du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*)
dimension jdltyp(mdload,*),adlmag(mdload,*),ddlmag(mdload,*)
dimension predef(2,npredf,nnode),lflags(*),jprops(*)



! explicitly define all local variables



! element K matrix and F vector
real(dp),allocatable    :: Kmat(:,:), Fvec(:)

! variables of element nodes, to be extracted from uel passed-in variables
! u, to be updated to global node array 'global_node_list' according to
! the element nodal connectivity
real(dp)                :: uj(mcrd,nnode)

! create this element
type(fBrickLam_element) :: elem

! elem nodes
type(xnode)             :: nodes(nnode)

! elem edge status
integer                 :: edge_status(nedge)

! nodal cnc of this elem
integer, allocatable    :: node_cnc(:), edge_cnc(:)

! counters
integer :: i,j,l

uj=zero
i=0; j=0; l=0

amatrx   = zero
rhs(:,1) = zero

!---------------------------------------------------------------------------------!
!           Abaqus UEL interface to FNM codes for variable pass-in
!---------------------------------------------------------------------------------!

! copy this element definition from glb list, using the element key jelem
elem = global_fBrickLam_list(jelem)

! extract nodal connec of this elem
call extract(elem, node_connec=node_cnc, edge_connec=edge_cnc)

! extract nodal values from passed in variables
do j=1, nnode
  uj(1:mcrd,j) = u( (j-1)*mcrd+1 : j*mcrd )
end do

! update nodal values to glb node library lib_node according to cnc
do j=1, nnode
  call update(global_node_list(node_cnc(j)),u=uj(:,j))
end do

! extract info needed for integration
nodes       = global_node_list(node_cnc)
edge_status = global_edge_list(edge_cnc)

! integrate this element
call integrate (elem, nodes, edge_status, plylam_mat, plycoh_mat, interf_mat, &
& K_matrix, F_vector, istat, emsg)
if (istat == STAT_FAILURE) then
  emsg = 'error in integrate, abaqusUEL_fBrickLam,f'
  return
end if

! update intent inout global lists
global_fBrickLam_list(jelem) = elem
global_node_list(node_cnc)   = nodes
global_edge_list(edge_cnc)   = edge_status

! in the end, pass Kmat and Fvec to Abaqus UEL amatrx and rhs
amatrx   =  Kmat
rhs(:,1) = -Fvec(:)

if(allocated(Kmat)) deallocate(Kmat)
if(allocated(Fvec)) deallocate(Fvec)
if(allocated(node_cnc)) deallocate(node_cnc)
if(allocated(edge_cnc)) deallocate(edge_cnc)

end subroutine uel
