module output_module
use parameter_module, only: DIRLENGTH,MSGLENGTH,ELTYPELENGTH,DP,ZERO
use node_list_module, only: node_list
use elem_list_module, only: elem_list
use fnode_module,     only: extract
use brickPly_elem_module, only: brickPly_elem, extract
use abstPly_elem_module,  only: abstPly_elem,  extract
use coh8Crack_elem_module, only: coh8Crack_elem, extract
use coh8Delam_elem_module, only: coh8Delam_elem, extract
use abstDelam_elem_module, only: abstDelam_elem, extract
use fBrickPly_elem_module, only: fBrickPly_elem, extract
use fCoh8Delam_subelem_module, only: fCoh8Delam_subelem, extract
use fCoh8Delam_elem_module, only: fCoh8Delam_elem, extract

implicit none
private

! define global variable for output directory
character(len=DIRLENGTH), save :: outdir


public :: outdir, output


contains


subroutine output(kstep,kinc,outdir,istat,emsg)
   
  ! passed-in variables
  integer,                  intent(in)  :: kstep    ! current step number
  integer,                  intent(in)  :: kinc     ! increment number of current step
  character(len=DIRLENGTH), intent(in)  :: outdir   ! output directory name
  integer,                  intent(out) :: istat
  character(len=MSGLENGTH), intent(out) :: emsg
  
  ! local variables
  
  character(len=MSGLENGTH)    :: msgloc
  
  ! output file variables
  integer                     :: outunit  ! output file unit
  character(len=DIRLENGTH)    :: outfile  ! output file name
  character(len=DIRLENGTH)    :: outnum   ! output increment number (embedded in the outfile name)
  character(len=10)           :: FMAT, FMATKINC, FMATNNODE

  ! no. of nodes & elem in the mesh
  integer                     :: nnode, nelem
  
  ! no. of elems of each elem type
  integer                     :: nwedge, nbrick, ncoH6, ncoh8
  integer                     :: nfLam, nplyblk, ninterf, nsubBulk, nsubinterf
  
  ! sub element type variable
  character(len=ELTYPELENGTH) :: eltype
  
  ! nsize: size of vtk element output; elnode: no. nodes in each elem
  integer                     :: nsize

  ! connectivity of an elem extracted from lib_elem
  integer, allocatable        :: connec(:)
  
  ! nodal coordinates
  real(kind=DP), allocatable  :: x(:)     ! coordinates of nodes extracted from lib_node
  
  ! nodal displacements
  real(kind=DP), allocatable  :: disp(:)  ! displacements of nodes extracted from lib_node

  ! Laminate sub elems
  type(fBrickPly_elem),     allocatable :: plyblks(:)
  type(fCoh8Delam_elem),    allocatable :: interfs(:)
  ! plyblks sub elems
  type(brickPly_elem),      allocatable :: intactplyblk
  type(abstPly_elem),       allocatable :: subBulks(:)
  type(coh8Crack_elem),     allocatable :: cohCrack
  ! interfs sub elems
  type(coh8Delam_elem),     allocatable :: intactinterf
  type(fCoh8Delam_subelem), allocatable :: top_interf
  type(fCoh8Delam_subelem), allocatable :: bot_interf
  type(abstDelam_elem),     allocatable :: subinterfs(:)
  
  ! stress and strain vectors
  real(DP) :: str(6)
  real(DP) :: str0(6)
  real(DP) :: scalar
  
  ! counters
  integer  :: i, j, l, m, n, nfl



  ! -----------------------------------------------------------------!
  !                       initialize variables
  ! -----------------------------------------------------------------!
  istat  = STAT_SUCCESS
  emsg   = ''
  msgloc = ', output module'
  
  outunit   = 0
  outfile   = ''
  outnum    = ''
  FMAT      = ''
  FMATKINC  = ''
  FMATNNODE = ''
  eltype    = ''
  nnode     = 0
  nfLam     = 0
  nplyblk   = 0
  ninterf   = 0
  nsubBulk  = 0 
  nsubinterf= 0
  nwedge    = 0
  nbrick    = 0 
  ncoh6     = 0 
  ncoh8     = 0
  nelem     = 0
  nsize     = 0
  str       = ZERO
  str0      = ZERO
  scalar    = ZERO

  nfl=0; i=0; j=0; l=0; m=0; n=0
  

  ! format parameters
  FMATKINC  = '(i5.5)' ! for increment no.
  FMATNNODE = '(i10)'  ! for connec node no.

  ! write the increment number as a character and store in outnum
  write(outnum,FMATKINC) kinc
  
  ! create the output file name
  !~outfile=trim(outdir)//'/outputs/'//trim(outnum)//'.vtk'
  outfile = trim(outdir)//'fnm-'//trim(outnum)//'.vtk'
  outfile = trim(outfile)
  
  ! open the outfile
  open(newunit(outunit), file=outfile, status="replace", action="write")
  
  ! write header
  write(outunit,'1X, (a)')'# vtk DataFile Version 3.1'
  write(outunit,'1X, (a)')'for Floating Node Method output'
  
  ! write vtk format
  write(outunit,'1X, (a)')'ASCII'      
  
  ! write vtk data type
  write(outunit,'1X, (a)')'DATASET UNSTRUCTURED_GRID'      
  


  ! -----------------------------------------------------------------!
  !                     write nodes
  ! -----------------------------------------------------------------!
  ! obtain nnode value
  nnode = size(node_list)
  ! write node summary line
  write(outunit,'(1X, a, '//trim(FMATNNODE)//', a)')'POINTS ',nnode,' FLOAT'
  ! write all nodes
  do i = 1,nnode
      ! extract nodal coords from lib_node
      call extract(node_list(i),x=x)
      ! write x3d array into output file
      write(outunit,*) x(1),x(2),x(3)
  end do            
  write(outunit,'(1X)')

  
  ! -----------------------------------------------------------------!
  !                     write element connecs
  ! -----------------------------------------------------------------!  
  ! wflag = 'nelem' would allow the wfLam subroutine to update
  ! nwedge, nbrick, ncoh6 and ncoh8
  call wfLam('nelem')

  ! total no. of elems
  nelem = nwedge + nbrick + ncoh6 + ncoh8
  
  ! calculate total no. of nodes to print; each row has 1+elnode no. of indices to print
  nsize = (nwedge+ncoh6) * (1+6) + (nbrick+ncoh8) * (1+8)
  
  ! write a summary of output
  write(outunit,'(1X, a, '//trim(FMATNNODE)//', '//trim(FMATNNODE)//')')'CELLS ', nelem, nsize
  ! wflag = 'connec' would allow the wfLam subroutine to write out all
  ! the nodal connec of all sub elems in the mesh
  call wfLam('connec')


  ! -----------------------------------------------------------------!
  !                     write element types (order matters)
  ! -----------------------------------------------------------------!  
  write(outunit,'(1X, a, i10)')'CELL_TYPES ', nelem
  ! wflag = 'eltype' would allow the wfLam subroutine to write out all
  ! the elem type code (used by VTK) of all sub elems in the mesh
  call wfLam('eltype')



  ! -----------------------------------------------------------------!
  ! BEGIN WRITE NODAL DATA
  ! -----------------------------------------------------------------!
  write(outunit,'(1X, a, i10)')'POINT_DATA ', nnode

  
  ! -----------------------------------------------------------------!
  !                     write displacements
  ! -----------------------------------------------------------------!
  write(outunit,'(1X, a)')'VECTORS displacement float'
  
  do i=1, nnode
      call extract(node_list(i),u=disp)
      write(outunit,*) disp(1),disp(2),disp(3)
      deallocate(disp)
  end do  
          
  write(outunit,'(1X)')  


  ! -----------------------------------------------------------------!
  ! BEGIN WRITE CELL DATA
  ! -----------------------------------------------------------------!
  write(outunit,'(1X, a, i10)')'CELL_DATA ', nelem
  
  
  ! -----------------------------------------------------------------!
  !                     write stress (order matters)
  ! -----------------------------------------------------------------!     
  write(outunit,'(1X, a)')'TENSORS stress float'
  ! wflag = 'stress' would allow the wfLam subroutine to write out all
  ! stress tensors of all ply sub elems in the mesh
  call wfLam('stress')

  
  ! -----------------------------------------------------------------!
  !                     write strain (order matters)
  ! -----------------------------------------------------------------!  
  write(outunit,'(1X, a)')'TENSORS strain float'
  ! wflag = 'strain' would allow the wfLam subroutine to write out all
  ! strain tensors of all ply sub elems in the mesh
  call wfLam('strain')
 

  ! -----------------------------------------------------------------!
  !                     write failure status
  ! -----------------------------------------------------------------! 
  write(outunit,'(1X, a)')'SCALARS fstat float'
  write(outunit,'(1X, a)')'LOOKUP_TABLE default'
  ! wflag = 'fstat' would allow the wfLam subroutine to write out all
  ! failure status var. of all sub elems in the mesh
  call wfLam('fstat')
  

  ! -----------------------------------------------------------------!
  !                     write fibre damage variable
  ! -----------------------------------------------------------------!
  write(outunit,'(1X, a)')'SCALARS df float'
  write(outunit,'(1X, a)')'LOOKUP_TABLE default'
  ! wflag = 'df' would allow the wfLam subroutine to write out all
  ! fibre degradation factors of all ply sub elems in the mesh
  call wfLam('df')


  ! -----------------------------------------------------------------!
  !                     write matrix damage variable
  ! -----------------------------------------------------------------! 
  write(outunit,'(1X, a)')'SCALARS dm float'
  write(outunit,'(1X, a)')'LOOKUP_TABLE default'
  ! wflag = 'dm' would allow the wfLam subroutine to write out all 
  ! cohesive degradation factors of all cohesive sub elems in the mesh
  call wfLam('dm')
  
  
  ! -----------------------------------------------------------------!
  !             deallocate local dynamic arrays and return
  ! -----------------------------------------------------------------!   
  if(allocated(connec))         deallocate(connec)
  if(allocated(x))              deallocate(x)
  if(allocated(disp))           deallocate(disp)
  if(allocated(plyblks))        deallocate(plyblks)
  if(allocated(interfs))        deallocate(interfs)  
  if(allocated(intactplyblk))   deallocate(intactplyblk)
  if(allocated(subBulks))       deallocate(subBulks)
  if(allocated(cohCrack))       deallocate(cohCrack)
  if(allocated(intactinterf))   deallocate(intactinterf)
  if(allocated(top_interf))     deallocate(top_interf)
  if(allocated(bot_interf))     deallocate(bot_interf)
  if(allocated(subinterfs))     deallocate(subinterfs)
  
  close(outunit)
  
  return
  
  
  
  contains
  
  
    subroutine wfLam(wflag)
      character(len=*), intent(in)  :: wflag
      
      nfLam = size(elem_list)
    
      do nfl = 1, nfLam
      
          ! extract the plyblks and interfs elems of the fLam elem
          call extract(elem_list(nfl), plyblks=plyblks, interfs=interfs)
          
          ! read the plyblks elems
          if(allocated(plyblks)) then
              ! read the no. of plyblks
              nplyblk = size(plyblks)
              ! loop over each plyblk subelem
              do m = 1, nplyblk
                  ! extract the subelems of this plyblk
                  call extract(plyblks(m), intact_elem=intactplyblk, &
                  & subBulks=subBulks, cohCrack=cohCrack)
                  ! write intact elem
                  if (allocated(intactplyblk)) then                      
                      call wbrickply(intactplyblk,wflag)
                  ! write sub bulks and coh crack
                  else
                      ! write sub bulks
                      nsubBulk = size(subBulks)
                      do n = 1, nsubBulk
                        call wabstply(subBulk(n),wflag)
                      end do
                      deallocate(subBulks)
                      ! write coh crack
                      call wcohcrack(cohCrack,wflag)
                  end if
              end do
              deallocate(plyblks) 
          end if
          
          ! read the interfs elems
          if(allocated(interfs)) then
              ! read the no. of interfs
              ninterf = size(interfs)
              ! loop over all interfs
              do m = 1, ninterf
                  ! extract the subinterfs
                  call extract(interfs(m), intact_elem=intactinterf, &
                  & top_subelem=top_interf, bot_subelem=bot_interf)
                  
                  ! if intact elem is still present, then it is a coh8 elem type
                  if (allocated(intactinterf)) then
                      call wcoh8Delam(intactinterf,wflag) 
                  ! if not, then it has decomposed into sub elems
                  else
                      ! if top interf is present
                      if (allocated(top_interf)) then
                          call extract(top_interf, subelems=subinterfs)
                          ! find the subinterf type and increase the respective type count
                          nsubinterf = size(subinterfs)
                          do n = 1, nsubinterf
                              call wabstDelam(subinterfs(n), wflag)
                          end do
                          deallocate(subinterfs)
                      end if
                      ! if bot interf is present
                      if (allocated(bot_interf)) then
                          call extract(bot_interf, subelems=subinterfs)
                          ! find the subinterf type and increase the respective type count
                          nsubinterf = size(subinterfs)
                          do n = 1, nsubinterf
                              call wabstDelam(subinterfs(n), wflag)
                          end do
                          deallocate(subinterfs)
                      end if
                  end if
              end do
              deallocate(interfs)
          end if
          
      end do

      write(outunit,'(1X)')
      
    
    end subroutine wfLam
  
    subroutine wbrickply(brickply,wflag)
      type(brickPly_elem), intent(in) :: brickply
      character(len=*),    intent(in) :: wflag
        select case(wflag)
        case('nelem')
          nbrick = nbrick + 1
        case('connec')
          call extract(brickply, connec=connec)
          call wconnec(connec)
        case('eltype')
          call weltype('brick')
        case('stress')
          call extract(brickply, stress=str)
          call wtensor(str)
        case('strain')
          call extract(brickply, strain=str)
          call wtensor(str)
        case('fstat')
          call extract(brickply, fstat=scalar)
          call wscalar(scalar)
        case('df')
          call extract(brickply, df=scalar)
          call wscalar(scalar)
        case('dm')
          call wscalar(ZERO)
        end select  
    end subroutine wbrickply
    
    subroutine wabstply(abstply,wflag)
      type(abstPly_elem), intent(in) :: abstply
      character(len=*),   intent(in) :: wflag
        select case(wflag)
        case('nelem')
          call extract(abstply, eltype=eltype)
          if(eltype=='brick') then
            nbrick = nbrick + 1
          else if(eltype=='wedge') then
            nwedge = nwedge + 1
          end if
        case('connec')
          call extract(abstply, connec=connec)
          call wconnec(connec)
        case('eltype')
          call extract(abstply, eltype=eltype)
          call weltype(eltype)
        case('stress')
          call extract(abstply, stress=str)
          call wtensor(str)
        case('strain')
          call extract(abstply, strain=str)
          call wtensor(str)
        case('fstat')
          call extract(abstply, fstat=scalar)
          call wscalar(scalar)
        case('df')
          call extract(abstply, df=scalar)
          call wscalar(scalar)
        case('dm')
          call wscalar(ZERO)
        end select
    end subroutine wabstply
    
    subroutine wcohcrack(cohCrack,wflag)
      type(cohCrack_elem), intent(in) :: cohCrack
      character(len=*),    intent(in) :: wflag
        select case(wflag)
        case('nelem')
          ncoh8 = ncoh8 + 1
        case('connec')
          call extract(cohCrack, connec=connec)
          call wconnec(connec)
        case('eltype')
          call weltype('coh8')
        case('stress')
          call wtensor(str0)
        case('strain')
          call wtensor(str0)
        case('fstat')
          call extract(cohCrack, fstat=scalar)
          call wscalar(scalar)
        case('df')
          call wscalar(ZERO)
        case('dm')
          call extract(cohCrack, dm=scalar)
          call wscalar(scalar)
        end select
    end subroutine wcohcrack
    
    subroutine wcoh8Delam(coh8Delam,wflag)
      type(coh8Delam_elem), intent(in) :: coh8Delam
      character(len=*),    intent(in) :: wflag
        select case(wflag)
        case('nelem')
          ncoh8 = ncoh8 + 1
        case('connec')
          call extract(coh8Delam, connec=connec)
          call wconnec(connec)
        case('eltype')
          call weltype('coh8')
        case('stress')
          call wtensor(str0)
        case('strain')
          call wtensor(str0)
        case('fstat')
          call extract(coh8Delam, fstat=scalar)
          call wscalar(scalar)
        case('df')
          call wscalar(ZERO)
        case('dm')
          call extract(coh8Delam, dm=scalar)
          call wscalar(scalar)
        end select
    end subroutine wcoh8Delam
   
    subroutine wabstDelam(abstDelam,wflag)
      type(abstDelam_elem), intent(in) :: abstDelam
      character(len=*),     intent(in) :: wflag
        select case(wflag)
        case('nelem')
          call extract(abstDelam, eltype=eltype)
          if (eltype=='coh6Delam') then
            ncoh6 = ncoh6 + 1
          else if (eltype=='coh8Delam') then
            ncoh8 = ncoh8 + 1
          end if
        case('connec')
          call extract(abstDelam, connec=connec)
          call wconnec(connec)
        case('eltype')
          call extract(abstDelam, eltype=eltype)
          call weltype(eltype)
        case('stress')
          call wtensor(str0)
        case('strain')
          call wtensor(str0)
        case('fstat')
          call extract(abstDelam, fstat=scalar)
          call wscalar(scalar)
        case('df')
          call wscalar(ZERO)
        case('dm')
          call extract(abstDelam, dm=scalar)
          call wscalar(scalar)
        end select
    end subroutine wabstDelam

    subroutine wconnec(connec)
      integer, allocatable, intent(inout) :: connec(:)
      ! print connec in vtk; note that in vtk node no. starts from 0
      connec=connec-1
      write(outunit,'1X'//FMATNNODE,advance="no") size(connec)
      do j=1,size(connec)
          write(outunit,FMATNNODE,advance="no") connec(j)
      end do
      write(outunit,'(a)')''
      deallocate(connec)
    end subroutine wconnec
  
    subroutine weltype(eltype)
      character(len=*), intent(in) :: eltype
    
      select case(eltype)
          case('wedge','coh6Delam','coh6')
              write(outunit,'(1X,i2)') 13 ! 13 for wedge/coh6
          case('brick','coh8Delam','coh8')
              write(outunit,'(1X,i2)') 12 ! 12 for brick/coh8
      end select
    
    end subroutine weltype
   
    subroutine wtensor(x)
      real(DP), intent(in) :: x(6)
      real(DP) :: tensor(3,3)
      integer  :: l
    
      tensor = ZERO
      tensor(1,1)=x(1)
      tensor(2,2)=x(2)
      tensor(3,3)=x(3)
      tensor(1,2)=x(4)
      tensor(1,3)=x(5)
      tensor(2,3)=x(6)
      tensor(2,1)=x(4)
      tensor(3,1)=x(5)
      tensor(3,2)=x(6)
    
      do l=1,3
          write(outunit,*) tensor(1,l), tensor(2,l), tensor(3,l)
      end do
      write(outunit,'(a)')''
    
    end subroutine wtensor
       
    subroutine wscalar(x)
      real(DP), intent(in) :: x
      write(outunit,*) x
    end subroutine wscalar
  
  
  
end subroutine output
  
  
  
  
  
  
  
  
  
  
integer function newunit(unit) result(n)
  ! returns lowest i/o unit number not in use
  integer, intent(out), optional :: unit
  logical inuse
  integer, parameter :: nmin=101   ! avoid lower numbers which are sometimes reserved
  integer, parameter :: nmax=999  ! may be system-dependent
  do n = nmin, nmax
      inquire(unit=n, opened=inuse)
      if (.not. inuse) then
          if (present(unit)) unit=n
          return
      end if
  end do
  write(msg_file,*)'newunit error: available unit not found.'
  call exit_function
end function




end module output_module
