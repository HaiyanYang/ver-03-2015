module output_module
use parameter_module, only: DIRLENGTH, ELTYPELENGTH, DP, ZERO, MSG_FILE, EXIT_FUNCTION

implicit none
private

! define global variable for output directory
character(len=DIRLENGTH), save :: outdir


public :: outdir, output


contains


subroutine output(kstep,kinc,outdir)
use node_list_module, only: node_list
use elem_list_module, only: elem_list, elem_node_connec
use fnode_module,     only: extract
use fBrickLam_elem_module, only: extract
   
  ! passed-in variables
  integer,                  intent(in)  :: kstep    ! current step number
  integer,                  intent(in)  :: kinc     ! increment number of current step
  character(len=DIRLENGTH), intent(in)  :: outdir   ! output directory name
  
  ! local variables
  
  ! private parameter, length of format string
  integer, parameter          :: FMTLENGTH = 10

  
  ! output file variables
  integer                     :: outunit  ! output file unit
  character(len=DIRLENGTH)    :: outfile  ! output file name
  character(len=DIRLENGTH)    :: outnum   ! output increment number (embedded in the outfile name)
  character(len=FMTLENGTH)    :: FMATKINC, FMATNNODE, FMATNELEM, FMATFSTAT, FMATFLOAT

  ! no. of nodes & elem in the mesh
  integer                     :: nnode, nelem
  
  ! no. of elems of each elem type
  integer                     :: nwedge, nbrick
  
  ! sub element type variable
  character(len=ELTYPELENGTH) :: eltype
  
  ! nsize: size of vtk element output; elnode: no. nodes in each elem
  integer                     :: nsize
  
  ! nodal coordinates
  real(kind=DP), allocatable  :: x(:)     ! coordinates of nodes extracted from lib_node
  
  ! nodal displacements
  real(kind=DP), allocatable  :: disp(:)  ! displacements of nodes extracted from lib_node
  
  ! counters
  integer  :: i



  ! -----------------------------------------------------------------!
  !                       initialize variables
  ! -----------------------------------------------------------------!
  
  outunit   = 0
  outfile   = ''
  outnum    = ''
  FMATKINC  = ''
  FMATNNODE = ''
  FMATNELEM = ''
  FMATFLOAT = ''
  eltype    = ''
  nnode     = 0
  nwedge    = 0
  nbrick    = 0
  nelem     = 0
  nsize     = 0

  i=0

  ! format parameters
  FMATKINC  = 'i5.5' ! for increment no.
  FMATNNODE = 'i10'  ! for node no.
  FMATNELEM = 'i10'  ! for elem no.
  FMATFSTAT = 'i2'   ! for fstat variable format (2 digits would suffice)
  FMATFLOAT = 'ES20.3E9' ! scientific notation, repeat=1, min. width=10, digits=3

  ! write the increment number as a character and store in outnum
  write(outnum,'('//trim(FMATKINC)//')') kinc
  
  ! create the output file name
  !~outfile=trim(outdir)//'/outputs/'//trim(outnum)//'.vtk'
  outfile = trim(outdir)//'fnm-'//trim(outnum)//'.vtk'
  outfile = trim(outfile)
  
  ! open the outfile
  open(newunit(outunit), file=outfile, status="replace", action="write")
  
  ! write header
  write(outunit,'(a)')'# vtk DataFile Version 3.1'
  write(outunit,'(a)')'for Floating Node Method output'
  
  ! write vtk format
  write(outunit,'(a)')'ASCII'      
  
  ! write vtk data type
  write(outunit,'(a)')'DATASET UNSTRUCTURED_GRID'


  ! -----------------------------------------------------------------!
  !                     write nodes
  ! -----------------------------------------------------------------!
  ! obtain nnode value
  nnode = size(node_list)
  ! write node summary line
  write(outunit,'(a, '//trim(FMATNNODE)//', a)')'POINTS ',nnode,' FLOAT'
  ! write all nodes
  do i = 1,nnode
      ! extract nodal coords from lib_node
      call extract(node_list(i),x=x)
      ! write x3d array into output file
      write(outunit,'(3'//trim(FMATFLOAT)//')') x(1),x(2),x(3)
  end do            
  write(outunit,'(a)')''

  
  ! -----------------------------------------------------------------!
  !                     write element connecs
  ! -----------------------------------------------------------------!  
  ! wflag = 'nelem' would allow the wfLam subroutine to update
  ! nwedge, nbrick, ncoh6 and ncoh8
  call wfLam('nelem')

  ! total no. of elems
  nelem = nwedge + nbrick
  
  ! calculate total no. of nodes to print; each row has 1+elnode no. of indices to print
  nsize = nwedge * 7 + nbrick * 9
  
  ! write a summary of output
  write(outunit,'(a, 2'//trim(FMATNNODE)//')')'CELLS ', nelem, nsize
  
  ! wflag = 'connec' would allow the wfLam subroutine to write out all
  ! the nodal connec of all sub elems in the mesh
  call wfLam('connec')


  ! -----------------------------------------------------------------!
  !                     write element types (order matters)
  ! -----------------------------------------------------------------!  
  write(outunit,'(a, '//trim(FMATNELEM)//')')'CELL_TYPES ', nelem
  
  ! wflag = 'eltype' would allow the wfLam subroutine to write out all
  ! the elem type code (used by VTK) of all sub elems in the mesh
  call wfLam('eltype')



  ! -----------------------------------------------------------------!
  ! BEGIN WRITE NODAL DATA
  ! -----------------------------------------------------------------!
  write(outunit,'(a, '//trim(FMATNNODE)//')')'POINT_DATA ', nnode

  
  ! -----------------------------------------------------------------!
  !                     write displacements
  ! -----------------------------------------------------------------!
  write(outunit,'(a)')'VECTORS displacement float'
  
  do i=1, nnode
      call extract(node_list(i),u=disp)
      write(outunit,'(3'//trim(FMATFLOAT)//')') disp(1),disp(2),disp(3)
      deallocate(disp)
  end do  
          
  write(outunit,'(a)')'' 
  
  
  !~! -----------------------------------------------------------------!
  !~!                     write node status
  !~! -----------------------------------------------------------------!
  !~write(outunit,'(a)')'SCALARS nstat int'
  !~write(outunit,'(a)')'LOOKUP_TABLE default'
  !~do i=1, nnode
  !~    call extract(node_list(i),nstat=intvar)
  !~    write(outunit,'('//trim(FMATFSTAT)//')') intvar
  !~end do  
  !~        
  !~write(outunit,'(a)')''
  


  ! -----------------------------------------------------------------!
  ! BEGIN WRITE CELL DATA
  ! -----------------------------------------------------------------!
  write(outunit,'(a, '//trim(FMATNELEM)//')')'CELL_DATA ', nelem
  


  ! -----------------------------------------------------------------!
  !                     write stress (order matters)
  ! -----------------------------------------------------------------!     
  write(outunit,'(a)')'TENSORS stress float'
  
  ! wflag = 'stress' would allow the wfLam subroutine to write out all
  ! stress tensors of all ply sub elems in the mesh
  call wfLam('stress')

  
  ! -----------------------------------------------------------------!
  !                     write strain (order matters)
  ! -----------------------------------------------------------------!  
  write(outunit,'(a)')'TENSORS strain float'
  
  ! wflag = 'strain' would allow the wfLam subroutine to write out all
  ! strain tensors of all ply sub elems in the mesh
  call wfLam('strain')


  
  ! -----------------------------------------------------------------!
  !                     write stress (order matters)
  ! -----------------------------------------------------------------!     
  write(outunit,'(a)')'VECTORS traction float'
  
  ! wflag = 'traction' would allow the wfLam subroutine to write out all
  ! traction vectors of all coh elems in the mesh
  call wfLam('traction')

  
  ! -----------------------------------------------------------------!
  !                     write strain (order matters)
  ! -----------------------------------------------------------------!  
  write(outunit,'(a)')'VECTORS separation float'
  
  ! wflag = 'strain' would allow the wfLam subroutine to write out all
  ! separation vectors of all coh elems in the mesh
  call wfLam('separation')
  

  ! -----------------------------------------------------------------!
  !                     write fibre damage variable
  ! -----------------------------------------------------------------!
  write(outunit,'(a)')'SCALARS df float'
  write(outunit,'(a)')'LOOKUP_TABLE default'
  
  ! wflag = 'df' would allow the wfLam subroutine to write out all
  ! fibre degradation factors of all ply sub elems in the mesh
  call wfLam('df')


  ! -----------------------------------------------------------------!
  !                     write matrix damage variable
  ! -----------------------------------------------------------------! 
  write(outunit,'(a)')'SCALARS dm float'
  write(outunit,'(a)')'LOOKUP_TABLE default'
  
  ! wflag = 'dm' would allow the wfLam subroutine to write out all 
  ! cohesive degradation factors of all cohesive sub elems in the mesh
  call wfLam('dm')
  
  
  ! -----------------------------------------------------------------!
  !                     write delamination damage variable
  ! -----------------------------------------------------------------! 
  write(outunit,'(a)')'SCALARS dd float'
  write(outunit,'(a)')'LOOKUP_TABLE default'
  
  ! wflag = 'dm' would allow the wfLam subroutine to write out all 
  ! cohesive degradation factors of all cohesive sub elems in the mesh
  call wfLam('dd')
  
  
  close(outunit)
  
  return
  
  
  
  contains
  
  
    subroutine wfLam(wflag)
      character(len=*), intent(in)  :: wflag
      
      integer,  allocatable :: elnodes(:,:)
      real(DP), allocatable :: stress(:,:)
      real(DP), allocatable :: strain(:,:)
      real(DP), allocatable :: tau(:,:)
      real(DP), allocatable :: delta(:,:)
      real(DP), allocatable :: df(:)
      real(DP), allocatable :: dm(:)
      real(DP), allocatable :: dd(:)
      
      integer :: nfLam, nfl, nsub, i
      
      nfLam = size(elem_list)
    
      do nfl = 1, nfLam
    
          call extract(elem_list(nfl), elnodes, stress, strain, tau, delta, df, dm, dd)
          
          !*** debug check ***
          if (.not.allocated(elnodes)) then
            write(MSG_FILE,*)'elnodes is NOT allocated in output'
            CALL EXIT_FUNCTION
          end if
          
          if (.not.allocated(stress)) then
            write(MSG_FILE,*)'stress is NOT allocated in output'
            CALL EXIT_FUNCTION
          end if
          
          if (.not.allocated(strain)) then
            write(MSG_FILE,*)'strain is NOT allocated in output'
            CALL EXIT_FUNCTION
          end if
          
          if (.not.allocated(tau)) then
            write(MSG_FILE,*)'tau is NOT allocated in output'
            CALL EXIT_FUNCTION
          end if
          
          if (.not.allocated(delta)) then
            write(MSG_FILE,*)'delta is NOT allocated in output'
            CALL EXIT_FUNCTION
          end if
          
          if (.not.allocated(df)) then
            write(MSG_FILE,*)'df is NOT allocated in output'
            CALL EXIT_FUNCTION
          end if
          
          if (.not.allocated(dm)) then
            write(MSG_FILE,*)'dm is NOT allocated in output'
            CALL EXIT_FUNCTION
          end if
          
          if (.not.allocated(dd)) then
            write(MSG_FILE,*)'dd is NOT allocated in output'
            CALL EXIT_FUNCTION
          end if
          
          nsub = size(elnodes(1,:))
          
          if (nsub < 1) then
            write(MSG_FILE,*)'elnodes size is incorrect in output'
            CALL EXIT_FUNCTION
          end if
          !***********
          
          do i = 1, nsub
            select case(trim(adjustl(wflag)))
            case('nelem')
                if (count(elnodes(:,i)>0) == 8) then
                  nbrick = nbrick + 1
                else if (count(elnodes(:,i)>0) == 6) then
                  nwedge = nwedge + 1
                else
                  write(MSG_FILE,*) 'unexpected no. of nodes in output!'
                  call EXIT_FUNCTION
                end if
            case('connec')
                call wconnec( elem_node_connec( pack(elnodes(:,i),elnodes(:,i)>0) , nfl) )
            case('eltype')
                if (count(elnodes(:,i)>0) == 8) then
                  write(outunit,'(i2)') 12 ! 12 for brick
                else if (count(elnodes(:,i)>0) == 6) then
                  write(outunit,'(i2)') 13 ! 13 for wedge
                else
                  write(MSG_FILE,*) 'unexpected no. of nodes in output!'
                  call EXIT_FUNCTION
                end if
            case('stress')
                call wtensor(stress(:,i))
            case('strain')
                call wtensor(strain(:,i))
            case('traction')
                write(outunit,'(3'//trim(FMATFLOAT)//')') tau(1,i),tau(2,i),tau(3,i)
            case('separation')
                write(outunit,'(3'//trim(FMATFLOAT)//')') delta(1,i),delta(2,i),delta(3,i)
            case('df')
                call wscalar(df(i))
            case('dm')
                call wscalar(dm(i))
            case('dd')
                call wscalar(dd(i))
            end select
          end do
          
      end do
      
      write(outunit,'(a)')''
    
    end subroutine wfLam
  

    subroutine wconnec(connec)
      integer, intent(in) :: connec(:)
      integer :: j
      ! print connec in vtk; note that in vtk node no. starts from 0
      write(outunit,'('//FMATNNODE//')',advance="no") count(connec>0)
      do j = 1, count(connec>0)
          write(outunit,'('//FMATNNODE//')',advance="no") connec(j)-1
      end do
      write(outunit,'(a)')''
    end subroutine wconnec
    
   
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
      do l = 1, 3
          write(outunit,'(3'//trim(FMATFLOAT)//')') tensor(1,l), tensor(2,l), tensor(3,l)
      end do
      write(outunit,'(a)')''
    end subroutine wtensor
    
       
    subroutine wscalar(x)
      !~class(*), intent(in) :: x
      !~select type(x)
      !~type is (real(DP))
      !~  write(outunit,'('//trim(FMATFLOAT)//')') x
      !~type is (integer)
      !~  write(outunit,'('//trim(FMATFSTAT)//')') x
      !~end select
      real(DP), intent(in) :: x
      write(outunit,'('//trim(FMATFLOAT)//')') x
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
  write(MSG_FILE,*)'newunit error: available unit not found.'
  call EXIT_FUNCTION
end function




end module output_module
