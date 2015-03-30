      module output_module
      use parameter_module
      use lib_node_module
      use lib_elem_module
      use integration_point_module
      
      implicit none
      private
      
      ! define global variable for output directory
      character(len=dirlength), save :: outdir
      
      
      public :: outdir, output
      
      
      contains
      
      
      subroutine output(kstep,kinc,outdir)
         
        ! passed-in variables
        integer,intent(in)          :: kstep    ! current step number
        integer,intent(in)          :: kinc     ! increment number of current step
        character(len=dirlength),intent(in) :: outdir   ! output directory name
        
        ! output file variables
        integer                     :: outunit  ! output file unit
        character(len=dirlength)    :: outfile  ! output file name
        character(len=dirlength)    :: outnum   ! output increment number (embedded in the outfile name)
        character(len=10)           :: fmat, fmatcnc     ! format specs
      
        ! no. of nodes & elem in the mesh
        integer                     :: nnode, nelem
        
        ! no. of elems of each elem type
        integer                     :: nsub3d, nsubwedge, nsubbrick, nsubcoh3d6, nsubcoh3d8
        integer                     :: nxbrick
        integer                     :: nxlam,nplyblk,ninterf
        
        ! sub element type variable
        character(len=eltypelength) :: subtype
        
        ! nsize: size of vtk element output; elnode: no. nodes in each elem
        integer                     :: nsize, elnode

        ! connectivity of an elem extracted from lib_elem
        integer, allocatable        :: connec(:)
        
        ! nodal coordinates
        real(kind=dp), allocatable  :: x(:)     ! coordinates of nodes extracted from lib_node
        real(kind=dp)               :: x3d(3)   ! coordinates of a node in 3D
        
        ! nodal displacements
        real(kind=dp), allocatable  :: disp(:)  ! displacements of nodes extracted from lib_node
        real(kind=dp)               :: disp3d(3)! displacements of a node in 3D
        
        
        ! sub element arrays
        type(wedge_element),allocatable :: subwedge(:)
        type(brick_element),allocatable :: subbrick(:)
        type(coh3d6_element),allocatable :: subcoh3d6(:)
        type(coh3d8_element),allocatable :: subcoh3d8(:)
        type(sub3d_element),allocatable :: sub3d(:)
        type(xbrick_element),allocatable :: subplyblk(:)
        !type(coh3d8_element),allocatable :: subinterf(:)
        type(xcoh_element),allocatable :: subinterf(:)
        type(coh3d8_element),allocatable :: mainelem(:)
        type(subxcoh_element),allocatable :: subelem(:)
        
        ! integration points
        type(integration_point), allocatable :: igpnt(:) ! intg point array
        
        ! stress and strain
        real(kind=dp), allocatable  :: sig(:), eps(:) ! stress & strain arrays extracted from lib_elem ig pnt
        real(kind=dp)               :: sigtsr(3,3), epstsr(3,3) ! stress & strain tensors for vtk output
        
        ! damage/failure variables
        real(kind=dp)               :: fvar    ! temporary failure variable for scalar outpt
        type(sdv_array),allocatable :: fsdv(:)  ! failure variables extracted from ig point sdv array
        integer                     :: ifvar,ifvar1,ifvar2
        
        ! counters
        integer                     :: i, j, l, m, n, nxl

      
      
      
      




      
        ! -----------------------------------------------------------------!
        !                       initialize variables
        ! -----------------------------------------------------------------!
        
        outunit=0; outfile=''; outnum=''; fmat=''; subtype=''
        nnode=0; nelem=0
        nsub3d=0; nsubwedge=0; nsubbrick=0; nsubcoh3d6=0; nsubcoh3d8=0
        nxbrick=0; nxlam=0; nplyblk=0; ninterf=0
        nsize=0; elnode=0
        x3d=zero; disp3d=zero
        sigtsr=zero; epstsr=zero
        fvar=zero; ifvar=0; ifvar1=0; ifvar2=0
        i=0; j=0; l=0; m=0; n=0; nxl=0
        
        ! obtain nnode value from glb libraries
        nnode=size(lib_node)
      
        ! set format for integer output
        fmat='(i5.5)'   ! for increment no.
        fmatcnc='(i10)'  ! for connec node no.
      
        ! write the increment number as a character and store in outnum
        write(outnum,fmat) kinc
        
        
        ! create the output file name
        !~outfile=trim(outdir)//'/outputs/'//trim(outnum)//'.vtk'
        outfile=trim(outdir)//'fnm-'//trim(outnum)//'.vtk'
        outfile=trim(outfile)
        
        ! open the outfile
        open(newunit(outunit), file=outfile,status="replace",action="write")
        
        
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
        
        write(outunit,'(a, i10, a)')'POINTS ',nnode,' FLOAT'    
        do i=1,nnode
            ! empty x3d array
            x3d=zero
            ! extract nodal coords from lib_node
            call extract(lib_node(i),x=x)
            ! pass nodal coords to x3d array
            if(allocated(x)) then 
                x3d(1:size(x))=x(:)
                deallocate(x)
            end if
            ! write x3d array into output file
            write(outunit,*) x3d(1),x3d(2),x3d(3)
        end do            
        write(outunit,'(a)')''
        
        
        
        





        
        ! -----------------------------------------------------------------!
        !                     write elements 
        !       (order matters: new elem type must join the queue)
        ! -----------------------------------------------------------------!
        
        
        !                   write elements' connec (order matters)
        
        ! get no. of elems of each type
              
        if(allocated(lib_xlam)) then
            nxlam=size(lib_xlam)
            do nxl=1,nxlam
                call extract(lib_xlam(nxl),plyblk=subplyblk,interf=subinterf)
                
                if(allocated(subplyblk)) then
                    nplyblk=size(subplyblk)
                    do m=1,nplyblk
                        call extract(subplyblk(m),subelem=sub3d)
                        nsub3d=size(sub3d)
                        do j=1,nsub3d
                            call extract(sub3d(j),eltype=subtype)
                            select case(subtype)
                                case('wedge')
                                    nsubwedge=nsubwedge+1
                                case('brick')
                                    nsubbrick=nsubbrick+1
                                case('coh3d6')
                                    nsubcoh3d6=nsubcoh3d6+1
                                case('coh3d8')
                                    nsubcoh3d8=nsubcoh3d8+1
                                case default
                                    continue
                            end select
                        end do
                        deallocate(sub3d)
                    end do
                    deallocate(subplyblk)
                end if
                
                if(allocated(subinterf)) then
                    ninterf=ninterf+size(subinterf)
                    deallocate(subinterf)
                end if
                
            end do
        end if
        ! .... and other elem types ....
        
        ! total no. of elems
        nelem=nsubwedge+nsubbrick+nsubcoh3d6+nsubcoh3d8+ninterf
        
        ! calculate total no. of nodes to print; each row has 1+elnode no. of indices to print
        nsize=nsubwedge*(1+6)+nsubbrick*(1+8)+nsubcoh3d6*(1+6)+nsubcoh3d8*(1+8)+ninterf*(1+8)
        
        ! write a summary of output
        write(outunit,'(a, i10, i10)')'CELLS ', nelem, nsize
        
        
        ! ------ write x elements ---------
        
        
        if(nxlam > 0) then
            ! loop over all xlam elems
            do nxl=1,nxlam
                call extract(lib_xlam(nxl),plyblk=subplyblk,interf=subinterf)
                
                
                ! write nodal cnc of plyblk elems
                
                if(allocated(subplyblk)) then
                    nplyblk=size(subplyblk)                   
                    ! loop over all plyblks (type xbrick)
                    do m=1,nplyblk
                        call extract(subplyblk(m),subelem=sub3d)
                        nsub3d=size(sub3d)
                        
                        if(nsub3d > 0) then                        
                            ! loop over all sub elems of xbrick and write each sub elem's connec individually
                            do i=1,nsub3d
                                ! extract connec
                                call extract(sub3d(i),glbcnc=connec) 
                                ! print connec in vtk; note that in vtk node no. starts from 0
                                connec=connec-1
                                write(outunit,fmatcnc,advance="no") size(connec)
                                do j=1,size(connec)
                                    write(outunit,fmatcnc,advance="no") connec(j)
                                end do
                                write(outunit,'(a)')''
                                deallocate(connec)
                            end do                       
                        end if
                        
                        deallocate(sub3d)
                    end do 
                    deallocate(subplyblk)
                end if


                ! write nodal cnc of interface elems
                
                if(allocated(subinterf)) then
                    ninterf=size(subinterf)                
                    if(ninterf > 0) then
                        do i=1,ninterf
                            call extract(subinterf(i),nodecnc=connec)
                            connec=connec-1
                            write(outunit,fmatcnc,advance="no") 8   ! only write real nodes
                            do j=1,8
                                write(outunit,fmatcnc,advance="no") connec(j)
                            end do
                            write(outunit,'(a)')''
                            deallocate(connec)
                        end do
                    end if 
                    deallocate(subinterf)
                end if
            
                
                
            
            ! loop over all xlam elems
            end do
            
        end if
        
        write(outunit,'(a)')''


         
        
        !                   write elements' types (order matters)
        
        write(outunit,'(a, i10)')'CELL_TYPES ', nelem
            

        if(nxlam > 0) then
            ! loop over all xlam elems
            do nxl=1,nxlam
                call extract(lib_xlam(nxl),plyblk=subplyblk,interf=subinterf)
                
                
                ! write typekeys of plyblk elems
                
                if(allocated(subplyblk)) then
                    nplyblk=size(subplyblk)                   
                    ! loop over all plyblks (type xbrick)
                    do m=1,nplyblk
                        call extract(subplyblk(m),subelem=sub3d)
                        nsub3d=size(sub3d)
                        
                        if(nsub3d > 0) then                        
                            ! loop over all sub elems of xbrick and write each sub elem's typekey individually
                            do i=1,nsub3d
                                call extract(sub3d(i),eltype=subtype)
                                select case(subtype)
                                    case('wedge')
                                        write(outunit,'(i2)') 13 ! 13 for wedge
                                    case('brick')
                                        write(outunit,'(i2)') 12 ! 12 for brick
                                    case('coh3d6')
                                        write(outunit,'(i2)') 13 ! 13 for coh3d6
                                    case('coh3d8')
                                        write(outunit,'(i2)') 12 ! 12 for coh3d8
                                    case default
                                        continue
                                end select
                            end do                       
                        end if
                        
                        deallocate(sub3d)
                    end do 
                    deallocate(subplyblk)
                end if


                ! write typekeys of interface elems
                
                if(allocated(subinterf)) then
                    ninterf=size(subinterf)                
                    if(ninterf > 0) then
                        do i=1,ninterf
                            write(outunit,'(i2)') 12 ! 12 for coh3d8
                        end do
                    end if   
                    deallocate(subinterf)
                end if
            
                
                
            
            ! loop over all xlam elems
            end do
            
        end if
       
        write(outunit,'(a)')''










        ! -----------------------------------------------------------------!
        !                     write displacements
        ! -----------------------------------------------------------------!


        ! write nodal varibales (disp.)
        write(outunit,'(a, i10)')'POINT_DATA ', nnode
        write(outunit,'(a)')'VECTORS displacement float'
        !~write(outunit,'(a)')'LOOKUP_TABLE default' ! ** this is only for scalar **
        
        do i=1, nnode
            disp3d=zero ! empty disp3d for reuse
            call extract(lib_node(i),u=disp)
            if(allocated(disp)) disp3d(1:size(disp))=disp(:)
            write(outunit,*) disp3d(1),disp3d(2),disp3d(3)
            deallocate(disp)
        end do            
        write(outunit,'(a)')''
        
        
        
        
   





        ! -----------------------------------------------------------------!
        !                     write stress (order matters)
        ! -----------------------------------------------------------------!     
        
        ! write element stress
        write(outunit,'(a, i10)')'CELL_DATA ', nelem
        write(outunit,'(a)')'TENSORS stress float'



        if(nxlam > 0) then
            ! loop over all xlam elems
            do nxl=1,nxlam
                call extract(lib_xlam(nxl),plyblk=subplyblk,interf=subinterf)
                
                
                ! write stress tensors of plyblk elems
                
                if(allocated(subplyblk)) then
                    nplyblk=size(subplyblk)                   
                    ! loop over all plyblks (type xbrick)
                    do m=1,nplyblk
                        call extract(subplyblk(m),subelem=sub3d)
                        nsub3d=size(sub3d)
                        
                        if(nsub3d > 0) then                        
                            ! loop over all sub elems of xbrick and write each sub elem's stress individually
                            do i=1,nsub3d
                                call extract(sub3d(i),eltype=subtype)
                                select case(subtype)                
                                    case('wedge')
                                        ! extract wedge sub elem
                                        call extract(sub3d(i),wedge=subwedge)
                                        sigtsr=zero ! empty sig & eps tensor for reuse  
                                        call extract(subwedge(1),stress=sig)
                                        sigtsr(1,1)=sigtsr(1,1)+sig(1)
                                        sigtsr(2,2)=sigtsr(2,2)+sig(2)
                                        sigtsr(3,3)=sigtsr(3,3)+sig(3)
                                        sigtsr(1,2)=sigtsr(1,2)+sig(4)
                                        sigtsr(1,3)=sigtsr(1,3)+sig(5)
                                        sigtsr(2,3)=sigtsr(2,3)+sig(6)
                                        sigtsr(2,1)=sigtsr(2,1)+sig(4)
                                        sigtsr(3,1)=sigtsr(3,1)+sig(5)
                                        sigtsr(3,2)=sigtsr(3,2)+sig(6)
                                        do l=1,3
                                            write(outunit,*) sigtsr(1,l), sigtsr(2,l), sigtsr(3,l)
                                        end do
                                        write(outunit,'(a)')'' ! separate from next element
                                        
                                        deallocate(sig)
                                        deallocate(subwedge)
                                        
                                    case('brick')
                                        ! extract brick sub elem
                                        call extract(sub3d(i),brick=subbrick)
                                        sigtsr=zero ! empty sig & eps tensor for reuse 
                                        call extract(subbrick(1),stress=sig)
                                        sigtsr(1,1)=sigtsr(1,1)+sig(1)
                                        sigtsr(2,2)=sigtsr(2,2)+sig(2)
                                        sigtsr(3,3)=sigtsr(3,3)+sig(3)
                                        sigtsr(1,2)=sigtsr(1,2)+sig(4)
                                        sigtsr(1,3)=sigtsr(1,3)+sig(5)
                                        sigtsr(2,3)=sigtsr(2,3)+sig(6)
                                        sigtsr(2,1)=sigtsr(2,1)+sig(4)
                                        sigtsr(3,1)=sigtsr(3,1)+sig(5)
                                        sigtsr(3,2)=sigtsr(3,2)+sig(6)
                                        do l=1,3
                                            write(outunit,*) sigtsr(1,l), sigtsr(2,l), sigtsr(3,l)
                                        end do
                                        write(outunit,'(a)')'' ! separate from next element 
                                        
                                        deallocate(sig)
                                        deallocate(subbrick)
                    
                                    case('coh3d6')
                                        ! extract coh3d6 sub elem
                                        !call extract(sub3d(i),coh3d6=subcoh3d6)
                                        sigtsr=zero ! empty sig & eps tensor for reuse
                                        do l=1,3
                                            write(outunit,*) sigtsr(1,l), sigtsr(2,l), sigtsr(3,l)
                                        end do
                                        write(outunit,'(a)')'' ! separate from next element 
                                        
                                        !deallocate(subcoh3d6)
                                        
                                    case('coh3d8')
                                        ! extract coh3d8 sub elem
                                        !call extract(sub3d(i),coh3d8=subcoh3d8)
                                        sigtsr=zero ! empty sig & eps tensor for reuse
                                        do l=1,3
                                            write(outunit,*) sigtsr(1,l), sigtsr(2,l), sigtsr(3,l)
                                        end do
                                        write(outunit,'(a)')'' ! separate from next element
                                        
                                        !deallocate(subcoh3d8)
                   
                                    case default
                                        continue                       
                                end select    
                            end do                       
                        end if
                        
                        deallocate(sub3d)
                    end do  
                    deallocate(subplyblk)
                end if


                ! write stress tensors of interface elems (all zero; coh elems only output failure)
                
                if(allocated(subinterf)) then
                    ninterf=size(subinterf)                
                    if(ninterf > 0) then
                        do i=1,ninterf
                            sigtsr=zero ! empty sig & eps tensor for reuse
                            do l=1,3
                                write(outunit,*) sigtsr(1,l), sigtsr(2,l), sigtsr(3,l)
                            end do
                            write(outunit,'(a)')'' ! separate from next element
                        end do
                    end if 
                    deallocate(subinterf)
                end if
            
                
                
            
            ! loop over all xlam elems
            end do
            
        end if
        
        
               
        
        
        ! -----------------------------------------------------------------!
        !                     write strain (order matters)
        ! -----------------------------------------------------------------!  
        
        ! write element strain
        write(outunit,'(a)')'TENSORS strain float'





        if(nxlam > 0) then
            ! loop over all xlam elems
            do nxl=1,nxlam
                call extract(lib_xlam(nxl),plyblk=subplyblk,interf=subinterf)
                
                
                ! write strain tensors of plyblk elems
                
                if(allocated(subplyblk)) then
                    nplyblk=size(subplyblk)                   
                    ! loop over all plyblks (type xbrick)
                    do m=1,nplyblk
                        call extract(subplyblk(m),subelem=sub3d)
                        nsub3d=size(sub3d)
                        
                        if(nsub3d > 0) then                        
                            ! loop over all sub elems of xbrick and write each sub elem's strains individually
                            do i=1,nsub3d
                                call extract(sub3d(i),eltype=subtype)
                                select case(subtype)                
                                    case('wedge')
                                        ! extract wedge sub elem
                                        call extract(sub3d(i),wedge=subwedge)
                                        epstsr=zero ! empty eps & eps tensor for reuse  
                                        call extract(subwedge(1),strain=eps)
                                        epstsr(1,1)=epstsr(1,1)+eps(1)
                                        epstsr(2,2)=epstsr(2,2)+eps(2)
                                        epstsr(3,3)=epstsr(3,3)+eps(3)
                                        epstsr(1,2)=epstsr(1,2)+eps(4)
                                        epstsr(1,3)=epstsr(1,3)+eps(5)
                                        epstsr(2,3)=epstsr(2,3)+eps(6)
                                        epstsr(2,1)=epstsr(2,1)+eps(4)
                                        epstsr(3,1)=epstsr(3,1)+eps(5)
                                        epstsr(3,2)=epstsr(3,2)+eps(6)
                                        do l=1,3
                                            write(outunit,*) epstsr(1,l), epstsr(2,l), epstsr(3,l)
                                        end do
                                        write(outunit,'(a)')'' ! separate from next element
                                        
                                        deallocate(subwedge)
                                        deallocate(eps)
                                        
                                    case('brick')
                                        ! extract brick sub elem
                                        call extract(sub3d(i),brick=subbrick)
                                        epstsr=zero ! empty eps & eps tensor for reuse  
                                        call extract(subbrick(1),strain=eps)
                                        epstsr(1,1)=epstsr(1,1)+eps(1)
                                        epstsr(2,2)=epstsr(2,2)+eps(2)
                                        epstsr(3,3)=epstsr(3,3)+eps(3)
                                        epstsr(1,2)=epstsr(1,2)+eps(4)
                                        epstsr(1,3)=epstsr(1,3)+eps(5)
                                        epstsr(2,3)=epstsr(2,3)+eps(6)
                                        epstsr(2,1)=epstsr(2,1)+eps(4)
                                        epstsr(3,1)=epstsr(3,1)+eps(5)
                                        epstsr(3,2)=epstsr(3,2)+eps(6)
                                        do l=1,3
                                            write(outunit,*) epstsr(1,l), epstsr(2,l), epstsr(3,l)
                                        end do
                                        write(outunit,'(a)')'' ! separate from next element 
                                        
                                        deallocate(subbrick)
                                        deallocate(eps)
                    
                                    case('coh3d6')
                                        ! extract coh3d6 sub elem
                                        !call extract(sub3d(i),coh3d6=subcoh3d6)
                                        epstsr=zero ! empty eps & eps tensor for reuse
                                        do l=1,3
                                            write(outunit,*) epstsr(1,l), epstsr(2,l), epstsr(3,l)
                                        end do
                                        write(outunit,'(a)')'' ! separate from next element 
                                        
                                        !deallocate(subcoh3d6)
                                        
                                    case('coh3d8')
                                        ! extract coh3d8 sub elem
                                        !call extract(sub3d(i),coh3d8=subcoh3d8)
                                        epstsr=zero ! empty eps & eps tensor for reuse
                                        do l=1,3
                                            write(outunit,*) epstsr(1,l), epstsr(2,l), epstsr(3,l)
                                        end do
                                        write(outunit,'(a)')'' ! separate from next element
                                        
                                        !deallocate(subcoh3d8)
                   
                                    case default
                                        continue                       
                                end select    
                            end do                       
                        end if
                        
                        deallocate(sub3d)
                    end do
                    deallocate(subplyblk)
                end if


                ! write strain tensors of interface elems (all zero; coh elems only output failure)
                
                if(allocated(subinterf)) then
                    ninterf=size(subinterf)                
                    if(ninterf > 0) then
                        do i=1,ninterf
                            epstsr=zero ! empty eps & eps tensor for reuse
                            do l=1,3
                                write(outunit,*) epstsr(1,l), epstsr(2,l), epstsr(3,l)
                            end do
                            write(outunit,'(a)')'' ! separate from next element
                        end do
                    end if  
                    deallocate(subinterf)
                end if
                
                
            
            ! loop over all xlam elems
            end do
            
        end if
        
       
       
       
       

        ! -----------------------------------------------------------------!
        !                     write failure status
        ! -----------------------------------------------------------------! 
        
        ! write element failure status
        write(outunit,'(a)')'SCALARS fstat float'
        write(outunit,'(a)')'LOOKUP_TABLE default'
        


        if(nxlam > 0) then
            ! loop over all xlam elems
            do nxl=1,nxlam
                call extract(lib_xlam(nxl),plyblk=subplyblk,interf=subinterf)
                
                
                ! write failure status of plyblk elems
                
                if(allocated(subplyblk)) then
                    nplyblk=size(subplyblk)                   
                    ! loop over all plyblks (type xbrick)
                    do m=1,nplyblk
                        call extract(subplyblk(m),subelem=sub3d)
                        nsub3d=size(sub3d)
                                               
                        ! loop over all sub elems of xbrick and write each sub elem's failure status individually
                        do i=1,nsub3d
                            call extract(sub3d(i),curr_status=ifvar)
                            write(outunit,*) ifvar  
                        end do
                        
                        deallocate(sub3d)
                    end do  

                    deallocate(subplyblk)
                end if


                ! write failure status of interface elems
                
                if(allocated(subinterf)) then
                    ninterf=size(subinterf)    
            
                    ! loop over all xcoh elems
                    do i=1,ninterf
                        call extract(subinterf(i),curr_status=ifvar)
                        write(outunit,*) ifvar 
                    end do

                    deallocate(subinterf)
                end if
                
            
            ! loop over all xlam elems
            end do
            
        end if
        


        ! -----------------------------------------------------------------!
        !                     write damage variable
        ! -----------------------------------------------------------------! 
        
        ! write element damage variable
        write(outunit,'(a)')'SCALARS dm float'
        write(outunit,'(a)')'LOOKUP_TABLE default'
        


        if(nxlam > 0) then
            ! loop over all xlam elems
            do nxl=1,nxlam
                call extract(lib_xlam(nxl),plyblk=subplyblk,interf=subinterf)
                
                
                ! write damage variable of plyblk elems
                
                if(allocated(subplyblk)) then
                    nplyblk=size(subplyblk)                   
                    ! loop over all plyblks (type xbrick)
                    do m=1,nplyblk
                        call extract(subplyblk(m),subelem=sub3d)
                        nsub3d=size(sub3d)
                                               
                        ! loop over all sub elems of xbrick and write each sub elem's damage variable individually
                    
                        do i=1,nsub3d
                            call extract(sub3d(i),eltype=subtype)
                            select case(subtype)
                                case('wedge')
                                    call extract(sub3d(i),wedge=subwedge)
                                    call extract(subwedge(1),ig_point=igpnt)
                                    fvar=zero
                                    do j=1,size(igpnt)
                                        call extract(igpnt(j),sdv=fsdv)
                                        if(allocated(fsdv).and.allocated(fsdv(2)%r)) fvar=fvar+fsdv(2)%r(1)
                                        if(allocated(fsdv)) deallocate(fsdv)
                                    end do 
                                    ! average damage variable in the element
                                    fvar=fvar/size(igpnt)
                                    write(outunit,*) fvar
                                    
                                    deallocate(igpnt)
                                    deallocate(subwedge)
                                    
                                case('brick')
                                    call extract(sub3d(i),brick=subbrick)
                                    call extract(subbrick(1),ig_point=igpnt)
                                    fvar=zero
                                    do j=1,size(igpnt)
                                        call extract(igpnt(j),sdv=fsdv)
                                        if(allocated(fsdv).and.allocated(fsdv(2)%r)) fvar=fvar+fsdv(2)%r(1)
                                        if(allocated(fsdv)) deallocate(fsdv)
                                    end do 
                                    ! average damage variable in the element
                                    fvar=fvar/size(igpnt)
                                    write(outunit,*) fvar
                                    
                                    deallocate(igpnt)
                                    deallocate(subbrick)
                                    
                                case('coh3d6')
                                    call extract(sub3d(i),coh3d6=subcoh3d6)
                                    call extract(subcoh3d6(1),ig_point=igpnt)
                                    fvar=zero
                                    do j=1,size(igpnt)
                                        call extract(igpnt(j),sdv=fsdv)
                                        if(allocated(fsdv).and.allocated(fsdv(2)%r)) fvar=fvar+fsdv(2)%r(1)
                                        if(allocated(fsdv)) deallocate(fsdv)
                                    end do 
                                    ! average damage variablen in the element
                                    fvar=fvar/size(igpnt)
                                    write(outunit,*) fvar
                                    
                                    deallocate(igpnt)
                                    deallocate(subcoh3d6)
                                    
                                case('coh3d8')
                                    call extract(sub3d(i),coh3d8=subcoh3d8)
                                    call extract(subcoh3d8(1),ig_point=igpnt)
                                    fvar=zero
                                    do j=1,size(igpnt)
                                        call extract(igpnt(j),sdv=fsdv)
                                        if(allocated(fsdv).and.allocated(fsdv(2)%r)) fvar=fvar+fsdv(2)%r(1)
                                        if(allocated(fsdv)) deallocate(fsdv)
                                    end do 
                                    ! average damage variable in the element
                                    fvar=fvar/size(igpnt)
                                    write(outunit,*) fvar
                                    
                                    deallocate(igpnt)
                                    deallocate(subcoh3d8)
                                    
                                case default
                                    continue
                            end select    
                        end do
     
                        
                        deallocate(sub3d)
                    end do  

                    deallocate(subplyblk)
                end if


                ! write damage variable of interface elems
                
                if(allocated(subinterf)) then
                    ninterf=size(subinterf)                

                    do i=1,ninterf
                        
                        fvar=zero
                        
                        call extract(subinterf(i),mainelem=mainelem,subelem=subelem)
                        
                        if(allocated(mainelem)) then
                        ! a coh3d8 elem
                            call extract(mainelem(1),ig_point=igpnt)
                            
                            do j=1,size(igpnt)
                                call extract(igpnt(j),sdv=fsdv)
                                if(allocated(fsdv).and.allocated(fsdv(2)%r)) fvar=fvar+fsdv(2)%r(1)
                                if(allocated(fsdv)) deallocate(fsdv)
                            end do 
                            ! average failure status in the element
                            fvar=fvar/size(igpnt)
                            write(outunit,*) fvar
                            
                            deallocate(igpnt)
                            deallocate(mainelem)
                        else 
                        ! 
                            if(.not.allocated(subelem)) then
                                write(msg_file,*)'subelem not allocated for output!'
                                call exit_function
                            end if
                            
                            ! extract sdv of the 2 subxcoh elems, and update to ifvar
                            do j=1, 2
                                call extract(subelem(j),sdv=fsdv)
                                if(allocated(fsdv)) then 
                                    if(allocated(fsdv(1)%r)) fvar=max(fvar,fsdv(1)%r(1))
                                    deallocate(fsdv)
                                end if
                            end do
                            
                            write(outunit,*) fvar
                            
                            deallocate(subelem)
                        
                        end if
                        
                        
                    end do

                    deallocate(subinterf)
                end if
                
            
            ! loop over all xlam elems
            end do
            
        end if
        


        


        if(allocated(connec)) deallocate(connec)
        if(allocated(x)) deallocate(x)
        if(allocated(disp)) deallocate(disp)
        if(allocated(subwedge)) deallocate(subwedge)
        if(allocated(subbrick)) deallocate(subbrick)
        if(allocated(subcoh3d6)) deallocate(subcoh3d6)
        if(allocated(subcoh3d8)) deallocate(subcoh3d8)
        if(allocated(sub3d)) deallocate(sub3d)
        if(allocated(subplyblk)) deallocate(subplyblk)
        if(allocated(subinterf)) deallocate(subinterf)
        if(allocated(mainelem)) deallocate(mainelem)
        if(allocated(subelem)) deallocate(subelem)
        if(allocated(igpnt)) deallocate(igpnt)
        if(allocated(sig)) deallocate(sig)
        if(allocated(eps)) deallocate(eps)   
        if(allocated(fsdv)) deallocate(fsdv)        
        
        
        close(outunit)
        
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
