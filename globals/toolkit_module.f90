    module toolkit_module
    use parameter_module
    
    implicit none
    
    contains
    
    
!********************************************************************************************
!                           function glb_dee
!               to transform d matrix from local to global coordinates 
!                           (for composite lamina)
!********************************************************************************************
      function glb_dee(dee,theta)
    
      real(kind=dp),intent(in)   :: dee(:,:), theta
      real(kind=dp)              :: glb_dee(size(dee(:,1)),size(dee(1,:)))

      ! local variables
      real(kind=dp) :: c, s
      integer       :: nst, nrow, ncol
      
      glb_dee=zero
      c=zero; s=zero
      nst=0; nrow=0; ncol=0

      c=cos(pi*theta/halfcirc)
      s=sin(pi*theta/halfcirc)
      
      nrow=size(dee(:,1))
      ncol=size(dee(1,:))
      
      if(nrow/=ncol) then
        write(msg_file,*)"error: d matrix must be square!"
        call exit_function
      else
        nst=nrow
      end if
      
      ! d matrix in global coords stored first in local array glb_dee
      if (nst == 3) then
        glb_dee(1,1) = c*c*c*c*dee(1,1) + two*c*c*s*s*(dee(1,2) &
     &            + two*dee(3,3)) + s*s*s*s*dee(2,2)
        glb_dee(1,2) = s*s*c*c*(dee(1,1) + dee(2,2) - four*dee(3,3)) &
     &            + (s*s*s*s+c*c*c*c)*dee(1,2)
        glb_dee(2,1) = glb_dee(1,2)
        glb_dee(2,2) = s*s*s*s*dee(1,1) + two*c*c*s*s*(dee(1,2) &
     &            + two*dee(3,3)) + c*c*c*c*dee(2,2)
        glb_dee(1,3) = s*c*(c*c*(dee(1,1) - dee(1,2) - two*dee(3,3)) &
     &            + s*s*(dee(1,2) - dee(2,2) + two*dee(3,3)))
        glb_dee(3,1) = glb_dee(1,3)
        glb_dee(2,3) = s*c*(s*s*(dee(1,1) - dee(1,2) - two*dee(3,3)) &
     &            + c*c*(dee(1,2) - dee(2,2) + two*dee(3,3)))
        glb_dee(3,2) = glb_dee(2,3)
        glb_dee(3,3) = c*c*s*s*(dee(1,1)+dee(2,2)-2*dee(1,2)) &
     &			  +(c*c-s*s)**2*dee(3,3)
      else if (nst == 6) then
        glb_dee(1,1) = c*c*c*c*dee(1,1) + two*c*c*s*s*(dee(1,2) &
     &            + two*dee(4,4)) + s*s*s*s*dee(2,2)
        glb_dee(1,2) = s*s*c*c*(dee(1,1) + dee(2,2) - four*dee(4,4)) &
     &            + (s*s*s*s+c*c*c*c)*dee(1,2)
        glb_dee(2,1) = glb_dee(1,2)
        glb_dee(2,2) = s*s*s*s*dee(1,1) + two*c*c*s*s*(dee(1,2) &
     &            + two*dee(4,4)) + c*c*c*c*dee(2,2)
        glb_dee(1,3) = c*c*dee(1,3) + s*s*dee(2,3)
        glb_dee(3,1) = glb_dee(1,3)
        glb_dee(1,4) = s*c*(c*c*(dee(1,1) - dee(1,2) - two*dee(4,4)) &
     &            + s*s*(dee(1,2) - dee(2,2) + two*dee(4,4)))
        glb_dee(4,1) = glb_dee(1,4)
        glb_dee(2,3) = s*s*dee(1,3) + c*c*dee(2,3)
        glb_dee(3,2) = glb_dee(2,3)
        glb_dee(2,4) = s*c*(s*s*(dee(1,1) - dee(1,2) - two*dee(4,4)) &
     &            + c*c*(dee(1,2) - dee(2,2) + two*dee(4,4)))
        glb_dee(4,2) = glb_dee(2,4)
        glb_dee(3,3) = dee(3,3)
        glb_dee(3,4) = c*s*(dee(1,3) - dee(2,3))
        glb_dee(4,3) = glb_dee(3,4)
        glb_dee(5,5) = c*c*dee(5,5) + s*s*dee(6,6)
        glb_dee(5,6) = c*s*(dee(6,6) - dee(5,5))
        glb_dee(6,5) = glb_dee(5,6)
        glb_dee(6,6) = s*s*dee(5,5) + c*c*dee(6,6)
        glb_dee(4,4) = s*s*c*c*(dee(1,1) - two*dee(1,2) + dee(2,2)) &
     &            + (s*s - c*c)*(s*s - c*c)*dee(4,4)
      
      else
       write(msg_file,*) 'error: no. of strains not supported for glb_dee function!'
       call exit_function
      end if

      end function glb_dee
!********************************************************************************************
!********************************************************************************************



!********************************************************************************************
!********************* function determinant **********************************************
!********* returns the determinant of a jacobian matrix	*************************************
!********************************************************************************************
      function determinant(jacob)
!
      real(kind=dp),intent(in) ::   jacob(:,:)
      real(kind=dp)            ::   determinant
      
      integer :: ndim, nshape(2)
      
      determinant=zero
      ndim=0; nshape=0
      
      nshape=shape(jacob)
      
      if(nshape(1)/=nshape(2)) then
        write(msg_file,*)"error: jacobian matrix must be square!"
        call exit_function
      else
        ndim=nshape(1)
      end if
!
      if (ndim .eq. 2) then
        determinant= jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1)
      else if(ndim .eq. 3) then
        determinant= jacob(1,1)*(jacob(2,2)*jacob(3,3)-jacob(3,2)*jacob(2,3))
        determinant= determinant- jacob(1,2)*(jacob(2,1)*jacob(3,3)-jacob(3,1)*jacob(2,3))
        determinant= determinant+ jacob(1,3)*(jacob(2,1)*jacob(3,2)-jacob(3,1)*jacob(2,2))
      else
        write(msg_file,*) "error:dimension not supported for determinant function!"
        call exit_function
      end if
!
      return
      end function determinant
!********************************************************************************************
!********************************************************************************************



!********************************************************************************************    
!********************* function beemat ***************************************************
!******* strain-displacement matrix for infinitesimal deformation ***************************
!********************************************************************************************
! convention of order: eps1, eps2, eps3, gamma12, gamma13, gamma23
!
      function beemat(gn)

      real(kind=dp),intent(in)  :: gn(:,:)
      real(kind=dp) :: beemat(size(gn(1,:))*(size(gn(1,:))+1)/2, size(gn(1,:))*size(gn(:,1)))

      ! local variables  
      real(kind=dp) :: x, y, z
      integer :: nnode, ndim, nst, ndof
      integer :: k,l,m,n
      
      beemat=zero
      x=zero; y=zero; z=zero
      nnode=0; ndim=0; nst=0; ndof=0
      k=0; l=0; m=0; n=0
      
      nnode=size(gn(:,1))
      ndim=size(gn(1,:))
      nst=ndim*(ndim+1)/2

      if (nst == 3) then
        do m=1,nnode
          k= 2*m
          l=k-1
          x=gn(m,1)
          y=gn(m,2)
          beemat(1,l)= x
          beemat(3,k)= x
          beemat(2,k)= y
          beemat(3,l)= y
        end do
      else if (nst == 6) then
        do m=1,nnode
          n= 3*m
          k=n-1
          l=k-1
          x=gn(m,1)
          y=gn(m,2)
          z=gn(m,3)
          beemat(1,l)=x
          beemat(4,k)=x
          beemat(5,n)=x
          beemat(2,k)=y
          beemat(4,l)=y
          beemat(6,n)=y
          beemat(3,n)=z
          beemat(5,l)=z
          beemat(6,k)=z
        end do
      else
       write(msg_file,*) 'no. of strains not supported for kbeemat!'
       call exit_function  
      end if

      return
      end function beemat
!********************************************************************************************
!********************************************************************************************
!
!
!
!********************************************************************************************
!********************* subroutine kjac_inv **************************************************
!************** inverts a jacobian matrix onto itself ***************************************
!********************************************************************************************
      function inverse(jacob,detj)
!
      real(kind=dp),intent(in)          :: jacob(:,:)
      real(kind=dp),optional,intent(in) :: detj  
      real(kind=dp)                     :: inverse(size(jacob(:,1)),size(jacob(1,:)))

      real(kind=dp) :: det
      integer :: ndim, nshape(2)
      
      inverse=zero
      det=zero
      ndim=0; nshape=0
      
      nshape=shape(jacob)
      
      if(nshape(1)/=nshape(2)) then
        write(msg_file,*)"error: jacobian matrix must be square!"
        call exit_function
      else
        ndim=nshape(1)
      end if  

      if(ndim == 2) then
        inverse(1,1)=jacob(2,2)
        inverse(2,1)=-jacob(2,1)
        inverse(1,2)=-jacob(1,2)
        inverse(2,2)=jacob(1,1)
      else if(ndim == 3) then
        inverse(1,1)=jacob(2,2)*jacob(3,3)-jacob(3,2)*jacob(2,3)
        inverse(2,1)=jacob(3,1)*jacob(2,3)-jacob(2,1)*jacob(3,3)
        inverse(3,1)=jacob(2,1)*jacob(3,2)-jacob(3,1)*jacob(2,2)
        inverse(1,2)=jacob(3,2)*jacob(1,3)-jacob(1,2)*jacob(3,3)
        inverse(2,2)=jacob(1,1)*jacob(3,3)-jacob(3,1)*jacob(1,3)
        inverse(3,2)=jacob(3,1)*jacob(1,2)-jacob(1,1)*jacob(3,2)
        inverse(1,3)=jacob(1,2)*jacob(2,3)-jacob(2,2)*jacob(1,3)
        inverse(2,3)=jacob(2,1)*jacob(1,3)-jacob(1,1)*jacob(2,3)
        inverse(3,3)=jacob(1,1)*jacob(2,2)-jacob(2,1)*jacob(1,2)
      else
       write(msg_file,*) "error: dimension not yet supported for inverse!"
       call exit_function
      end if
!
      if(present(detj)) then
        det=detj
      else
        det=determinant(jacob)
      end if
      
      if(det .gt. tiny(one)) then
        inverse=inverse/det
      else
        write(msg_file,*) "error: zero or negative detj in inverse!"
        call exit_function
      end if

      end function inverse
!********************************************************************************************
!********************************************************************************************



!********************************************************************************************
!******************* subroutine ktransfer_strain ********************************************
!********** transfer strains from global to local coordinate systems *********************
!                   (for composite lamina)
!********************************************************************************************
      function lcl_strain(strain,theta)

      real(kind=dp),intent(in)  :: strain(:), theta
      real(kind=dp)             :: lcl_strain(size(strain))

      real(kind=dp) :: c, s, t(size(strain),size(strain))
      integer :: nst
      
      lcl_strain=zero
      c=zero; s=zero; t=zero
      nst=0
      
      nst=size(strain)

      c=cos(pi*theta/halfcirc)
      s=sin(pi*theta/halfcirc)
!
      if (nst == 3) then
        t(1,1)=c*c
        t(1,2)=s*s
        t(1,3)=c*s
        t(2,1)=s*s
        t(2,2)=c*c
        t(2,3)=-c*s
        t(3,1)=-two*c*s
        t(3,2)=two*c*s
        t(3,3)=c*c-s*s
      else if (nst == 6) then
        t(1,1)=c*c
        t(1,2)=s*s
        t(1,4)=c*s
        t(2,1)=s*s
        t(2,2)=c*c
        t(2,4)=-c*s
        t(3,3)=one
        t(4,1)=-two*c*s
        t(4,2)=two*c*s
        t(4,4)=c*c-s*s
        t(5,5)=c
        t(5,6)=-s
        t(6,5)=s
        t(6,6)=c
      else
       write(msg_file,*) 'no. of strains not supported for lcl_strain!'
       call exit_function
      end if
      
      lcl_strain=matmul(t,strain)

      end function lcl_strain
!********************************************************************************************
!********************************************************************************************



!********************************************************************************************
!******************* subroutine normalize ***************************************************
!********** normalize vector and return its magnitude  **************************************
!********************************************************************************************
    subroutine normalize(a,amag)
       
        real(kind=dp), intent(inout) :: a(:)
        real(kind=dp), optional, intent(out)   :: amag
       
        real(dp) :: amag2
        
        if(present(amag)) amag=zero
        amag2=zero
       
        amag2=sqrt(dot_product(a,a))
        
        if(amag2 .gt. tiny(one)) then
            a=a/amag2
        else
            write(msg_file,*) 'warning:zero vector for kunitv!'
!           do nothing
        end if
        
        if(present(amag)) amag=amag2

       return
    end subroutine normalize
!********************************************************************************************
!********************************************************************************************








!********************************************************************************************
!               function to compute the cross-product of two 3D vectors
!********************************************************************************************
    function CrossProduct3D(a, b)
      real(dp), intent(in)  :: a(3), b(3)
      real(dp)              :: CrossProduct3D(3)

      CrossProduct3D(1) = a(2) * b(3) - a(3) * b(2)
      CrossProduct3D(2) = a(3) * b(1) - a(1) * b(3)
      CrossProduct3D(3) = a(1) * b(2) - a(2) * b(1)
    end function CrossProduct3D
!********************************************************************************************
!********************************************************************************************    
    
    
    
    
    
!********************************************************************************************
!               function to compute the distance between two vectors
!********************************************************************************************
    function distance(a, b)
      real(dp), intent(in)  :: a(:), b(:)
      real(dp)              :: distance

      if(size(a)/=size(b)) then
        write(msg_file,*) 'two vector sizes must be the same in distance function!'
        call exit_function
      end if
      
      distance=sqrt(dot_product(a-b,a-b))
      
    end function distance
!********************************************************************************************
!********************************************************************************************    
   



    
!********************************************************************************************
!           subroutine to checks if two lines cross, and locate the crossing point
!********************************************************************************************    
      subroutine klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct,detlc)

      
      ! feed in variables
      real(kind=dp),intent(in) :: x1,y1,x2,y2       ! nodal coords of line 1 (elem edge)
      real(kind=dp),intent(in) :: xp1,yp1,xp2,yp2   ! nodal coords of line 2 (precrack)
      ! update variables
      real(kind=dp),intent(inout) :: xct,yct        ! coords of intersection
      integer,intent(inout) :: iscross              ! status of intersection
      real(kind=dp),optional,intent(out)::detlc
      ! local variables
      real(kind=dp) :: a1,b1,c1,a2,b2,c2,det
      real(kind=dp) :: xmin,ymin,xpmin,ypmin,xmax,ymax,xpmax,ypmax
      ! initialize local/update variables (don't use data in subroutines)     
      a1=zero;b1=zero;c1=zero
      a2=zero;b2=zero;c2=zero
      xmin=zero;ymin=zero;xmax=zero;ymax=zero
      xpmin=zero;ypmin=zero;xpmax=zero;ypmax=zero
      det=zero
!
!     **************** algorithm for the determination of precrack intersecting with element edges: *****************
!     equation of a line: a*x+b*y=c
!     the intersection of line 1 (an edge) and line 2 (a precrack) is to solve the system equation:
!     a1*x+b1*y=c1
!     a2*x+b2*y=c2
!     first, det=a1*b2-a2*b1; 
!     if det=0: lines are parallel, no intersection possible
!     else if det/=0, then the intersection point is:
!       xct=(b2*c1-b1*c2)/det
!       yct=(a1*c2-a2*c1)/det
!       check xct, yct with range of x,y of the edge and of the precrack
!       if (xmin<xct<xmax) or (ymin<yct<ymax): point on edge
!       if (xpmin<xct<xpmax) or (ypmin<yct<ypmax): point on precrack
!       if the point is both on edge and on precrack, then store the coords of this point
!           xp(edge no.)=xct
!           yp(edge no.)=yct
!           fedg(edge no.)=1 !edge crossed by precrack
!      
! line 1 equation
      a1=y2-y1
      b1=x1-x2
      c1=a1*x1+b1*y1
! line 2 equation
      a2=yp2-yp1
      b2=xp1-xp2
      c2=a2*xp1+b2*yp1
! range of line 1
      xmin=min(x1,x2)
      ymin=min(y1,y2)
      xmax=max(x1,x2)
      ymax=max(y1,y2)
! range of line 2
      xpmin=min(xp1,xp2)
      ypmin=min(yp1,yp2)
      xpmax=max(xp1,xp2)
      ypmax=max(yp1,yp2)
      ! check intersection
      det=a1*b2-a2*b1
      if(abs(det)<=tiny(one)) then ! lines are parallel; no intersection possible
!      if(abs(det)==zero) then ! lines are parallel; no intersection possible
        ! do nothing
        
      else ! intersection possible
      
        xct=(b2*c1-b1*c2)/det !-find intersection point
        yct=(a1*c2-a2*c1)/det

        if ((xmin.gt.xct).or.(xmax.lt.xct).or.(ymin.gt.yct).or.(ymax.lt.yct)) then ! intersection not on edge
        ! do nothing
          !write(6,*)'intersection out of range of edge'
            iscross=-1
        else if (xct.eq.x1.and.yct.eq.y1) then
            !write(6,*)'precrack line passes node 1'
            xct=max(xct,xmin+tolerance*(xmax-xmin))
            xct=min(xct,xmax-tolerance*(xmax-xmin))
            yct=max(yct,ymin+tolerance*(ymax-ymin))
            yct=min(yct,ymax-tolerance*(ymax-ymin))
            if ((xpmin.gt.xct).or.(xpmax.lt.xct) .or. (ypmin.gt.yct).or.(ypmax.lt.yct)) then                    
                iscross=21  ! intersection pass node 1 but not on precrack
            else
                iscross=11  ! intersection pass node 1 and lie on precrack
            end if
        else if (xct.eq.x2.and.yct.eq.y2) then
            !write(6,*)'precrack line passes node 2'
            xct=max(xct,xmin+tolerance*(xmax-xmin))
            xct=min(xct,xmax-tolerance*(xmax-xmin))
            yct=max(yct,ymin+tolerance*(ymax-ymin))
            yct=min(yct,ymax-tolerance*(ymax-ymin))
            if ((xpmin.gt.xct).or.(xpmax.lt.xct) .or. (ypmin.gt.yct).or.(ypmax.lt.yct)) then
                iscross=22  ! intersection pass node 2 but not on precrack
            else
                iscross=12  ! intersection pass node 2 and lie on precrack
            end if
        else !-intersection on this edge
            xct=max(xct,xmin+tolerance*(xmax-xmin))
            xct=min(xct,xmax-tolerance*(xmax-xmin))
            yct=max(yct,ymin+tolerance*(ymax-ymin))
            yct=min(yct,ymax-tolerance*(ymax-ymin))
            if ((xpmin.gt.xct).or.(xpmax.lt.xct) .or. (ypmin.gt.yct).or.(ypmax.lt.yct)) then                    
                iscross=2   ! intersection on edge but not on precrack
            else
                iscross=1   ! intersection on edge and also on precrack
            end if
        end if
      end if 

      if(present(detlc)) detlc=det
             
      end subroutine klinecross
      
      
      
      




!********************************************************************************************
!           subroutine to find two crack tips on the edges of an element, passing origin
!******************************************************************************************** 
    subroutine elem_ctips_origin(eltype,theta,coords,edg,nedge,ctip)
        character(len=*),   intent(in)  :: eltype
        real(dp),           intent(in)  :: theta,coords(:,:)
        integer,            intent(in)  :: edg(:,:),nedge
        real(dp),           intent(out) :: ctip(2,2)
    
    
    ! variables used for calculating clength
        real(kind=dp)   :: xo,yo,xp1,yp1,xp2,yp2,x1,y1,x2,y2,xct,yct,xct1,yct1,xct2,yct2
        real(kind=dp)   :: a1,b1,a2,b2,xmid,ymid,xmid1,ymid1,xmid2,ymid2,detlc
        integer :: nfailedge,iscross,ifedg(nedge)
        
        integer :: i,j
        
        
        !-----------------------------------------------------------!
        !           calculate approximate clength
        !-----------------------------------------------------------!
        ! initialize output variables
        ctip=zero
        ! initialize local variables
        xo=zero; yo=zero; xp1=zero; yp1=zero; xp2=zero; yp2=zero
        x1=zero; y1=zero; x2=zero; y2=zero; xct=zero; yct=zero
        xct1=zero; yct1=zero; xct2=zero; yct2=zero
        a1=zero; b1=zero; a2=zero; b2=zero
        xmid=zero; ymid=zero; xmid1=zero; ymid1=zero; xmid2=zero; ymid2=zero
        detlc=zero
        nfailedge=0; iscross=0; ifedg=0
        i=0; j=0
        
        ! find centroid
        select case(eltype)
            case('brick')
                xo=quarter*(coords(1,1)+coords(1,2)+coords(1,3)+coords(1,4))
                yo=quarter*(coords(2,1)+coords(2,2)+coords(2,3)+coords(2,4))
            case('wedge')
                xo=one_third*(coords(1,1)+coords(1,2)+coords(1,3))
                yo=one_third*(coords(2,1)+coords(2,2)+coords(2,3))
            case default
                write(msg_file,*)'unsupported elem type for crack partition from origin'
                call exit_function
        end select
        
        
        xp1=xo
        yp1=yo
        ! find xp2, yp2
        xp2=xp1+cos(theta/halfcirc*pi)
        yp2=yp1+sin(theta/halfcirc*pi)
        do i=1,nedge 
            ! tip corods of edge i
            x1=coords(1,edg(1,i))
            y1=coords(2,edg(1,i))
            x2=coords(1,edg(2,i))
            y2=coords(2,edg(2,i))
            iscross=0
            xct=zero
            yct=zero
            call klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct)
            if (iscross.gt.0) then
                nfailedge=nfailedge+1
                ifedg(nfailedge)=i
                ctip(:,nfailedge)=[xct,yct]
             endif
             if(nfailedge==2) exit ! found 2 broken edges already
        end do


        ! badly shaped element may have large angles; e.g.: 3 or all edges almost parallel to crack, then numerical error may
        ! prevent the algorithm from finding any broken edge or only one broken edge

        if(nfailedge==0) then
        ! use trial lines: connecting midpoints of two edges and find the most parrallel-to-crack one

            ! the following algorithm will find two broken edges
            nfailedge=2
            ! line equation constants of the crack
            a2=sin(theta/halfcirc*pi)
            b2=-cos(theta/halfcirc*pi)
            ! connecting midpoints of two edges and find the most parrallel one
            do i=1,nedge-1
                ! tip coords of edge i
                x1=coords(1,edg(1,i))
                y1=coords(2,edg(1,i))
                x2=coords(1,edg(2,i))
                y2=coords(2,edg(2,i))
                ! mid point of edge i
                xmid1=half*(x1+x2)
                ymid1=half*(y1+y2)
                ! loop over midpoints of other edges
                do j=i+1,nedge
                    ! tip coords of edge j
                    x1=coords(1,edg(1,j))
                    y1=coords(2,edg(1,j))
                    x2=coords(1,edg(2,j))
                    y2=coords(2,edg(2,j))
                    ! mid point of edge j
                    xmid2=half*(x1+x2)
                    ymid2=half*(y1+y2)
                    ! line equation constants of midpoint1-midpoint2
                    a1=ymid1-ymid2
                    b1=xmid2-xmid1
                    ! initialize detlc and intersection info
                    if(i==1 .and. j==2) then
                        detlc=a1*b2-a2*b1
                        xct1=xmid1
                        yct1=ymid1
                        xct2=xmid2
                        yct2=ymid2
                    end if
                    ! find the most parallel trial line and update stored info
                    if(abs(a1*b2-a2*b1)<abs(detlc)) then
                        xct1=xmid1
                        yct1=ymid1
                        xct2=xmid2
                        yct2=ymid2
                    endif
                end do
            end do
            ctip(:,1)=[xct1,yct1] 
            ctip(:,2)=[xct2,yct2]

        else if(nfailedge==1) then                
        ! use trial lines: connecting the existing crack tip (xp1,yp1) to the midpoints of the other 3 edges

            ! the following algorithm will find the second broken edge
            nfailedge=2
            ! existing crack tip coords
            xp1=xct
            yp1=yct
            ! crack line equation constants
            a2=sin(theta/halfcirc*pi)
            b2=-cos(theta/halfcirc*pi)
            ! find the midpoint which forms the most parrallel-to-crack line with (xp1,yp1)
            do i=1,nedge
                if (i==ifedg(1)) cycle
                ! tip coords of edge i
                x1=coords(1,edg(1,i))
                y1=coords(2,edg(1,i))
                x2=coords(1,edg(2,i))
                y2=coords(2,edg(2,i))
                ! mid point of edge i
                xmid=half*(x1+x2)
                ymid=half*(y1+y2)
                ! line equation constants of midpoint-(xp1,yp1)
                a1=ymid-yp1
                b1=xp1-xmid
                ! initialize detlc and intersection info
                if(i==1 .or. (ifedg(1)==1 .and. i==2)) then
                    detlc=a1*b2-a2*b1
                    xct=xmid
                    yct=ymid
                end if
                ! find the most parallel trial line and update stored info
                if(abs(a1*b2-a2*b1)<abs(detlc)) then
                    xct=xmid
                    yct=ymid
                endif
            end do
            ctip(:,2)=[xct,yct]
        end if
        

    end subroutine elem_ctips_origin





!********************************************************************************************
!           subroutine to assemble elem/subelem K and F to system/xelem K and F
!******************************************************************************************** 
    subroutine assembleKF(Kmat,Fvec,Ki,Fi,cnc)
    
        real(dp),intent(inout) :: Kmat(:,:),Fvec(:)
        real(dp),intent(in)    :: Ki(:,:), Fi(:)
        integer, intent(in)    :: cnc(:)
        
        ! local variable
        integer :: i, j, nsize
        i=0; j=0; nsize=0
        
        nsize=size(cnc)
        
        do i=1, nsize
            do j=1, nsize
                Kmat(cnc(j),cnc(i))=Kmat(cnc(j),cnc(i))+Ki(j,i)
            end do
            Fvec(cnc(i))=Fvec(cnc(i))+Fi(i)
        end do    
    
    end subroutine assembleKF



    
    
    end module toolkit_module
