module global_toolkit_module

implicit none

contains



  pure function glb_dee(dee, theta)
  ! Purpose :
  ! to transform d matrix from local to global coordinates 
  ! (for composite lamina, 2D or 3D, with fibre angle = theta)  
  
    use parameter_module, only : DP, ZERO, PI, HALFCIRC, TWO, NST=>NST_STANDARD

    real(DP), intent(in) :: dee(NST,NST), theta
    real(DP)             :: glb_dee(NST,NST)

    ! local variables
    real(DP) :: c, s
    
    glb_dee = ZERO
    c = ZERO
    s = ZERO

    c = cos(PI*theta/HALFCIRC)
    s = sin(PI*theta/HALFCIRC)
    
    ! d matrix in global coords stored first in local array glb_dee
    if (NST == 3) then
      glb_dee(1,1) = c*c*c*c*dee(1,1) + TWO*c*c*s*s*(dee(1,2)       &
   &               + TWO*dee(3,3)) + s*s*s*s*dee(2,2)
      glb_dee(1,2) = s*s*c*c*(dee(1,1) + dee(2,2) - four*dee(3,3))  &
   &               + (s*s*s*s+c*c*c*c)*dee(1,2)
      glb_dee(2,1) = glb_dee(1,2)
      glb_dee(2,2) = s*s*s*s*dee(1,1) + TWO*c*c*s*s*(dee(1,2)       &
   &               + TWO*dee(3,3)) + c*c*c*c*dee(2,2)
      glb_dee(1,3) = s*c*(c*c*(dee(1,1) - dee(1,2) - TWO*dee(3,3))  &
   &               + s*s*(dee(1,2) - dee(2,2) + TWO*dee(3,3)))
      glb_dee(3,1) = glb_dee(1,3)
      glb_dee(2,3) = s*c*(s*s*(dee(1,1) - dee(1,2) - TWO*dee(3,3))  &
   &               + c*c*(dee(1,2) - dee(2,2) + TWO*dee(3,3)))
      glb_dee(3,2) = glb_dee(2,3)
      glb_dee(3,3) = c*c*s*s*(dee(1,1)+dee(2,2)-2*dee(1,2))         &
   &			         + (c*c-s*s)**2*dee(3,3)
   
    else if (NST == 6) then
      glb_dee(1,1) = c*c*c*c*dee(1,1) + TWO*c*c*s*s*(dee(1,2)       &
   &               + TWO*dee(4,4)) + s*s*s*s*dee(2,2)
      glb_dee(1,2) = s*s*c*c*(dee(1,1) + dee(2,2) - four*dee(4,4))  &
   &               + (s*s*s*s+c*c*c*c)*dee(1,2)
      glb_dee(2,1) = glb_dee(1,2)
      glb_dee(2,2) = s*s*s*s*dee(1,1) + TWO*c*c*s*s*(dee(1,2)       &
   &               + TWO*dee(4,4)) + c*c*c*c*dee(2,2)
      glb_dee(1,3) = c*c*dee(1,3) + s*s*dee(2,3)
      glb_dee(3,1) = glb_dee(1,3)
      glb_dee(1,4) = s*c*(c*c*(dee(1,1) - dee(1,2) - TWO*dee(4,4))  &
   &               + s*s*(dee(1,2) - dee(2,2) + TWO*dee(4,4)))
      glb_dee(4,1) = glb_dee(1,4)
      glb_dee(2,3) = s*s*dee(1,3) + c*c*dee(2,3)
      glb_dee(3,2) = glb_dee(2,3)
      glb_dee(2,4) = s*c*(s*s*(dee(1,1) - dee(1,2) - TWO*dee(4,4))  &
   &               + c*c*(dee(1,2) - dee(2,2) + TWO*dee(4,4)))
      glb_dee(4,2) = glb_dee(2,4)
      glb_dee(3,3) = dee(3,3)
      glb_dee(3,4) = c*s*(dee(1,3) - dee(2,3))
      glb_dee(4,3) = glb_dee(3,4)
      glb_dee(5,5) = c*c*dee(5,5) + s*s*dee(6,6)
      glb_dee(5,6) = c*s*(dee(6,6) - dee(5,5))
      glb_dee(6,5) = glb_dee(5,6)
      glb_dee(6,6) = s*s*dee(5,5) + c*c*dee(6,6)
      glb_dee(4,4) = s*s*c*c*(dee(1,1) - TWO*dee(1,2) + dee(2,2))   &
   &               + (s*s - c*c)*(s*s - c*c)*dee(4,4)
   
    end if

  end function glb_dee




  pure function determinant_jacob(jacob)
  ! Purpose:
  ! returns the determinant of a jacobian matrix, 2D or 3D
  
    use parameter_module, only : DP, ZERO, NDIM
  
    real(DP), intent(in) ::   jacob(NDIM,NDIM)
    real(DP)             ::   determinant_jacob
    
    determinant_jacob = ZERO

    if (NDIM == 2) then
      determinant_jacob = jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1)
    else if (NDIM == 3) then
      determinant_jacob = jacob(1,1) * &
                  & (jacob(2,2)*jacob(3,3)-jacob(3,2)*jacob(2,3))
      determinant_jacob = determinant_jacob - jacob(1,2) * &
                  & (jacob(2,1)*jacob(3,3)-jacob(3,1)*jacob(2,3))
      determinant_jacob = determinant_jacob + jacob(1,3) * &
                  & (jacob(2,1)*jacob(3,2)-jacob(3,1)*jacob(2,2))
    end if

  end function determinant_jacob




  pure subroutine invert_jacob(jacob, istat, emsg, detj)
  ! Purpose:
  ! calculate the inverse of a jacobian matrix and update onto itself
  
    use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                          & ZERO, NDIM, SMALLNUM

    real(DP),                 intent(inout) :: jacob(NDIM,NDIM)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg
    real(DP),       optional, intent(in)    :: detj  
    
    ! local variables:
    ! local copy of intent inout arg
    real(DP) :: jacob_lcl(NDIM,NDIM)
    ! inverse of jacob matrix
    real(DP) :: inv_jacob(NDIM,NDIM)
    ! local copy of optional arg.
    real(DP) :: det
    
    ! initialize intent out and local variables
    istat = STAT_SUCCESS
    emsg  = ''
    jacob_lcl = ZERO
    inv_jacob = ZERO
    det     = ZERO
    
    ! assign values to local variables
    jacob_lcl = jacob
    
    ! get the determinant of jacob matrix; pass local copy of jacob for safety
    if(present(detj)) then
      det = detj
    else
      det = determinant_jacob(jacob_lcl)
    end if
    
    ! check to see if the matrix is singular
    if(det < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'the jacobian matrix is singular, invert_jacob, &
            & global_toolkit_module'
      return
    end if

    ! calculate the inverse of jacob matrix
    if(NDIM == 2) then
      inv_jacob(1,1)=jacob(2,2)
      inv_jacob(2,1)=-jacob(2,1)
      inv_jacob(1,2)=-jacob(1,2)
      inv_jacob(2,2)=jacob(1,1)
    else if(NDIM == 3) then
      inv_jacob(1,1)=jacob(2,2)*jacob(3,3)-jacob(3,2)*jacob(2,3)
      inv_jacob(2,1)=jacob(3,1)*jacob(2,3)-jacob(2,1)*jacob(3,3)
      inv_jacob(3,1)=jacob(2,1)*jacob(3,2)-jacob(3,1)*jacob(2,2)
      inv_jacob(1,2)=jacob(3,2)*jacob(1,3)-jacob(1,2)*jacob(3,3)
      inv_jacob(2,2)=jacob(1,1)*jacob(3,3)-jacob(3,1)*jacob(1,3)
      inv_jacob(3,2)=jacob(3,1)*jacob(1,2)-jacob(1,1)*jacob(3,2)
      inv_jacob(1,3)=jacob(1,2)*jacob(2,3)-jacob(2,2)*jacob(1,3)
      inv_jacob(2,3)=jacob(2,1)*jacob(1,3)-jacob(1,1)*jacob(2,3)
      inv_jacob(3,3)=jacob(1,1)*jacob(2,2)-jacob(2,1)*jacob(1,2)
    end if
    inv_jacob = inv_jacob/det
    
    ! update to jacob its inverse before successful return
    jacob = inv_jacob

  end subroutine invert_jacob




  pure function beemat(gn, nnode)
  ! Purpose:
  ! inputs the gradient of shape functions and
  ! returns the strain-displacement matrix for infinitesimal deformation
  ! order convention of strain terms in the strain vector: 
  ! eps1, eps2, eps3, gamma12, gamma13, gamma23
  
    use parameter_module, only : DP, ZERO, NDIM, NST=>NST_STANDARD

    integer,  intent(in) :: nnode
    real(DP), intent(in) :: gn(nnode, NDIM)
    real(DP)             :: beemat(NST, nnode*NDIM)

    ! local variables  
    real(DP) :: x, y, z
    integer  :: k, l, m, n
    
    beemat = ZERO
    x = ZERO
    y = ZERO
    z = ZERO
    k = 0; l = 0; m = 0; n = 0

    if (NST == 3) then
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
      
    else if (NST == 6) then
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
      
    end if

  end function beemat




  pure function lcl_strain(strain, theta)
  ! Purpose:
  ! transfer strains from global to local coordinate systems
  ! (for composite lamina, ply angle = theta)
  
    use parameter_module, only : DP, ZERO, PI, HALFCIRC, ONE, TWO, &
                          & NST=>NST_STANDARD

    real(DP),intent(in)  :: strain(NST), theta
    real(DP)             :: lcl_strain(NST)

    real(DP) :: c, s, Q_matrix(NST,NST)
    
    lcl_strain = ZERO
    c = ZERO
    s = ZERO
    Q_matrix = ZERO

    c = cos(PI*theta/HALFCIRC)
    s = sin(PI*theta/HALFCIRC)

    ! calculate rotational matrix Q
    if (NST == 3) then
      Q_matrix(1,1)=c*c
      Q_matrix(1,2)=s*s
      Q_matrix(1,3)=c*s
      Q_matrix(2,1)=s*s
      Q_matrix(2,2)=c*c
      Q_matrix(2,3)=-c*s
      Q_matrix(3,1)=-TWO*c*s
      Q_matrix(3,2)=TWO*c*s
      Q_matrix(3,3)=c*c-s*s
      
    else if (NST == 6) then
      Q_matrix(1,1)=c*c
      Q_matrix(1,2)=s*s
      Q_matrix(1,4)=c*s
      Q_matrix(2,1)=s*s
      Q_matrix(2,2)=c*c
      Q_matrix(2,4)=-c*s
      Q_matrix(3,3)=ONE
      Q_matrix(4,1)=-TWO*c*s
      Q_matrix(4,2)=TWO*c*s
      Q_matrix(4,4)=c*c-s*s
      Q_matrix(5,5)=c
      Q_matrix(5,6)=-s
      Q_matrix(6,5)=s
      Q_matrix(6,6)=c
      
    end if
    
    ! rotate global strain to obtain local strain
    lcl_strain = matmul(Q_matrix,strain)

  end function lcl_strain




  pure subroutine normalize(a, is_zero_vect, amag)
  ! Purpose:
  ! normalize vector and return its magnitude
  
  use parameter_module, only : DP, ZERO, SMALLNUM
     
    real(DP), intent(inout) :: a(:)
    logical,  intent(out)   :: is_zero_vect
    real(DP), optional, intent(out) :: amag
   
    ! local copy of amag
    real(DP) :: mag
    
    ! initialize intent out and local variables
    is_zero_vect = .false.
    if(present(amag)) amag = ZERO
    mag = ZERO
   
    ! calculate the length of a
    mag = sqrt(dot_product(a,a))
    
    if (mag > SMALLNUM) then
    ! a is NOT a zero-length vector;
    ! normalize a, update amag (if present) and return
      is_zero_vect = .false.
      a = a / mag
      if(present(amag)) amag = mag
      return
    else
    ! a is a zero-length vector;
    ! do NOT update a, just return
      is_zero_vect = .true.
      return
    end if
  
  end subroutine normalize




  pure function cross_product3d(a, b)
  ! Purpose:
  ! function to compute the cross-product of TWO 3D vectors, a * b
  
  use parameter_module, only : DP
  
    real(DP), intent(in)  :: a(3), b(3)
    real(DP)              :: cross_product3d(3)

    cross_product3d(1) = a(2) * b(3) - a(3) * b(2)
    cross_product3d(2) = a(3) * b(1) - a(1) * b(3)
    cross_product3d(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product3d




  pure function distance(a, b, n)
  ! Purpose:
  ! function to compute the distance between TWO vectors
  
  use parameter_module, only : DP, ZERO
  
    integer,  intent(in)  :: n
    real(DP), intent(in)  :: a(n), b(n)
    real(DP)              :: distance
    
    real(DP) :: c(n)
    
    c = ZERO
    
    c = a - b
    
    distance=sqrt(dot_product(c,c))
    
  end function distance




  
  pure subroutine edge_crack_cross2d (edge, crack, cross_stat, cross_point)
  ! Purpose:
  ! to checks if TWO line sections cross each other in 2D, 
  ! return the cross status and locate the crossing point if they cross
  
  use parameter_module, only : DP, ZERO, 
                        & CROSS_ON_EDGE_ON_CRACK,  CROSS_ON_EDGE_OFF_CRACK,  &
                        & CROSS_OFF_EDGE_ON_CRACK, CROSS_OFF_EDGE_OFF_CRACK, &
                        & EDGE_CRACK_PARALLEL, ENDNODE_CLEARANCE_RATIO
  
    ! dummy args list:
    ! edge        : nodal coords of two end nodes of edge  line segment
    ! crack       : nodal coords of two end nodes of crack line segment
    ! cross_stat  : status of crossing
    ! cross_point : coords of the cross point
    real(DP), intent(in)  :: edge(2,2)
    real(DP), intent(in)  :: crack(2,2)
    integer,  intent(out) :: cross_stat
    real(DP), intent(out) :: cross_point(2)
    
    ! local variables
    ! a/b/cLj,      j=1,2 : line equation constants a, b, c of Line j
    ! det                 : determinant of the system equations of two lines
    real(DP) :: aL1, bL1, cL1
    real(DP) :: aL2, bL2, cL2
    real(DP) :: det
    
    ! initialize intent out and local variables
    ! intent out
    cross_stat  = 0
    cross_point = ZERO
    ! line equation constants and determinant
    aL1 = ZERO; bL1 = ZERO; cL1 = ZERO
    aL2 = ZERO; bL2 = ZERO; cL2 = ZERO
    det = ZERO
    
    !
    !**** algorithm for finding the intersection of two lines ****
    ! equation of a line: aL * x + bL * y = cL
    ! to find the cross point of line 1 (edge) and line 2 (crack) 
    ! is to solve the system equations of two lines:
    ! aL1 * x + bL1 * y = cL1 (equation of line 1)
    ! aL2 * x + bL2 * y = cL2 (equation of line 2)
    ! calculate det = aL1 * bL2 - aL2 * bL1 
    ! if det == ZERO, lines are parallel, no intersection possible;
    ! if det /= ZERO, lines are NOT parallel, the cross point coordinates are:
    ! x_crosspoint = (bL2 * cL1 - bL1 * cL2) / det
    ! y_crosspoint = (aL1 * cL2 - aL2 * cL1) / det
    !  
    
    ! line 1 equation (edge)
    call line_equation(edge, aL1, bL1, cL1)
    
    ! line 2 equation (crack)
    call line_equation(crack, aL2, bL2, cL2)
    
    ! calculate determinant
    det   = aL1 * bL2 - aL2 * bL1
   
   
    !**** MAIN CALCULATIONS ****
    
    ! select what to do based on the value of det:
    ! if det is NONZERO: intersection is possible
    !   - calculate cross point coordinates
    !   - if it is on edge:
    !       update the cross_stat:
    !         * if it is on crack : cross_stat = CROSS_ON_EDGE_ON_CRACK
    !         * if it is off crack: cross_stat = CROSS_ON_EDGE_OFF_CRACK
    !       adjust cross point coordinates away from edge end nodes,
    !       then exit the program
    !   - if it is off edge:
    !       update the cross_stat:
    !         * if it is on crack : cross_stat = CROSS_OFF_EDGE_ON_CRACK
    !         * if it is off crack: cross_stat = CROSS_OFF_EDGE_OFF_CRACK
    !       ZERO cross point coordinates,
    !       then exit the program
    !
    ! if det is ZERO: edge and crack are parallel, intersection NOT possible
    !   update the cross_stat: 
    !     * cross_stat = EDGE_CRACK_PARALLEL
    !   ZERO cross point coordinates,
    !   then exit the program
    
    det_if: if (abs(det) > SMALLNUM) then 
    ! intersection possible
    
        ! calculate cross point coordinates
        cross_point(1) = (bL2 * cL1 - bL1 * cL2) / det 
        cross_point(2) = (aL1 * cL2 - aL2 * cL1) / det

        ! update cross status
        cross_edge_if: if ( in_range(cross_point, edge) ) then
        ! the cross point is on this edge
            if ( in_range(cross_point, crack) ) then
            ! cross point on edge and also on crack
              cross_stat = CROSS_ON_EDGE_ON_CRACK   ! = 2   
            else
            ! cross point on edge but not on crack
              cross_stat = CROSS_ON_EDGE_OFF_CRACK  ! = 1   
            end if
            ! adjust cross point positions away from edge end nodes
            call adjust_cross_point(cross_point, edge)
        else cross_edge_if 
        ! the cross point is out of the range of the edge
            if ( in_range(cross_point, crack) ) then
            ! cross point off edge but on crack
              cross_stat = CROSS_OFF_EDGE_ON_CRACK  ! = -1   
            else
            ! cross point off edge and off crack
              cross_stat = CROSS_OFF_EDGE_OFF_CRACK ! = -2   
            end if
            ! zero cross point coordinates
            cross_point = ZERO
        end if cross_edge_if
        
        ! exit program
        return
        
    else det_if
    ! lines are parallel; no intersection possible
    
        ! just update output variables and return
        cross_stat  = EDGE_CRACK_PARALLEL           ! = 0
        cross_point = ZERO
        ! exit program
        return
        
    end if det_if
    
    !**** END MAIN CALCULATIONS ****
    
    
    contains
    
    ! internal procedures
    
    
    pure subroutine line_equation(line_seg, a, b, c)
    ! Purpose:
    ! obtain line equation constants a, b, c with passed-in line segment
    
      ! dummy args
      real(DP), intent(in)  :: line_seg(2,2)
      real(DP), intent(out) :: a, b, c
      
      ! local variables:
      real(DP) :: xN1L,  yN1L,  xN2L,  yN2L

      ! end node coords of line seg
      xN1L = line_seg(1,1)
      yN1L = line_seg(2,1)
      xN2L = line_seg(1,2)
      yN2L = line_seg(2,2)
      a   = yN2L - yN1L
      b   = xN1L - xN2L
      c   = a * xN1L + b * yN1L 
    
    end subroutine line_equation
    
    
    pure function in_range(point, line_seg)
    ! Purpose:
    ! return true of the point lies within the range of the line segment
      
      ! dummy args:
      real(DP), intent(in) :: point(2), line_seg(2,2)
      logical              :: in_range
      
      ! local variables:
      real(DP) :: xp, yp
      real(DP) :: xN1L,  yN1L,  xN2L,  yN2L
      real(DP) :: xminL, yminL, xmaxL, ymaxL
      
      ! coords of point
      xp   = point(1)
      yp   = point(2)
      ! end node coords of line seg
      xN1L = line_seg(1,1)
      yN1L = line_seg(2,1)
      xN2L = line_seg(1,2)
      yN2L = line_seg(2,2)
      ! range of line_seg
      xminL = min(xN1L, xN2L)
      yminL = min(yN1L, yN2L)
      xmaxL = max(xN1L, xN2L)
      ymaxL = max(yN1L, yN2L)
      
      ! point is in range of line seg if
      ! xminL <= xp <= xmaxL and
      ! yminL <= yp <= ymaxL
      if (xminL-SMALLNUM < xp .and. xp < xmaxL+SMALLNUM .and. &
      &   yminL-SMALLNUM < yp .and. yp < ymaxL+SMALLNUM) then
        in_range = .true.
      else
        in_range = .false.
      end if
    
    end function in_range
    
    
    pure subroutine adjust_cross_point(point, line_seg)
    ! Purpose:
    ! adjust the coords of cross_point to be away from the end nodes of edge
    
      real(DP), intent(inout) :: point(2)
      real(DP), intent(in)    :: line_seg(2,2)
      
      ! local variables:
      real(DP) :: xp, yp
      real(DP) :: xN1L,  yN1L,  xN2L,  yN2L
      real(DP) :: xminL, yminL, xmaxL, ymaxL
      
      ! coords of point
      xp   = point(1)
      yp   = point(2)
      ! end node coords of line seg
      xN1L = line_seg(1,1)
      yN1L = line_seg(2,1)
      xN2L = line_seg(1,2)
      yN2L = line_seg(2,2)
      ! range of line_seg
      xminL = min(xN1L, xN2L)
      yminL = min(yN1L, yN2L)
      xmaxL = max(xN1L, xN2L)
      ymaxL = max(yN1L, yN2L)
      
      ! adjust point positions away from line seg end nodes
      xp = max(xp, xminL + ENDNODE_CLEARANCE_RATIO * (xmaxL-xminL))
      xp = min(xp, xmaxL - ENDNODE_CLEARANCE_RATIO * (xmaxL-xminL))
      yp = max(yp, yminL + ENDNODE_CLEARANCE_RATIO * (ymaxL-yminL))
      yp = min(yp, ymaxL - ENDNODE_CLEARANCE_RATIO * (ymaxL-yminL))
      
      ! update point
      point = [xp, yp]
    
    end subroutine adjust_cross_point
    
         
  end subroutine edge_crack_cross2d
  
  
  
  




!********************************************************************************************
!           subroutine to find TWO crack tips on the edges of an element, passing origin
!******************************************************************************************** 
subroutine elem_ctips_origin(eltype,theta,coords,edg,nedge,ctip)
    character(len=*),   intent(in)  :: eltype
    real(dp),           intent(in)  :: theta,coords(:,:)
    integer,            intent(in)  :: edg(:,:),nedge
    real(dp),           intent(out) :: ctip(2,2)


! variables used for calculating clength
    real(DP)   :: xo,yo,xp1,yp1,xp2,yp2,x1,y1,x2,y2,xct,yct,xct1,yct1,xct2,yct2
    real(DP)   :: a1,b1,a2,b2,xmid,ymid,xmid1,ymid1,xmid2,ymid2,detlc
    integer :: nfailedge,iscross,ifedg(nedge)
    
    integer :: i,j
    
    
    !-----------------------------------------------------------!
    !           calculate approximate clength
    !-----------------------------------------------------------!
    ! initialize output variables
    ctip=ZERO
    ! initialize local variables
    xo=ZERO; yo=ZERO; xp1=ZERO; yp1=ZERO; xp2=ZERO; yp2=ZERO
    x1=ZERO; y1=ZERO; x2=ZERO; y2=ZERO; xct=ZERO; yct=ZERO
    xct1=ZERO; yct1=ZERO; xct2=ZERO; yct2=ZERO
    a1=ZERO; b1=ZERO; a2=ZERO; b2=ZERO
    xmid=ZERO; ymid=ZERO; xmid1=ZERO; ymid1=ZERO; xmid2=ZERO; ymid2=ZERO
    detlc=ZERO
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
        xct=ZERO
        yct=ZERO
        call klinecross(x1,y1,x2,y2,xp1,yp1,xp2,yp2,iscross,xct,yct)
        if (iscross.gt.0) then
            nfailedge=nfailedge+1
            ifedg(nfailedge)=i
            ctip(:,nfailedge)=[xct,yct]
         endif
         if(nfailedge==2) exit ! found 2 broken edges already
    end do


    ! badly shaped element may have large angles; e.g.: 3 or all edges almost parallel to crack, then numerical error may
    ! prevent the algorithm from finding any broken edge or only ONE broken edge

    if(nfailedge==0) then
    ! use trial lines: connecting midpoints of TWO edges and find the most parrallel-to-crack ONE

        ! the following algorithm will find TWO broken edges
        nfailedge=2
        ! line equation constants of the crack
        a2=sin(theta/halfcirc*pi)
        b2=-cos(theta/halfcirc*pi)
        ! connecting midpoints of TWO edges and find the most parrallel ONE
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





end module global_toolkit_module
