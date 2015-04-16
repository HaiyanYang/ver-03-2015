module global_toolkit_module
!
!  Purpose:
!   this module contains a collection of useful procedures used in several other
!   modules.
!
!   **** PROGRAMMING PRINCIPLES: ****
!   all procedures are PURE.
!   FIX the sizes of dummy args whenever possible.
!     - suffix '2d'/'3d' is added at the end of procedure name to indicate
!       its supported dimension
!     - the sizes can be set by private parameters in the procedure.
!   AVOID support of variable sizes of dummy args
!     - Keep It Simple and Stupid (KISS)
!     - forbid unexpected input arg. sizes, avoid explicit checking
!   AVOID optional arg. KISS
!   if error status and message are needed, make them COMPULSARY INPUTS
!     - do NOT give the option to not check for error status
!   ALWAYS CHECK dummy arg sizes if variable sizing has to be allowed,
!       and flag error when their sizes are not expected
!   ALWAYS CHECK the validity of intent in/inout dummy args if they are expected
!       to be within a certain range of values
!       EXCEPTION: if the input arg is a size parameter of other input args, 
!       then its checking can be left to the compilor
!   ALWAYS USE LOCAL COPY OF INTENT(INOUT) and OPTIONAL dummy args. for 
!       calculation, only update to the original arg. before successful return
!   
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    15/04/15  B. Y. Chen            Original code
!

implicit none

contains



  pure function glb_dee3d (dee, theta)
  ! Purpose :
  ! to transform d matrix from local to global coordinates
  ! (for composite lamina, 3D, with fibre angle = theta)

    use parameter_module, only : DP, ZERO, PI, HALFCIRC, TWO, FOUR

    ! define private parameters
    ! NST: no. of stress/strain terms, used to set the size of dummy args
    integer, parameter   :: NST = 6

    real(DP), intent(in) :: dee(NST,NST), theta
    real(DP)             :: glb_dee3d(NST,NST)

    ! local variables
    real(DP) :: c, s

    glb_dee3d = ZERO
    c = ZERO
    s = ZERO

    c = cos(PI*theta/HALFCIRC)
    s = sin(PI*theta/HALFCIRC)

    ! d matrix in global coords stored first in local array glb_dee3d
    glb_dee3d(1,1) = c*c*c*c*dee(1,1) + TWO * c*c*s*s*(dee(1,2)       &
                 & + TWO * dee(4,4)) + s*s*s*s*dee(2,2)
    glb_dee3d(1,2) = s*s*c*c*(dee(1,1) + dee(2,2) - FOUR * dee(4,4))  &
                 & + (s*s*s*s+c*c*c*c)*dee(1,2)
    glb_dee3d(2,1) = glb_dee3d(1,2)
    glb_dee3d(2,2) = s*s*s*s*dee(1,1) + TWO * c*c*s*s*(dee(1,2)       &
                 & + TWO * dee(4,4)) + c*c*c*c*dee(2,2)
    glb_dee3d(1,3) = c*c*dee(1,3) + s*s*dee(2,3)
    glb_dee3d(3,1) = glb_dee3d(1,3)
    glb_dee3d(1,4) = s*c*(c*c*(dee(1,1) - dee(1,2) - TWO * dee(4,4))  &
                 & + s*s*(dee(1,2) - dee(2,2) + TWO * dee(4,4)))
    glb_dee3d(4,1) = glb_dee3d(1,4)
    glb_dee3d(2,3) = s*s*dee(1,3) + c*c*dee(2,3)
    glb_dee3d(3,2) = glb_dee3d(2,3)
    glb_dee3d(2,4) = s*c*(s*s*(dee(1,1) - dee(1,2) - TWO * dee(4,4))  &
                 & + c*c*(dee(1,2) - dee(2,2) + TWO * dee(4,4)))
    glb_dee3d(4,2) = glb_dee3d(2,4)
    glb_dee3d(3,3) = dee(3,3)
    glb_dee3d(3,4) = c*s*(dee(1,3) - dee(2,3))
    glb_dee3d(4,3) = glb_dee3d(3,4)
    glb_dee3d(5,5) = c*c*dee(5,5) + s*s*dee(6,6)
    glb_dee3d(5,6) = c*s*(dee(6,6) - dee(5,5))
    glb_dee3d(6,5) = glb_dee3d(5,6)
    glb_dee3d(6,6) = s*s*dee(5,5) + c*c*dee(6,6)
    glb_dee3d(4,4) = s*s*c*c*(dee(1,1) - TWO * dee(1,2) + dee(2,2))   &
                 & + (s*s - c*c)*(s*s - c*c)*dee(4,4)

  end function glb_dee3d




  pure function determinant3d (jacob)
  ! Purpose:
  ! returns the determinant of a 3D jacobian matrix

    use parameter_module, only : DP, ZERO

    ! define private parameters
    ! NDIM: no. of dimensions, used to set the size of dummy arg
    integer, parameter   :: NDIM = 3

    real(DP), intent(in) ::   jacob(NDIM,NDIM)
    real(DP)             ::   determinant3d

    determinant3d = ZERO

    determinant3d = jacob(1,1) * &
                & (jacob(2,2)*jacob(3,3)-jacob(3,2)*jacob(2,3))
    determinant3d = determinant3d - jacob(1,2) * &
                & (jacob(2,1)*jacob(3,3)-jacob(3,1)*jacob(2,3))
    determinant3d = determinant3d + jacob(1,3) * &
                & (jacob(2,1)*jacob(3,2)-jacob(3,1)*jacob(2,2))

  end function determinant3d




  pure subroutine invert_self3d (jacob, istat, emsg, detj)
  ! Purpose:
  ! calculate the inverse of a 3D jacobian matrix and update onto itself
  ! the invertability of jacob must be checked before inverting it

    use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                          & ZERO, SMALLNUM

    ! define private parameters
    ! NDIM: no. of dimensions, used to set the size of dummy arg
    integer, parameter   :: NDIM = 3

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
    istat     = STAT_SUCCESS
    emsg      = ''
    jacob_lcl = ZERO
    inv_jacob = ZERO
    det       = ZERO

    ! assign values to local variables
    jacob_lcl = jacob

    ! get the determinant of jacob matrix; pass local copy of jacob for safety
    if(present(detj)) then
      det = detj
    else
      det = determinant3d(jacob_lcl)
    end if

    ! check to see if the matrix is singular; if so, flag error and exit program
    if(det < SMALLNUM) then
      istat = STAT_FAILURE
      emsg  = 'the jacobian matrix is singular, invert_self3d, &
            & global_toolkit_module'
      return
    end if

    ! calculate the inverse of jacob matrix
    inv_jacob(1,1)=jacob(2,2)*jacob(3,3)-jacob(3,2)*jacob(2,3)
    inv_jacob(2,1)=jacob(3,1)*jacob(2,3)-jacob(2,1)*jacob(3,3)
    inv_jacob(3,1)=jacob(2,1)*jacob(3,2)-jacob(3,1)*jacob(2,2)
    inv_jacob(1,2)=jacob(3,2)*jacob(1,3)-jacob(1,2)*jacob(3,3)
    inv_jacob(2,2)=jacob(1,1)*jacob(3,3)-jacob(3,1)*jacob(1,3)
    inv_jacob(3,2)=jacob(3,1)*jacob(1,2)-jacob(1,1)*jacob(3,2)
    inv_jacob(1,3)=jacob(1,2)*jacob(2,3)-jacob(2,2)*jacob(1,3)
    inv_jacob(2,3)=jacob(2,1)*jacob(1,3)-jacob(1,1)*jacob(2,3)
    inv_jacob(3,3)=jacob(1,1)*jacob(2,2)-jacob(2,1)*jacob(1,2)

    inv_jacob = inv_jacob/det

    ! update intent(inout) args before successful return
    jacob = inv_jacob

  end subroutine invert_self3d




  pure function beemat3d (gn, nnode)
  ! Purpose:
  ! inputs the gradient of shape functions and
  ! returns the strain-displacement matrix for infinitesimal deformation
  ! order convention of strain terms in the strain vector:
  ! eps1, eps2, eps3, gamma12, gamma13, gamma23

    use parameter_module, only : DP, ZERO

    ! define private parameters, used to set the size of dummy arg
    ! NDIM: no. of dimensions
    ! NST : no. of stress/strain terms
    integer, parameter   :: NDIM = 3, NST = 6

    integer,  intent(in) :: nnode
    real(DP), intent(in) :: gn(nnode, NDIM)
    real(DP)             :: beemat3d(NST, nnode*NDIM)

    ! local variables
    real(DP) :: x, y, z
    integer  :: k, l, m, n

    ! initialize intent out and local variables
    beemat3d = ZERO
    x = ZERO
    y = ZERO
    z = ZERO
    k = 0; l = 0; m = 0; n = 0

    ! nnode is size parameter, let compilor to check it
    ! gn is fixed sized and can take any value, no check

    do m=1, nnode
      n = 3*m
      k = n-1
      l = k-1
      x = gn(m,1)
      y = gn(m,2)
      z = gn(m,3)
      beemat3d(1,l) = x
      beemat3d(4,k) = x
      beemat3d(5,n) = x
      beemat3d(2,k) = y
      beemat3d(4,l) = y
      beemat3d(6,n) = y
      beemat3d(3,n) = z
      beemat3d(5,l) = z
      beemat3d(6,k) = z
    end do

  end function beemat3d




  pure function lcl_strain3d (strain, theta)
  ! Purpose:
  ! transfer strains from global to local coordinate systems
  ! (for 3D composite lamina, ply angle = theta)

    use parameter_module, only : DP, ZERO, PI, HALFCIRC, ONE, TWO

    ! define private parameters, used to set the size of dummy arg
    ! NST : no. of stress/strain terms
    integer, parameter   :: NST = 6

    real(DP), intent(in) :: strain(NST), theta
    real(DP)             :: lcl_strain3d(NST)

    real(DP) :: c, s, Q_matrix(NST,NST)

    lcl_strain3d = ZERO
    c = ZERO
    s = ZERO
    Q_matrix = ZERO

    ! strain and theta can take any value, nothing to check

    c = cos(PI*theta/HALFCIRC)
    s = sin(PI*theta/HALFCIRC)

    ! calculate rotational matrix Q
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

    ! rotate global strain to obtain local strain
    lcl_strain3d = matmul(Q_matrix,strain)

  end function lcl_strain3d




  pure subroutine normalize_vect (a, is_zero_vect, amag)
  ! Purpose:
  ! normalize vector and return its magnitude

  use parameter_module, only : DP, ZERO, SMALLNUM

    real(DP),           intent(inout) :: a(:)
    logical,            intent(out)   :: is_zero_vect
    real(DP), optional, intent(out)   :: amag

    ! local copy of amag
    real(DP) :: mag

    ! initialize intent out and local variables
    is_zero_vect = .false.
    if(present(amag)) amag = ZERO
    mag = ZERO

    ! a can take any value, nothing to check

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

  end subroutine normalize_vect




  pure function cross_product3d (a, b)
  ! Purpose:
  ! function to compute the cross-product of TWO 3D vectors, a * b

  use parameter_module, only : DP

    real(DP), intent(in)  :: a(3), b(3)
    real(DP)              :: cross_product3d(3)

    cross_product3d(1) = a(2) * b(3) - a(3) * b(2)
    cross_product3d(2) = a(3) * b(1) - a(1) * b(3)
    cross_product3d(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product3d




  pure function distance (a, b, n)
  ! Purpose:
  ! function to compute the distance between TWO vectors
  ! their common size n is set as a required input to avoid having to
  ! explicitly check for their size consistency

  use parameter_module, only : DP, ZERO

    integer,  intent(in)  :: n
    real(DP), intent(in)  :: a(n), b(n)
    real(DP)              :: distance
    ! local variable
    real(DP) :: c(n)
    ! initialize local variable
    c = ZERO

    ! intent in variables a and b are fixed-sized and may take any value
    ! intent in variable n must > 0; it is a dummy arg size parameter,
    ! so the checking of n is left to the compilor

    c = a - b

    distance=sqrt(dot_product(c,c))

  end function distance




  pure subroutine edge_crack_cross2d (edge, crack, cross_stat, cross_point)
  ! Purpose:
  ! to checks if TWO line sections cross each other in 2D,
  ! return the cross status and locate the crossing point if they cross

  use parameter_module, only : DP, ZERO, SMALLNUM, &
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

    ! intent in dummy args' sizes are fixed, and they may take on any value,
    ! no need to check their validity

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




  pure subroutine crack_elem_centroid2d (nedge, crack_angle, coords, &
  & nodes_on_edges, istat, emsg, edge_crack_points, crack_edge_indices)
  ! Purpose:
  ! to find TWO cross points between the crack line (passing the centroid) and
  ! the edges of an element, and the two edges crossed.

  use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE, &
                        & ZERO, HALF, HALFCIRC, PI, SMALLNUM,             &
                        & CROSS_ON_EDGE_ON_CRACK,  CROSS_ON_EDGE_OFF_CRACK

    ! list of dummy args:
    ! - nedge              : no. of edges in this element; = no. of nodes
    ! - crack_angle
    ! - coords             : nodal coordinates of this element
    ! - nodes_on_edges(:,i) : indices of the two end nodes of the edge i
    ! - edge_crack_points  : coords of the two crack points on two edges
    ! - crack_edge_indices : indices of the two edges passed by the crack
    integer,  intent(in)  :: nedge
    real(DP), intent(in)  :: crack_angle
    real(DP), intent(in)  :: coords(2, nedge)
    integer,  intent(in)  :: nodes_on_edges(2, nedge)
    integer,                  intent(out) :: istat
    character(len=MSGLENGTH), intent(out) :: emsg
    real(DP),                 intent(out) :: edge_crack_points(2,2)
    integer,        optional, intent(out) :: crack_edge_indices(2)

    ! local variables
    ! - centroid          : coordinates of element centroid
    ! - theta             : ply angle in radiant
    ! - crack_unit_vector : unit vector along the crack line
    ! - coords_crack      : coordinates of end nodes of crack line segment
    ! - coords_edge       : coordinates of end nodes of an element edge
    ! - cross_point       : coordinates of the cross point btw crack and edge
    ! - n_crack_edge      : no. of cracked edges
    ! - cross_stat        : status of crossing btw a crack and an edge
    ! - crack_edge_IDs    : local copy of optional arg crack_edge_indices
    real(DP) :: centroid(2)
    real(DP) :: theta, crack_unit_vect(2)
    real(DP) :: coords_crack(2,2), coords_edge(2,2), cross_point(2)
    integer  :: n_crack_edge, cross_stat, crack_edge_IDs(2)

    integer :: i, j

    ! initialize intent out variables
    istat             = STAT_SUCCESS
    emsg              = ''
    edge_crack_points = ZERO
    if (present(crack_edge_indices)) crack_edge_indices = 0
    ! initialize local variables
    centroid     = ZERO
    theta        = ZERO
    crack_unit_vect = ZERO
    coords_crack = ZERO
    coords_edge  = ZERO
    cross_point  = ZERO
    cross_stat     = 0
    n_crack_edge   = 0
    crack_edge_IDs = 0
    i=0; j=0

    ! check validity of inputs: nedge, crack_angle, coords, nodes_on_edges
    ! nedge must be > 0, this check is omitted as if it is not satisfied,
    ! compilor should flag error as it sets the size of coords and nodes_on_edges
    ! crack_angle and coords can take any values
    ! nodes_on_edges must be >= 1, this need to be checked
    if ( any(nodes_on_edges < 1) ) then
      istat = STAT_FAILURE
      emsg  = 'end node indices must be >=1, crack_elem_centroid2d, &
      &global_toolkit_module'
      return
    end if

    ! find centroid
    centroid(1) = sum(coords(1,:))/real(nedge,DP)
    centroid(2) = sum(coords(2,:))/real(nedge,DP)

    ! find crack line unit vector
    theta           = crack_angle / HALFCIRC * PI
    crack_unit_vect = [cos(theta), sin(theta)]

    ! set centroid as node 1 of crack line segment
    coords_crack(:,1) = centroid(:)
    ! set node 2 coordinates of crack line segment
    ! to be a unit distance away from node 1 along the crack line
    coords_crack(:,2) = coords_crack(:,1) + crack_unit_vect(:)



    !**** MAIN CALCULATIONS ****

    ! loop over all edges to find TWO cross points with the crack line
    ! theoretically, this should find EXACTLY TWO cross points, as the
    ! crack line passes the element centroid and it must cross TWO edges
    do i = 1, nedge
      ! the two end nodes' coords of edge i
      coords_edge(:, 1) = coords(:, nodes_on_edges(1,i))
      coords_edge(:, 2) = coords(:, nodes_on_edges(2,i))
      ! zero cross_stat and cross_point for reuse
      cross_stat  = 0
      cross_point = ZERO
      ! check if the edge and crack cross each other
      call edge_crack_cross2d (coords_edge, coords_crack, cross_stat, cross_point)
      ! check cross_stat, update only if cross point is on edge
      if (cross_stat == CROSS_ON_EDGE_ON_CRACK .or. &
      &   cross_stat == CROSS_ON_EDGE_OFF_CRACK) then
          ! increase the no. of cracked edges
          n_crack_edge = n_crack_edge + 1
          ! store the index of this cracked edge
          crack_edge_IDs(n_crack_edge) = i
          ! store the cross point coords of this cracked edge
          edge_crack_points(:, n_crack_edge) = cross_point(:)
       endif
       ! exit the loop when TWO broken edges are already found
       if (n_crack_edge == 2) exit
    end do

    !**** END MAIN CALCULATIONS ****

    ! check if indeed two cracked edges have been found; if not, then the elem
    ! is likely to be poorly shaped, flag an error, clean up and return
    if (n_crack_edge /= 2) then
      istat = STAT_FAILURE
      emsg  = 'two cracked edges cannot be found, element is likely to be &
      &poorly shaped, crack_elem_centroid2d, global_toolkit_module'
      ! clean up intent out variable before error exit
      edge_crack_points = ZERO
      return
    end if

    ! update optional intent out dummy arg only before successful return
    if (present(crack_edge_indices)) crack_edge_indices = crack_edge_IDs

  end subroutine crack_elem_centroid2d




  pure subroutine assembleKF (Kmat, Fvec, Ki, Fi, cnc, istat, emsg)
  ! Purpose:
  ! to assemble elem/subelem K and F to system/xelem K and F
  ! variable sizes of dummy args have to be allowed here.
  ! explicit checking of dummy args' sizes must be performed
  ! with informative error messages flagged if unexpected sizes
  ! are encountered
  !
  ! the inputs must safisty the following conditions:
  ! - Kmat is square matrix with size = size of Fvec
  ! - Ki   is square matrix with size = size of Fi
  ! - size of cnc = size of Fi
  ! - min element of cnc must > 0
  ! - max element of cnc must < size of Fvec

  use parameter_module, only : DP, MSGLENGTH, STAT_SUCCESS, STAT_FAILURE

    real(DP), intent(inout) :: Kmat(:,:), Fvec(:)
    real(DP), intent(in)    :: Ki(:,:),   Fi(:)
    integer,  intent(in)    :: cnc(:)
    integer,                  intent(out)   :: istat
    character(len=MSGLENGTH), intent(out)   :: emsg

    ! local variable
    integer :: i, j, n, ni, mincnc, maxcnc

    ! initialize intent out and local variables
    istat  = STAT_SUCCESS
    emsg   = ''
    i      = 0
    j      = 0
    n      = 0
    ni     = 0
    mincnc = 0
    maxcnc = 0

    n      = size(Fvec)
    ni     = size(Fi)
    mincnc = minval(cnc)
    maxcnc = maxval(cnc)

    ! check validity of inputs

    if ( any (shape(Kmat) /= [n,n]) ) then
      istat = STAT_FAILURE
      emsg  = 'Kmat shape is incorrect, assembleKF, global_toolkit_module'
    else if ( any (shape(Ki) /= [ni,ni]) ) then
      istat = STAT_FAILURE
      emsg  = 'Ki shape is incorrect, assembleKF, global_toolkit_module'
    else if (size(cnc) /= ni) then
      istat = STAT_FAILURE
      emsg  = 'cnc size is incorrect, assembleKF, global_toolkit_module'
    else if (mincnc < 1) then
      istat = STAT_FAILURE
      emsg  = 'cnc min element <1, assembleKF, global_toolkit_module'
    else if (maxcnc > n) then
      istat = STAT_FAILURE
      emsg  = 'cnc max element is too large for Kmat, assembleKF, &
      &global_toolkit_module'
    end if

    if (istat == STAT_FAILURE) return

    ! proceed with the assembly only when all dummy arg sizes are consistent
    do i = 1, ni
      do j = 1, ni
          Kmat(cnc(j),cnc(i)) = Kmat(cnc(j),cnc(i)) + Ki(j,i)
      end do
      Fvec(cnc(i)) = Fvec(cnc(i)) + Fi(i)
    end do

  end subroutine assembleKF





end module global_toolkit_module
