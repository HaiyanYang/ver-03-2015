    ! variables to be included in declaration:
    ! - max_proj          : max projection of line segments on the crack line
    ! - line_vect         : a line segment vector
    ! - line_proj         : projection of a line segment on the crack line
    real(DP) :: max_proj, line_vect(2), line_proj


    !**** RARE ADDITIONAL CALCULATIONS ****

    ! badly shaped element may have large angles. the numerical error may
    ! prevent the algorithm from finding any broken edge or only ONE broken edge

    ! if only ONE cracked edge is found by the above algorithm, then this
    ! indicates that the crack line is almost parallel to all of the rest
    ! of the edges. So you can just pick a midpoint of one edge to be the second
    ! cross point with any criterion.
    if (n_crack_edge == 1) then 
    ! find the edge whose midpoint and the exist edge crack point forms a
    ! line which has the largest projection on the crack line
 
        ! initialize max proj. on crack line
        max_proj = - SMALLNUM
        do i = 1, nedge
          ! cycle to next round if edge i is the existing cracked edge
          if (i == crack_edge_IDs(1)) cycle
          ! find the two end nodes' coords of edge i
          coords_edge(:, 1) = coords(:, nodes_on_edge(1,i))
          coords_edge(:, 2) = coords(:, nodes_on_edge(2,i))
          ! find mid point of edge i as a potential cross point
          cross_point(:)    = HALF * (coords_edge(:,1) + coords_edge(:,2))
          ! line vector from this cross point to the existing crack point
          line_vect(:)      = edge_crack_points(:,1) - cross_point(:)
          ! projection of this line vector on the crack line
          line_proj         = abs(dot_product(line_vect, crack_unit_vect))
          ! find the trial line with the most projection and update stored info
          if (line_proj > max_proj) then
            ! update no. of cracked edges
            n_crack_edge = 2
            ! update the index of the 2nd cracked edge
            crack_edge_IDs(2) = i
            ! update the coords of the 2nd crack point
            edge_crack_points(:,2) = cross_point(:)
            ! update max projection
            max_proj = line_proj
          end if
        end do

    ! if NO cracked edge is found by the above algorithm, then this
    ! indicates that the crack line is almost parallel to all edges
    ! So just pick the midpoints of two NEIGHBOUR edges to be the 
    ! two crack points.
    else if (n_crack_edge == 0) then
    ! find the midpoints of TWO NEIGHBOUR edges

        ! the following algorithm will guarantee to find TWO broken edges
        n_crack_edge = 2
        ! find midpoints of the first two edges
        do i = 1, 2
          ! find the two end nodes' coords of edge i
          coords_edge(:, 1) = coords(:, nodes_on_edge(1,i))
          coords_edge(:, 2) = coords(:, nodes_on_edge(2,i))
          ! find mid point of edge i 
          cross_point(:)    = HALF * (coords_edge(:,1) + coords_edge(:,2))
          ! store the mid point in edge_crack_points
          edge_crack_points(:, i) = cross_point(:)
        end do
        ! the indices of the two cracked edges are simply the first two
        crack_edge_IDs(1:2) = [1,2]

    end if

    !**** END RARE ADDITIONAL CALCULATIONS ****