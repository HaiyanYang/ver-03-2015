! fortran test program

! include useful modules

        include 'parameter_module.f90'
        include 'xgeometry_module.f90'

        program test
        use parameter_module
        use xgeometry_module

        implicit none

! create geometrical quantities
        type(xnode) :: xnodes(2)
        type(xedge) :: xedges(2)
        type(xtriangle) :: xtri(2)
        type(xquadrilateral) :: xquad(2)
        type(xtetrahedron) :: xtetra(2)
        type(xwedge) :: xwedges(2)
        type(xbrick) :: xbricks(2)


        integer :: i

        i=0


!! test xnode procedures
!
!        xnodes(2)=xnode(x=[1,2,3],u=0.1,&
!        & du=0.0,v=1.0,dof=[1.5,2.0],ddof=[0.0])
!
!
!        print*, 'xnode 2 x:',xnodes(2)%x
!        print*, 'xnode 2 u:',xnodes(2)%u
!        print*, 'xnode 2 du:',xnodes(2)%du
!        print*, 'xnode 2 v:',xnodes(2)%v
!        print*, 'xnode 2 dof:',xnodes(2)%dof
!        print*, 'xnode 2 ddof:',xnodes(2)%ddof
!
!        !call empty(xnodes(2))
!
!        call update(xnodes(2),x=[one,one,one],u=[one,one,one],&
!        & du=[zero,zero,zero],v=[zero,one,zero],dof=[ten,ten],ddof=[one,one])
!        ! must be exactly kind=dp, otherwise won't compile
!
!        print*, 'xnode 2 x:',xnodes(2)%x
!        print*, 'xnode 2 u:',xnodes(2)%u
!        print*, 'xnode 2 du:',xnodes(2)%du
!        print*, 'xnode 2 v:',xnodes(2)%v
!        print*, 'xnode 2 dof:',xnodes(2)%dof
!        print*, 'xnode 2 ddof:',xnodes(2)%ddof


!! test xedge procedures
!
!        xedges(2)=xedge(end_nodes=[1,3],floating_nodes=[1],&
!        & current_status=0,parent=1,children=[1,2])
!
!
!        print*, 'xedge 2 end nodes:',xedges(2)%end_nodes
!        print*, 'xedge 2 floating nodes:',xedges(2)%floating_nodes
!        print*, 'xedge 2 current status:',xedges(2)%current_status
!        print*, 'xedge 2 parent:',xedges(2)%parent
!        print*, 'xedge 2 children:',xedges(2)%children
!
!        !call empty(xedges(2))
!
!        call update(xedges(2),current_status=1,end_nodes=[10,25],floating_nodes=[3,7],&
!        & children=[2,6],parent=8)
!
!        print*, 'xedge 2 end nodes:',xedges(2)%end_nodes
!        print*, 'xedge 2 floating nodes:',xedges(2)%floating_nodes
!        print*, 'xedge 2 current status:',xedges(2)%current_status
!        print*, 'xedge 2 parent:',xedges(2)%parent
!        print*, 'xedge 2 children:',xedges(2)%children


!! test xtriangle procedures
!
!        xtri(2)=xtriangle(end_nodes=[1,3,9],edges=[2,5,8],floating_nodes=[1],&
!        & current_status=0,parent=1,children=[1,2])
!
!
!        print*, 'xtri 2 end nodes:',xtri(2)%end_nodes
!        print*, 'xtri 2 edges:',xtri(2)%edges
!        print*, 'xtri 2 floating nodes:',xtri(2)%floating_nodes
!        print*, 'xtri 2 current status:',xtri(2)%current_status
!        print*, 'xtri 2 parent:',xtri(2)%parent
!        print*, 'xtri 2 children:',xtri(2)%children
!
!        !call empty(xtri(2))
!
!        call update(xtri(2),current_status=1,end_nodes=[10,11,25],floating_nodes=[3,7],&
!        & children=[2,6,11],parent=8)
!
!        print*, 'xtri 2 end nodes:',xtri(2)%end_nodes
!        print*, 'xtri 2 edges:',xtri(2)%edges
!        print*, 'xtri 2 floating nodes:',xtri(2)%floating_nodes
!        print*, 'xtri 2 current status:',xtri(2)%current_status
!        print*, 'xtri 2 parent:',xtri(2)%parent
!        print*, 'xtri 2 children:',xtri(2)%children



!! test xquadrilateral procedures
!
!        xquad(2)=xquadrilateral(end_nodes=[1,3,9],edges=[2,5,8],floating_nodes=[1],&
!        & current_status=0,parent=1,children=[1,2])
!
!
!        print*, 'xquad 2 end nodes:',xquad(2)%end_nodes
!        print*, 'xquad 2 edges:',xquad(2)%edges
!        print*, 'xquad 2 floating nodes:',xquad(2)%floating_nodes
!        print*, 'xquad 2 current status:',xquad(2)%current_status
!        print*, 'xquad 2 parent:',xquad(2)%parent
!        print*, 'xquad 2 children:',xquad(2)%children
!
!        call empty(xquad(2))
!
!        call update(xquad(2),current_status=1,end_nodes=[10,11,25,1],floating_nodes=[3,7],&
!        & edges=[8,2,3,5],children=[2,6,11],parent=8)
!
!        print*, 'xquad 2 end nodes:',xquad(2)%end_nodes
!        print*, 'xquad 2 edges:',xquad(2)%edges
!        print*, 'xquad 2 floating nodes:',xquad(2)%floating_nodes
!        print*, 'xquad 2 current status:',xquad(2)%current_status
!        print*, 'xquad 2 parent:',xquad(2)%parent
!        print*, 'xquad 2 children:',xquad(2)%children



!! test xtetra procedures
!
!        xtetra(2)=xtetrahedron(end_nodes=[1,3,9],edges=[2,5,8],floating_nodes=[1],&
!        & tris=[11,32,65],current_status=0,parent=1,children=[1,2])
!
!
!        print*, 'xtetra 2 end nodes:',xtetra(2)%end_nodes
!        print*, 'xtetra 2 edges:',xtetra(2)%edges
!        print*, 'xtetra 2 tris:',xtetra(2)%tris
!        print*, 'xtetra 2 floating nodes:',xtetra(2)%floating_nodes
!        print*, 'xtetra 2 current status:',xtetra(2)%current_status
!        print*, 'xtetra 2 parent:',xtetra(2)%parent
!        print*, 'xtetra 2 children:',xtetra(2)%children
!
!        call empty(xtetra(2))
!
!        call update(xtetra(2),current_status=1,end_nodes=[10,11,25,1],floating_nodes=[3,7],&
!        & tris=[4,51,33,76],edges=[8,2,3,5,22,4],children=[2,6,11],parent=8)
!
!        print*, 'xtetra 2 end nodes:',xtetra(2)%end_nodes
!        print*, 'xtetra 2 edges:',xtetra(2)%edges
!        print*, 'xtetra 2 tris:',xtetra(2)%tris
!        print*, 'xtetra 2 floating nodes:',xtetra(2)%floating_nodes
!        print*, 'xtetra 2 current status:',xtetra(2)%current_status
!        print*, 'xtetra 2 parent:',xtetra(2)%parent
!        print*, 'xtetra 2 children:',xtetra(2)%children



!! test xwedge procedures
!
!        xwedges(2)=xwedge(end_nodes=[1,3,9],edges=[2,5,8],floating_nodes=[1],&
!        & tris=[11,32,65],quads=[21,62,72],current_status=0,parent=1,children=[1,2])
!
!
!        print*, 'xwedges 2 end nodes:',xwedges(2)%end_nodes
!        print*, 'xwedges 2 edges:',xwedges(2)%edges
!        print*, 'xwedges 2 tris:',xwedges(2)%tris
!        print*, 'xwedges 2 quads:',xwedges(2)%quads
!        print*, 'xwedges 2 floating nodes:',xwedges(2)%floating_nodes
!        print*, 'xwedges 2 current status:',xwedges(2)%current_status
!        print*, 'xwedges 2 parent:',xwedges(2)%parent
!        print*, 'xwedges 2 children:',xwedges(2)%children
!
!        call empty(xwedges(2))
!
!        call update(xwedges(2),current_status=1,end_nodes=[10,11,25,1,98,5],floating_nodes=[3,7],&
!        & tris=[4,51],edges=[8,2,3,5,22,4,3,12,7],children=[2,6,11],parent=8)
!
!        print*, 'xwedges 2 end nodes:',xwedges(2)%end_nodes
!        print*, 'xwedges 2 edges:',xwedges(2)%edges
!        print*, 'xwedges 2 tris:',xwedges(2)%tris
!        print*, 'xwedges 2 quads:',xwedges(2)%quads
!        print*, 'xwedges 2 floating nodes:',xwedges(2)%floating_nodes
!        print*, 'xwedges 2 current status:',xwedges(2)%current_status
!        print*, 'xwedges 2 parent:',xwedges(2)%parent
!        print*, 'xwedges 2 children:',xwedges(2)%children


! test xwedge procedures

        xbricks(2)=xbrick(end_nodes=[1,3,9],edges=[2,5,8],floating_nodes=[1],&
        & quads=[21,62,72],current_status=0,parent=1,children=[1,2])


        print*, 'xbricks 2 end nodes:',xbricks(2)%end_nodes
        print*, 'xbricks 2 edges:',xbricks(2)%edges
        print*, 'xbricks 2 quads:',xbricks(2)%quads
        print*, 'xbricks 2 floating nodes:',xbricks(2)%floating_nodes
        print*, 'xbricks 2 current status:',xbricks(2)%current_status
        print*, 'xbricks 2 parent:',xbricks(2)%parent
        print*, 'xbricks 2 children:',xbricks(2)%children

        call empty(xbricks(2))

        call update(xbricks(2),current_status=1,end_nodes=[10,11,25,1,98,5,9,12],floating_nodes=[3,7],&
        & edges=[8,2,3,5,22,4,3,12,7,12,4,9],children=[2,6,11],parent=8)

        print*, 'xbricks 2 end nodes:',xbricks(2)%end_nodes
        print*, 'xbricks 2 edges:',xbricks(2)%edges
        print*, 'xbricks 2 quads:',xbricks(2)%quads
        print*, 'xbricks 2 floating nodes:',xbricks(2)%floating_nodes
        print*, 'xbricks 2 current status:',xbricks(2)%current_status
        print*, 'xbricks 2 parent:',xbricks(2)%parent
        print*, 'xbricks 2 children:',xbricks(2)%children

        end program test
