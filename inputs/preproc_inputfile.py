################################################################
################# Preprocessing for the FNM ####################
##  writing fnm_nodes, fnm_edges and fnm_elems          ##
################################################################
##
##
################################################################
##  so far, only works on models with a single part, meshed
##  with 3D brick-type elements (C3D8(R) or SC8(R)); it prepares
##  a new input file and fnm modules with xlam element type only.
################################################################
##
##
################################################################
##  Applicable abaqus input file:
##  - single 3D part meshed with either C3D8(R) or SC8(R) elements
##  - part definition can only include *Node and *Element
##  - assembly definition can only include *Instance and *Nset
##  - boundary conditions can only be applied to assembly node sets
################################################################



# glb objects defined for FNM
from preproc_classes import*

import math
# operating system module
import os, sys
# shutil module has functions for different kinds of copying
import shutil


#***************************************************************
#   define parameters
#***************************************************************
# ndim    : dimension; currently only support 3D
# nprops  : no. of input material properties used in uel code (min 1)
# nsvars  : no. of sol. dpdnt var. used and ouput by uel code (min 1)
# uelcode : code for the user element fBrickLam, used in uel input file
# uellinelength : line length in writing uel input file (max. =256)
# fnmlinelength : line length in writing fnm input file (max. =132)
ndim    = 3
nprops  = 1
nsvars  = 1
uelcode = 308
uellinelength = 240
fnmlinelength = 120

#***************************************************************
#   fetch original input file and open uel files
#***************************************************************
# ask for job name
jobname      = raw_input('abaqus job name:')
# abaqus original input file
abqinputfile = jobname+'.inp'    
# uel input file     
uelinputfile = 'uel_'+jobname+'.inp'  
# uel nodes and elems files, to be included in uel input file
uelnodesfile = 'uel_nodes.inp'
uelelemsfile = 'uel_elems.inp'
# open files
abq_input = open(abqinputfile,'r')
uel_input = open(uelinputfile,'w')
uel_nodes = open(uelnodesfile,'w')
uel_elems = open(uelelemsfile,'w')

#***************************************************************
#   Open Fortran modules to be written during pre-processing
#***************************************************************
fnm_nodes = open('fnm_nodes.f90','w')  # list of all nodes
fnm_edges = open('fnm_edges.f90','w')  # list of all edges
fnm_elems = open('fnm_elems.f90','w')  # list of all elems


#***************************************************************
#       ask for layup from user
#***************************************************************
# rawlayup: layup of the laminate in a list of all plies
# blklayup: layup of the laminate in a list of plyblocks
rawlayup = []
blklayup = []

# ask for a list of ply angles
rawlayup = \
input('Enter fibre angles (int/float) of all plies, \
separated by comma, end with comma:')
# check if it is a list and if all elements in the list are numbers
while ( (not isinstance(rawlayup, (list,tuple))) or \
 (not all(isinstance(item, (int,float)) for item in rawlayup)) ):
    rawlayup = \
    input('Enter fibre angles (int/float) of all plies, \
    separated by comma, end with comma::')      

# ask for the thickness of a single ply
plythick = \
input('Enter the thickness of a single ply (positive real number):')
while ( not ( isinstance(plythick, float) and plythick > 0 ) ):
    plythick = \
    input('Enter the thickness of a single ply (positive real number):')
   
# find blocked plies and update blklayup
# initiate blklayup
blklayup.append( plyblk(angle=rawlayup[0], nplies=0, thickness=0.0) )
for plyangle in rawlayup:
    if (plyangle == blklayup[-1].angle):
        blklayup[-1].nplies += 1
        blklayup[-1].thickness += plythick
    else:
        blklayup.append( plyblk(angle=plyangle, nplies=1, thickness=plythick) )



#***************************************************************
#       Read Abaqus input file
#***************************************************************

# Store all lines first
All_lines = abq_input.readlines()

# Find the length of all lines
lenAll = len(All_lines)


#==================================================
# READ HEADER SECTION:
#==================================================
header  = []    # list of header lines
header.extend(All_lines[0:5])


#==================================================
# READ PARTS SECTION: 
# - Only ONE part is read
# - this part represents a single ply mesh
# - store part's real nodes, elems & its real nodes
# - find all breakable edges & create fl. nodes for each edge 
# - add fl. nodes to part's node list & elem's node list
# - add edges to part's edge list & elem's edge list
# (this piece of code works for input file with multiple parts)
# ** NOTE **: 
# - for now, only ONE type of elem is supported; all elems will be FNM-ized
# - for now, only reads nodes and elems. nset, elset & surfaces can be included
#   in the future 
#==================================================
parts   = []    # list of all parts in the model
botnds  = []    # list of all nodes on the bot surf.
topnds  = []    # list of all nodes on the bot surf.

# find the line no. of *Part and store in jparts
jparts = [j for j,line in enumerate(All_lines) if '*Part' in line]
if (len(jparts)>1):
    print("ERROR: more than one part is not yet supported!")
    sys.exit()

for jp in jparts:

    # read Part name
    pname = All_lines[jp][12:].rstrip()
    
    # create a new part in parts list
    parts.append(part(name=pname, nodes=[], NtN=[], edges=[], elems=[]))
    
    # find the line of *Node and *Element
    # ** NOTE: only ONE type of elem is supported; all elems will be FNM-ized
    jnode = next( j for j in range(jp,lenAll) if '*Node'    in All_lines[j] )
    jelem = next( j for j in range(jp,lenAll) if '*Element' in All_lines[j] )
    
    # read nodes of this part 
    for jn in range(jnode+1,lenAll):
        nline = All_lines[jn]
        # break out of for loop if end of node section is reached
        if('*' in nline):
            break
        # read the coords of this node into a list of float numbers
        coords = []
        for t in nline.split(','):
            try:
                coords.append(float(t))
            except ValueError:
                pass
        # store the coords in nodes list of this part
        parts[-1].nodes.append(node(x=coords[1], y=coords[2], z=coords[3]))
        # update the node-to-node (NtN) matrix row length
        parts[-1].NtN.append([])
    ##check node correctness
    #for nd in parts[-1].nodes:
    #    print(str(nd.x)+','+str(nd.y)+','+str(nd.z))
    
    # append the NtN matrix column to the correct length (nnode)
    for r in range(len(parts[-1].NtN)):
        parts[-1].NtN[r] = [0]*len(parts[-1].NtN)
        
    # read elems of this part
    for je in range(jelem+1,lenAll):
        eline = All_lines[je]
        # break out of for loop if end of elem section is reached
        if('*' in eline):
            break
        # read the index and real nodes of this elem into a list of int numbers
        el = []
        for t in eline.split(','):
            try:
                el.append(int(t))
            except ValueError:
                pass
        id  = el[0]  # elem index
        nds = el[1:] # elem real nodes
        # store the index and real nodes in elem list of this part
        parts[-1].elems.append(element(index=id, nodes=nds, edges=[]))
        # update the top and bot surf nodes lists
        halfnds = len(nds)/2
        for j in range(halfnds):
            # bot surf nodes
            if not(nds[j] in botnds):
                botnds.append(nds[j])
            if not(nds[halfnds+j] in topnds):
                topnds.append(nds[halfnds+j])
        # form edges
        # in 3D FNM for composites, only edges parallel to shell plane are breakable
        # j = 0: bottom edges; 
        # j = 1: top edges
        for j in range(2):
            # j = 0: loop over nodes on bot surf
            # j = 1: loop over nodes on top surf
            for i in range(halfnds):
                # ith node on surf j
                row = nds[ j*halfnds + i ] - 1
                # (i+1)th node on surf j; 
                # if i is the last node, i+1 is the first node
                if i == halfnds-1:
                    col = nds[ j*halfnds ] - 1
                else:
                    col = nds[ j*halfnds + i + 1 ] - 1
                # fill in the NtN matrix:
                # fill in the edge index composed of the two end nodes
                if parts[-1].NtN[row][col]==0:
                # this pair of nodes hasn't formed an edge
                # a new edge will be formed
                    # indices of 2 fl. nodes on this new edge
                    # they are two new nodes to be added in the part's node list
                    fn1 = len(parts[-1].nodes)
                    fn2 = fn1 + 1
                    parts[-1].nodes.append( node(0.0,0.0,0.0) )
                    parts[-1].nodes.append( node(0.0,0.0,0.0) )
                    # form a new edge and append to existing list of edges 
                    # the new edge has 4 nodes: 2 real, 2 floating
                    parts[-1].edges.append(edge(nodes=[row+1,col+1,fn1+1,fn2+1]))
                    # fill the new edge index in the NtN matrix
                    # nodes in rev. order makes the same edge in rev. dir.
                    parts[-1].NtN[row][col] =  len(parts[-1].edges)
                    parts[-1].NtN[col][row] = -len(parts[-1].edges)
                    # append this edge no. in this elem
                    parts[-1].elems[-1].edges.append(parts[-1].NtN[row][col]) 
                    # append the fl. nodes in this elem
                    parts[-1].elems[-1].nodes.extend([fn1+1,fn2+1])
                else:
                # this pair of nodes has already formed an edge
                    eg = parts[-1].NtN[row][col]
                    # append this edge no. in this elem 
                    parts[-1].elems[-1].edges.append(eg)
                    # find the two fl. nodes on this edge
                    fn1 = parts[-1].edges[abs(eg)].nodes[2]
                    fn2 = parts[-1].edges[abs(eg)].nodes[3]
                    # append the fl. nodes in this elem
                    # edge is in the same order as saved
                    if eg > 0:
                        parts[-1].elems[-1].nodes.extend([fn1,fn2])
                    # edge is in the rev. order as saved
                    else:
                        parts[-1].elems[-1].nodes.extend([fn2,fn1])
    
    #check elem correctness
    for el in parts[-1].elems:
        print(str(el.index)+','+str(el.nodes)+','+str(el.edges))
    ##check NtN
    #print(str(parts[-1].NtN))


#==================================================
# read assembly section:
# - only a single assembly is supported
# - only a single instance is supported
# - multiple nsets are supported
# - 'generate' operation in nset is supported
# - no translation or similar operations on a instance is supported
# - elset is currently not yet supported
#==================================================
assemblies = [] # list of all assemblies in the model

# find the line no. of *Assembly and store in jassemblies
jassemblies = [j for j,line in enumerate(All_lines) if '*Assembly' in line]
if (len(jassemblies)>1):
    print("ERROR: more than one assembly is not yet supported!")
    sys.exit()

for ja in jassemblies:

    # read assembly name
    aname = All_lines[ja].rstrip()
    
    # create a new assembly in the assembly list
    assemblies.append(assembly(name=aname, instances=[], nsets=[], elsets=[]))
    
    # find the lines of *Instance, *Nset and *Elset
    jinstances = [ j for j in range(ja,lenAll) if '*Instance' in All_lines[j] ]
    jnsets     = [ j for j in range(ja,lenAll) if '*Nset'     in All_lines[j] ]
    jelsets    = [ j for j in range(ja,lenAll) if '*Elset'    in All_lines[j] ]
    if (len(jinstances)>1):
        print("ERROR: more than one instance in an assembly is not yet supported!")
        sys.exit()
        
    # read instances in this assembly
    for jin in jinstances:
        # add this instance in the list of instances in this assembly
        assemblies[-1].instances.append( instance( name=All_lines[jin].rstrip() ) )
        
    # read nsets in this assembly
    for jns in jnsets:
        nline = All_lines[jns].rstrip()
        # remove 'generate' in the line if present
        if ('generate' in All_lines[jns]):
            nline = All_lines[jns][0:-11]
        # add this nset in the list of nsets in this assembly
        assemblies[-1].nsets.append( nset( name=nline, rnodes=[], edges=[] ) )
        # read nodes in the nset
        # if generate is used, then calculate all nodes;
        # otherwise, read all nodes directly
        if ('generate' in All_lines[jns]):
            nline = All_lines[jns+1]
            nl = []
            for t in nline.split(','):
                try:
                    nl.append(int(t))
                except ValueError:
                    pass
            nds = nl[0] # start node
            ndf = nl[1] # final node
            try:
                itv = nl[2] # interval
            except IndexError:
                itv = 1
            for n in range(nds,ndf+1,itv):
                assemblies[-1].nsets[-1].rnodes.append(n)
        else:
            # read the lines of nodes in this nset
            nl = [] # list of node to be filled
            for n in range(jns+1,lenAll):
                nline = All_lines[n]
                # break out of loop if end of section encountered
                if ('*' in nline):
                    break
                for t in nline.split(','):
                    try:
                        nl.append(int(t))
                    except ValueError:
                        pass
            assemblies[-1].nsets[-1].rnodes.extend(nl)
        # find the edges involved in this nset, 
        # and include the fl. nodes in the nset
        # extract this nset from assembly
        nst = assemblies[-1].nsets[-1]
        # find the part involved in this nset using part name
        for prt in parts:
            if (prt.name in nst.name):
                # loop over all node pairs in this nset
                for n1 in range(len(nst.rnodes)-1):
                    for n2 in range(n1+1,len(nst.rnodes)):
                        rnd = nst.rnodes[n1]-1
                        cnd = nst.rnodes[n2]-1
                        # if this node pair forms an edge
                        if (prt.NtN[rnd][cnd]!=0):
                            # get this edge number
                            eg = abs(prt.NtN[rnd][cnd])
                            #print(' node '+str(rnd)+' node '+str(cnd)+' forms edge '+str(eg))
                            # store this edge in the nset
                            nst.edges.append(eg)
                # after finding all the edges, break out of for loop
                break
        # update this nset in assembly
        assemblies[-1].nsets[-1] = nst 
    
    # read elsets in this assembly (NOT YET SUPPORTED)

    # check assembly correctness
    print(assemblies[-1].name)
    print(assemblies[-1].instances[-1].name)
    for m in assemblies[-1].nsets:
        print(m.name)
        print(str(m.rnodes))
        print(str(m.edges))


#==================================================
# read (fixed) boundary section:
#==================================================
bcds  = []    # list of all boundary conditions

# find the separation line '** -------------'
jdash = next(j for j,line in enumerate(All_lines) \
if '** ----------------------------------------------------------------' in line)

# find the lines with *boundary
jbcds = [j for j in range(0,jdash) if '*Boundary' in All_lines[j]]

# loop over all bcds, store them without modification
for jb in jbcds:
    for k in range(jb+1,jdash):
        bline = All_lines[k]
        if ('*' in bline):
            break
        bcds.append(bline)
        # find the nsets involved in this bline (future)
        # find the edges involved in this bline (future)
print(bcds)


#==================================================
# read step
# store all the lines from jdash(inc.)
#==================================================
step = []
step.extend(All_lines[jdash:])
print(step)




#***************************************************************
#       write nodes
#***************************************************************   
# find no. of plyblocks
nplyblk = len(blklayup)
# find no. of nodes and edges in a ply of this mesh
nnode_p = len(parts[0].nodes)
nedge_p = len(parts[0].edges)
# find internal nodes for an interface of this mesh
nnodein = nedge_p
# find the total no. of nodes in this mesh
nnodett = nplyblk * nnode_p + (nplyblk-1) * nnodein

# write fnm_nodes.f90 header
fnm_nodes.write('subroutine fnm_nodes()                \n')
fnm_nodes.write('use parameter_module, only: DP, ZERO  \n')
fnm_nodes.write('use node_list_module, only: node_list \n')
fnm_nodes.write('use fnode_module,     only: update    \n')
fnm_nodes.write('                                      \n')
fnm_nodes.write('  integer :: nnode=0                  \n')   
fnm_nodes.write('  integer :: i=0                      \n')
fnm_nodes.write('                                      \n')
fnm_nodes.write('  nnode='+str(nnodett)+'              \n')
fnm_nodes.write('  allocate(node_list(nnode))          \n')

for ipb in range(nplyblk):
    # calculate the bot and top real node z-coordinate
    # zbot = 0 if this is the 1st plyblk
    if (ipb == 0):
        zbot = 0
    # zbot = thickness of last plyblk
    else:
        zbot = zbot + blklayup[ipb-1].thickness
    # ztop = zbot + thickness of this plyblk
    ztop = zbot + blklayup[ipb].thickness
    # loop over all nodes in the single ply mesh
    for cntr0, nd in enumerate(parts[0].nodes):
        # current node id of the node on the ith plyblk
        cntr = ipb * nnode_p + (cntr0+1)
        # check if this node is a real node on the bot/top surf, or a fl. node
        # must use cntr0+1 as bot/topnds lists are for nodes in 1st plyblk
        if ((cntr0+1) in botnds):
            zz = zbot
        elif ((cntr0+1) in topnds):
            zz = ztop
        else:
            zz = 0.0
        # write this node coords in uel_nodes.inp
        uel_nodes.write\
        (str(cntr)+', '+str(nd.x)+', '+str(nd.y)+', '+str(zz)+'\n')
        # write this node coords in fnm node_list
        fnm_nodes.write\
        ('  call update(node_list('+str(cntr)+'),\
        x=['+str(nd.x)+'_DP,'+str(nd.y)+'_DP,'+str(zz)+'_DP],\
        u=[ZERO,ZERO,ZERO])\n')

# write the additional nodes of interfaces if they're present
if (nplyblk > 1):
    for jintf in range(nplyblk-1):
        # loop over all edges in the base ply mesh
        # each edge has one additional node 
        for cntr0 in range(nedge_p):
            cntr = nplyblk * nnode_p + jintf * nedge_p + cntr0 + 1
            # write this node coords in uel_nodes.inp
            uel_nodes.write\
            (str(cntr)+', 0.0, 0.0, 0.0\n')
            # write this node coords in fnm node_list
            fnm_nodes.write\
            ('  call update(node_list('+str(cntr)+'),\
            x=[ZERO,ZERO,ZERO],\
            u=[ZERO,ZERO,ZERO])\n')
            
fnm_nodes.write('\n')
fnm_nodes.write('end subroutine fnm_nodes\n')

#***************************************************************
#       write edges
#***************************************************************
# find the total no. of edges in this mesh
nedgett = nplyblk * nedge_p
# write fnm_edges.f90
fnm_edges.write('subroutine fnm_edges()            \n')
fnm_edges.write('use edge_list_module, only: edge_list \n')
fnm_edges.write('                                  \n')
fnm_edges.write('  integer :: nedge=0              \n')
fnm_edges.write('                                  \n')
fnm_edges.write('  nedge='+str(nedgett)+'          \n')
fnm_edges.write('  allocate(edge_list(nedge))      \n')
fnm_edges.write('  edge_list = 0                   \n')
fnm_edges.write('                                  \n')
fnm_edges.write('end subroutine fnm_edges\n')

#***************************************************************
#       write elems
#***************************************************************   
# find no. of elems in a ply mesh
nelem_p = len(parts[0].elems) 
# find total no. of elems in the laminate
# it is the same as nelem_p, as a fnm elem contains all plies&interfs
nelemtt = nelem_p
# find the no. of r+f nodes in an elem of a single plyblk
elnndrf_p = len(parts[0].elems[0].nodes)
# find the no. of edges in an elem of a single plyblk
elnedge_p = len(parts[0].elems[0].edges)
# find the no. of r+f nodes in an elem of the laminate
elnndrf_l = elnndrf_p * nplyblk
# find the no. of interface internal nodes in an elem of the laminate
elnndin_l = elnedge_p * (nplyblk-1)
# find the total no. of nodes in an elem of the laminate
elnndtt_l = elnndrf_l + elnndin_l
# find the no. of edges in an elem of the laminate
elnedge_l = elnedge_p * nplyblk

fnm_elems.write('subroutine fnm_elems()                                     \n')
fnm_elems.write('use parameter_module,      only: DP                        \n')
fnm_elems.write('use elem_list_module,      only: elem_list,&               \n') 
fnm_elems.write('                      & elem_node_connec, elem_edge_connec \n')
fnm_elems.write('use fBrickLam_elem_module, only: plyblock_layup, set       \n') 
fnm_elems.write('                                                           \n') 
fnm_elems.write('  integer :: nelem   = 0                                   \n') 
fnm_elems.write('  integer :: elnnode = 0                                   \n') 
fnm_elems.write('  integer :: elnedge = 0                                   \n')  
fnm_elems.write('  integer :: nplyblk = 0                                   \n')
fnm_elems.write('  integer, allocatable :: nodecnc(:), edgecnc(:)           \n')
fnm_elems.write('  type(plyblock_layup), allocatable :: layup(:)            \n')
fnm_elems.write('                                                           \n')
fnm_elems.write('  nelem   ='+str(nelemtt)+'                                \n')
fnm_elems.write('  elnnode ='+str(elnndtt_l)+'                              \n')
fnm_elems.write('  elnedge ='+str(elnedge_l)+'                              \n')
fnm_elems.write('  nplyblk ='+str(nplyblk)+'                                \n')
fnm_elems.write('  allocate(elem_list(nelem))                               \n')
fnm_elems.write('  allocate(elem_node_connec(elnnode,nelem))                \n')
fnm_elems.write('  allocate(elem_edge_connec(elnedge,nelem))                \n')
fnm_elems.write('  allocate(nodecnc(elnnode))                               \n')
fnm_elems.write('  allocate(edgecnc(elnedge))                               \n')
fnm_elems.write('  allocate(layup(nplyblk))                                 \n')
fnm_elems.write('  nodecnc = 0                                              \n')
fnm_elems.write('  edgecnc = 0                                              \n')
fnm_elems.write('                                                           \n')
# write layup array
for jpb in range(nplyblk):
    angle  = str(blklayup[jpb].angle)
    nplies = str(blklayup[jpb].nplies)
    if '.' in angle:
        angle = angle+'_DP'
    else:
        angle = angle+'._DP'
    fnm_elems.write(' layup('+str(jpb+1)+')=plyblock_layup(angle='+angle+',nplies='+nplies+') \n')
fnm_elems.write('                                                           \n')

for jel in range(nelemtt):
    elnds_p = []
    elegs_p = []
    elnds_l = []
    elegs_l = []
    # find the r+f nodes in this elem ply
    elnds_p = parts[0].elems[jel].nodes
    # find the edges in this elem ply
    elegs_p = parts[0].elems[jel].edges
    # find the r+f nodes & edges in this elem laminate; 
    # abs is needed to remove the negative sign on some of the edge no.
    for jpb in range(nplyblk):
        elnds_l.extend( [     x  + nnode_p * jpb  for x in elnds_p ] )
        elegs_l.extend( [ abs(x) + nedge_p * jpb  for x in elegs_p ] )
    # find the internal nodes of interfs in this elem laminate
    # intern. nodes are listed after r and f nodes in the node list
    # they are assigned to each edges of the interfaces
    # so the elem's edge connec is used for assignment of intern. nodes
    if nplyblk > 1:
        for jit in range(nplyblk-1):
            elnds_l.extend( [ x + nnode_p * nplyblk + nedge_p * jit for x in elegs_p ] )
            
    #**** write elem's nodal and edge connec to uel&fnm_elems ****
    
    #** node cnc 
    # start the line with elem index jel+1
    eline = [str(jel+1)+',']  # node cnc dataline for uel_elems
    fline = ['']              # node cnc dataline for fnm_elems
    # add the node no. to the line one by one
    for k in elnds_l:
        # if the uel line gets too long, continue on next line
        if (len(eline[-1]+str(k)) >= uellinelength):
            eline.append('')
        # add the node no. to the line
        eline[-1] = eline[-1]+str(k)+','
        # if the fnm line gets too long, continue on the next line
        if (len(fline[-1]+str(k)) >= fnmlinelength):
            fline.append('')
        # add the node no. to the line
        fline[-1] = fline[-1]+str(k)+','
    # remove the last comma from the eline
    eline[-1] = eline[-1][:-1]
    # remove the last comma from the fline
    fline[-1] = fline[-1][:-1]
    
    #** edge cnc
    gline = ['']  # edge cnc dataline for fnm_elems
    # add the node no. to the line one by one
    for k in elegs_l:
        # if the fnm line gets too long, continue on the next line
        if (len(gline[-1]+str(k)) >= fnmlinelength):
            gline.append('')
        # add the node no. to the line
        gline[-1] = gline[-1]+str(k)+','
    # remove the last comma from the fline
    gline[-1] = gline[-1][:-1]
    
    # write the line of elem node connec
    for l in eline:
        uel_elems.write(l+'\n')
        
    # set this elem in fnm_elems subroutine
    fnm_elems.write('\n')
    # write nodecnc array
    fnm_elems.write('  nodecnc=[ &\n')
    for l in fline:
        fnm_elems.write('& '+l+' &\n')
    fnm_elems.write('& ]\n')
    # write edgecnc array
    fnm_elems.write('  edgecnc=[ &\n')
    for l in gline:
        fnm_elems.write('& '+l+' &\n')
    fnm_elems.write('& ]\n')
    fnm_elems.write('  call set(elem_list('+str(jel+1)+'), NPLYBLKS='+str(nplyblk)+',& \n')
    fnm_elems.write('& node_connec=nodecnc, layup=layup)\n')
    # write elem_node_connec array
    fnm_elems.write('  elem_node_connec(:,'+str(jel+1)+')=nodecnc(:)\n')
    # write elem_edge_connec array
    fnm_elems.write('  elem_edge_connec(:,'+str(jel+1)+')=edgecnc(:)\n')
    fnm_elems.write('\n')

fnm_elems.write('end subroutine fnm_elems\n')


#***************************************************************
#       write uel input file
#*************************************************************** 
#**** write HEADER ****
for hline in header:
    uel_input.write(str(hline[0:]))

#**** write PART   ****
# part name 
uel_input.write('*Part, name='+parts[0].name+'\n')
# part nodes
uel_input.write('*NODE,INPUT=uel_nodes.inp    \n')
# user element definition
uel_input.write('*USER ELEMENT, TYPE=U'+str(uelcode)+\
', NODES='+str(elnndtt_l)+', COORDINATES='+str(ndim)+\
', PROPERTIES='+str(nprops)+', VARIABLES='+str(nsvars)+'\n')
uel_input.write('1,2,3\n')
# elements and fnm elset
uel_input.write('*ELEMENT, TYPE=U'+str(uelcode)+', ELSET=fnm, INPUT=uel_elems.inp \n')
# write the mandatory uel property line (not needed for calculation)
uel_input.write('*UEL PROPERTY, ELSET=fnm\n')
uel_input.write('1\n')
uel_input.write('*End Part\n')

#**** write ASSEMBLY ****
assbl = assemblies[0]
# assembly name line
uel_input.write(assbl.name+'\n')
# assembly instance lines
for inst in assbl.instances:
    uel_input.write(inst.name+'\n')
    uel_input.write('*End Instance\n')
# assembly nsets
for nst in assbl.nsets:
    # write the nst name
    uel_input.write(nst.name+'\n')
    # nst dataline for uel_input, to be filled
    nstline = ['']
    # find the part this nst is based on
    for prt in parts:
        if (prt.name in nst.name):
            nstprt = prt
    # if all the real nodes in this nst are on the bot surface, then
    # only store the bot plyblk nodes, DO NOT include the other plies
    if all(n in botnds for n in nst.rnodes):
        pstart = 0
        pend   = 1
    # if all the real nodes in this nst are on the top surface, then
    # only store the top plyblk nodes, DO NOT include the other plies
    elif all(n in topnds for n in nst.rnodes):
        pstart = nplyblk-1
        pend   = nplyblk
    else:
        pstart = 0
        pend   = nplyblk
    # find nst nodes in plyblks
    for jpb in range(pstart,pend):
        # add the real nodes to the list one by one
        for n in nst.rnodes:
            # find the corresponding node on the jpb-th plyblk
            k = n + jpb * nnode_p
            # if the uel line gets too long, continue on the next line
            if (len(nstline[-1]+str(k)) >= uellinelength):
                nstline.append('')
            # add the node no. to the line
            nstline[-1] = nstline[-1]+str(k)+','
        # add the fl. nodes to the list one by one
        for eg in nst.edges:
            k1 = nstprt.edges[eg-1].nodes[2] + jpb * nnode_p
            k2 = nstprt.edges[eg-1].nodes[3] + jpb * nnode_p
            # if the uel line gets too long, continue on the next line
            if (len(nstline[-1]+str(k1)+str(k2)) >= uellinelength):
                nstline.append('')
            # add the nodes to the line
            nstline[-1] = nstline[-1]+str(k1)+','+str(k2)+','
    # remove the last comma from the line
    nstline[-1] = nstline[-1][:-1]
    # write all original nodes of the nset
    for nstl in nstline:
        uel_input.write(nstl+'\n')
# end assembly writing
uel_input.write('*End Assembly\n')

#**** write fixed BCDs ****
for bline in bcds:
    uel_input.write('*Boundary\n')
    uel_input.write(str(bline[0:]))

#**** write step ****
for sline in step:
    uel_input.write(str(sline[0:]))



#*************************************************************** 
# close all open files
#*************************************************************** 
#   close input files
abq_input.close()
uel_input.close()
#   close nodes files
fnm_nodes.close()
uel_nodes.close()
#   close edges file
fnm_edges.close()
#   close elems files
fnm_elems.close()
uel_elems.close()


#*************************************************************** 
# copy fnm input file to main directory
#*************************************************************** 

# get current working directory
cwd = os.getcwd()
# parent directory of cwd
pwd = os.path.dirname(cwd)
# copy fnm input file to parent directory of preprocessing directory (which is assumed to be the working directory)
shutil.copy (uelinputfile,pwd)
