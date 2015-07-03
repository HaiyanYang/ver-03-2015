################################################################
################# Preprocessing for the FNM ####################
##  writing input_nodes, input_edges and input_elems          ##
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
from fnm_classes import*

import math
# operating system module
import os, sys
# shutil module has functions for different kinds of copying
import shutil




#***************************************************************
#   fetch Abaqus input files
#***************************************************************
jobname      = raw_input('abaqus job name:')
abqinputfile = jobname+'.inp'         # abaqus input file
fnminputfile = 'fnm-'+jobname+'.inp'  # fnm uel input file




#***************************************************************
#   define model information
#***************************************************************
# non-useful variables defined to conform with abaqus formats:
# nprops: no. of input material properties used in uel code (min 1)
# nsvars: no. of sol. dpdnt var. used and ouput by uel code (min 1)
nprops = 1
nsvars = 1




#***************************************************************
#   Open Fortran modules to be written during pre-processing
#***************************************************************
input_nodes = open('input_nodes.f90','w')  # list of all nodes
input_edges = open('input_edges.f90','w')  # list of all edges
input_elems = open('input_elems.f90','w')  # list of all elems




#***************************************************************
#       open input files for I/O
#***************************************************************
abqinp = open(abqinputfile,'r')
fnminp = open(fnminputfile,'w')




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
    
# find blocked plies and update blklayup
blklayup.append( plyblk(angle=rawlayup[0], nplies=0) ) # initiate blklayup
for plyangle in rawlayup:
    if (plyangle == blklayup[-1].angle):
        blklayup[-1].nplies += 1
    else:
        blklayup.append( plyblk(angle=plyangle, nplies=1) )

# put layup info in a string for later output purpose
layupstr=''
for lp in blklayup:
    if ('.' in str(lp.angle)):
        layupstr = layupstr+str(lp.angle)+'_dp,'
    else:
        layupstr = layupstr+str(lp.angle)+'._dp,'        
    layupstr = layupstr+str(lp.nplies)+','
# print blocked layup 
print("The layup is: "+layupstr)




#***************************************************************
#       Read Abaqus input file
#***************************************************************

# Store all lines first
All_lines = abqinp.readlines()

# Find the length of all lines
lenAll = len(All_lines)


#==================================================
# READ HEADER SECTION:
#==================================================
header  = []    # list of header lines
header.extend(All_lines[0:5])


#==================================================
# READ PARTS SECTION: 
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
NtN     = []    # matrix of Node-to-Node edges (auxilliary)

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
        NtN.append([])
    ##check node correctness
    #for nd in parts[-1].nodes:
    #    print(str(nd.x)+','+str(nd.y)+','+str(nd.z))
    
    # append the NtN matrix column to the correct length (nnode)
    for r in range(len(NtN)):
        NtN[r] = [0]*len(NtN)
        
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
        parts[-1].elems.append(element(index=id, rnodes=nds, edges=[]))
        # form edges
        # in 3D FNM for composites, only edges parallel to shell plane are breakable
        # j = 0: bottom edges; 
        # j = 1: top edges
        halfnds = len(nds)/2
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
                if NtN[row][col]==0:
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
                    NtN[row][col] =  len(parts[-1].edges)
                    NtN[col][row] = -len(parts[-1].edges)
                    # append this edge no. in this elem
                    parts[-1].elems[-1].edges.append(NtN[row][col]) 
                else:
                # this pair of nodes has already formed an edge
                    # append this edge no. in this elem 
                    parts[-1].elems[-1].edges.append(NtN[row][col])           
    
    # store NtN matrix to the part component NtN
    parts[-1].NtN = NtN
    
    ##check elem correctness
    #for el in parts[-1].elems:
    #    print(str(el.index)+','+str(el.rnodes)+','+str(el.edges))
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

    ## check assembly correctness
    #print(assemblies[-1].name)
    #print(assemblies[-1].instances[-1].name)
    #for m in assemblies[-1].nsets:
    #    print(m.name)
    #    print(str(m.rnodes))
    #    print(str(m.edges))


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
#print(bcds)


#==================================================
# read step
# store all the lines from jdash(inc.)
#==================================================
step = []
step.extend(All_lines[jdash:])
#print(step)



nplyblk = len(blklayup)
nndr_p  = len(parts[0].rnodes)
nedge_p = len(parts[0].edges) 
nelem_p = len(parts[0].elems)
nndf_p  = 2*nedge_p
nndi_p  = nedge_p
nndt_p  = nndr_p + nndf_p + nndi_p

#***************************************************************
#       write input_nodes.f90
#***************************************************************   
input_nodes.write('subroutine input_nodes()              \n')
input_nodes.write('use node_list_module, only: node_list \n')
input_nodes.write('use fnode_module,     only: update    \n')
input_nodes.write('                                      \n')
input_nodes.write('  integer :: nnode=0                  \n')   
input_nodes.write('  integer :: i=0                      \n')
input_nodes.write('\n')
input_nodes.write('  nnode='+str(nplyblk*nndt_p)+'       \n')
input_nodes.write('  allocate(node_list(nnode))          \n')
# duplicate the real and fl. nodes 
for iply in range(nplyblk):   
    for cntr0, nd in enumerate(parts[0].nodes):
        cntr = cntr0 + 1 + iply * ( nndr_p + nndf_p )
        fnminp.write(str(cntr)+', '+str(nd.x)+', '+str(nd.y)+', '+str(nd.z)+'\n')
        input_nodes.write('  call update(node_list('+str(cntr)+'),x=['+str(nd.x)+'_dp,'+str(nd.y)+'_dp,'+str(nd.z)+'_dp],u=[zero,zero,zero])\n')

input_nodes.write('\n')

#***************************************************************
#       write input_edges.f90 common codes
#***************************************************************
input_edges.write('subroutine input_edges()          \n')
input_edges.write('                                  \n')
input_edges.write('  integer :: nedge=0              \n')
input_edges.write('  integer :: i=0                  \n')
input_edges.write('\n')
#***************************************************************
#       write input_elems.f90 common codes
#***************************************************************    
input_elems.write('subroutine input_elems()          \n')
input_elems.write('                                  \n') 
input_elems.write('  integer :: nelem=0              \n')
input_elems.write('  integer :: i=0                  \n')
input_elems.write('\n')





#
##***************************************************************        
##       Write Edges
##***************************************************************
#input_edges.write('        nedge='+str(nplyblk*plynedge)+'   \n')
#input_edges.write('        allocate(input_edges(nedge))   \n')
#input_edges.write('        input_edges=0                  \n')
#
#    
##***************************************************************        
##       Write Elements
##***************************************************************
## get total no. of elements in the mesh
#nelem=sum(elcount)
## allocate input_elems accordingly
#input_elems.write('        nelem='+str(nelem)+'            \n')
#input_elems.write('        allocate(input_elems(nelem))       \n')
#
#
## determine no. of nodes and uel element integer type
## uel element integer type:
##   1st digit: ndim
##   2nd digit: 0: bulk elem; 1: coh elem; 2: x(brick/wedge) elem; 3: xlam elem
##   3rd onwards: no. of real nodes
#if ndim == 2:
#        print 'error: fnm 2D elem type not supported for use!'
#elif ndim == 3:      
#    if elt == 'xlam':
#        nndrl=8*nplyblk         # no. of (real)nodes in the element. (here, 8 for linear brick elmt)
#        nedge=8*nplyblk         # no. of breakable edges in this element
#        nndfl=2*nedge
#        nnode=nndrl+nndfl
#        eltuel=338
#    else:
#        print 'error: fnm 3D elem type not supported for use!'
#
## write user element definition in abaqus fnm input file
#fnminp.write('*USER ELEMENT, TYPE=U'+str(eltuel)+', NODES='+str(nnode)+', COORDINATES='+str(ndim)+
#    ', PROPERTIES='+str(nprops)+', VARIABLES='+str(nsvars)+'\n')
#if ndim==2:
#    #fnminp.write('1,2\n')
#    print 'error: ndim = 2 not supported!' 
#else:
#    fnminp.write('1,2,3\n')
#  
## allocate input_elems in fnm input_elems module
#input_elems.write('        n'+elt+'='+str(elcount[eli])+'     \n')
#input_elems.write('        allocate(lib_'+elt+'(n'+elt+'))   \n')
#
## write elem nodal connec in abaqus fnm input file
#fnminp.write('*Element, type=U'+str(eltuel)+', elset=uel_'+elt+'\n')
#
## find start and end elem indices for elements of the same type
#if eli == 0:
#    elstart=0
#else:
#    elstart=sum(elcount[:eli])
#    
#elend=elstart+elcount[eli]
#
#
## write elem nodal connec of elements of the same type
#for cntr0, el in enumerate(elems[elstart:elend]):
#    cntr = cntr0+1    # index of this elem in its own type library lib_elt
#    
#    # elem line to be printed in abaqus fnm input file
#    eline = str(el.index)+','
#    
#    # elem node and edge cnc to be printed in fnm input_elems module (in string format)
#    nodecnc = ''
#    edgecnc = ''
#    
#    # elem ply node and edge cnc stored in lists
#    plynodelist = []
#    plyedgelist = []
#    
#    # include real nodes in lists
#    for k in el.nodes:
#        eline   = eline+str(k)+','
#        nodecnc = nodecnc+str(k)+','
#        plynodelist.append(k)
#        
#    # include edge fl. nodes for x-version elements
#    if (elt[0] == 'x'):    
#        for m in el.edges:
#            # update edge cnc and plyedgelist
#            edgecnc = edgecnc+str(abs(m))+','
#            plyedgelist.append(abs(m))
#            # update edge fl. nodes in eline, nodecnc and plynodelist
#            if m>0:     # elem edge's endnodes are in the same order as the definition of the edge
#                eline   = eline+str(edges[m-1].nodes[2])+','+str(edges[m-1].nodes[3])+','
#                nodecnc = nodecnc+str(edges[m-1].nodes[2])+','+str(edges[m-1].nodes[3])+','
#                plynodelist.append(edges[m-1].nodes[2])
#                plynodelist.append(edges[m-1].nodes[3])
#            elif m<0:   # elem edge's endnodes are in the reverse order as the definition of the edge
#                eline   = eline+str(edges[-m-1].nodes[3])+','+str(edges[-m-1].nodes[2])+','
#                nodecnc = nodecnc+str(edges[-m-1].nodes[3])+','+str(edges[-m-1].nodes[2])+','
#                plynodelist.append(edges[-m-1].nodes[3])
#                plynodelist.append(edges[-m-1].nodes[2])
#            else:
#                print 'warning: edges in elem ',cntr,' are empty'
#    
#    # duplicate elem nodes and edges according to blklayup
#    for ip in range(1,nplyblk):
#        eline = eline+'\n'    # print in a newline to avoid line being too long
#        for nd in plynodelist:
#            # add nodes for the ip_th ply into eline and nodecnc strings
#            eline   = eline+str(nd+ip*plynnd)+','
#            nodecnc = nodecnc+str(nd+ip*plynnd)+','
#        for eg in plyedgelist:
#            # add edges for the ip_th ply into edgecnc string
#            edgecnc = edgecnc+str(eg+ip*plynedge)+','
#        
#
#    
#    # wrap up eline (omit the last comma and append \n) and print in abaqus fnm input file
#    fnminp.write(eline[:-1]+'\n')
#    
#    
#    # write node cnc and edge cnc (for x version elems only) to the respective type arrays in fnm input_elems module
#    if (elt == 'xlam'):
#        # write elem type and keys (key to its type library lib_xlam is cntr; key to its glb library input_elems is el.index)
#        input_elems.write('        call prepare(lib_'+elt+'('+str(cntr)+'),key='+str(el.index)+', & \n')
#        # write elem associated material keys (three matkeys are needed, for ply, matrix crack and delamination respectively)
#        input_elems.write('& bulkmat='+str(ply_mkey)+', cohmat='+str(mcrack_mkey)+', interfmat='+str(delam_mkey)+', & \n')    # update matkeys later accord. to mat section assignment
#        # write elem nodal connectivity
#        if len(nodecnc) <= 100:
#            input_elems.write('& nodecnc=['+nodecnc[:-1]+'], & \n')
#        else:
#            iclist=[0]
#            ic0=-1
#            for ic, c in enumerate(nodecnc):
#                ic0=ic0+1
#                if ic0 >= 80 and c==',':
#                    iclist.append(ic)
#                    ic0=-1
#            input_elems.write('& nodecnc=[ & \n')
#            for i, ic in enumerate(iclist[:-1]):
#                istart=ic
#                iend=iclist[i+1]
#                input_elems.write('& '+nodecnc[istart:iend]+' & \n')
#            input_elems.write('& '+nodecnc[iend:-1]+'], & \n')
#            
#        # write elem edge connectivity
#        if len(edgecnc) <= 100:
#            input_elems.write('& edgecnc=['+edgecnc[:-1]+'], & \n')
#        else:
#            iclist=[0]
#            ic0=-1
#            for ic, c in enumerate(edgecnc):
#                ic0=ic0+1
#                if ic0 >= 80 and c==',':
#                    iclist.append(ic)
#                    ic0=-1
#            input_elems.write('& edgecnc=[ & \n')
#            for i, ic in enumerate(iclist[:-1]):
#                istart=ic
#                iend=iclist[i+1]
#                input_elems.write('& '+edgecnc[istart:iend]+' & \n')
#            input_elems.write('& '+edgecnc[iend:-1]+'], & \n')
#            
#        # write elem layup
#        input_elems.write('& layup=reshape(['+layupstr[:-1]+'],[2,'+str(nplyblk)+']) )\n') # need to write layup into strings of corresponding format
#        
#    else:
#        print 'unsupported fnm elem type:', elt, 'for preprocessing!'
#        print 'currently supported fnm elem types: xlam'
#        #input_elems.write('        call prepare(lib_'+elt+'('+str(cntr)+'),key='+str(el.index)+', & \n')
#        #input_elems.write('& connec=['+nodecnc[:-1]+'], & \n')
#        #input_elems.write('& matkey=1 ) \n')   # update matkey later accord. to mat section assignment
#    
#    # write elem type and typekey (elem index in its own type array) 
#    input_elems.write('        call update(input_elems('+str(el.index)+'),elname="'+elt+'",eltype="'+elt+'",typekey='+str(cntr)+') \n')
#    input_elems.write('\n')
#
## print the mandatory uel property line (not needed for calculation)
#fnminp.write('*UEL PROPERTY, elset=uel_'+elt+'\n')
#fnminp.write('1\n')
#
#
##***************************************************************        
##       Write the rest
##***************************************************************
#
#for line in lines[ln:]:
#    fnminp.write(line)





#   close input_nodes_module.f90
input_nodes.write('end subroutine input_nodes\n')
input_nodes.close()
#   close input_nodes_module.f90
input_edges.write('end subroutine input_edges\n')
input_edges.close()
#   close input_elems_module.f90
input_elems.write('end subroutine input_elems\n')
input_elems.close()



#*************************************************************** 
# copy fnm input file to main directory
#*************************************************************** 

# get current working directory
cwd = os.getcwd()
# parent directory of cwd
pwd = os.path.dirname(cwd)
# copy fnm input file to parent directory of preprocessing directory (which is assumed to be the working directory)
shutil.copy (fnminputfile,pwd)
