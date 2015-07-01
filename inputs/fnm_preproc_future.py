################################################################
################# Preprocessing for the FNM ####################
##  writing input_nodes, input_edges and input_elems modules           ##
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
# material definitions for this analysis (currently for laminates only)
from fnm_materials import ply_mkey,mcrack_mkey,delam_mkey
import math
# operating system module
import os
# shutil module has functions for different kinds of copying
import shutil


#***************************************************************
#   fetch Abaqus input files
#***************************************************************

jobname = raw_input('abaqus job name:')

abqinputfile = jobname+'.inp'         # abaqus input file
fnminputfile = 'fnm-'+jobname+'.inp'  # fnm uel input file



#***************************************************************
#   define model information
#***************************************************************

ndim = 3          # dimension; currently only support 3D

# non-useful variables defined to conform with abaqus formats
nprops = 1        # no. of input material properties used in uel code (min 1)
nsvars = 1        # no. of sol. dpdnt var. used and ouput by uel code to be determined (min 1)


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
#       define usefule variables
#***************************************************************
parts   = []    # list of all parts in the model
nodes   = []    # list of all nodes in the mesh
edges   = []    # list of all breakable edges in the mesh
elems   = []    # list of all elements in the mesh
ndedg   = []    # matrix of node_i-to-node_j edge no. 
eltype  = []    # list of all element types in the mesh
elcount = []    # count of elements of each type in the mesh
bcd     = []    # list of all edges with b.c.d

rawlayup = []    # layup of the laminate in its basic format (list of fibre angles of all plies)
blklayup = []    # layup of ply blocks, used in xlam element; an xlayup class is needed


#***************************************************************
#       ask for layup from user
#***************************************************************

# ask for a list of ply angles
rawlayup = input('layup of laminate is (fibre angle (int/float) of each ply, separated by comma):')

# check if it is a list and if all elements in the list are numbers
while ((not isinstance(rawlayup, (list,tuple))) or (not all(isinstance(item, (int,float)) for item in rawlayup))):
    rawlayup=input('layup of laminate is (fibre angle (int/float) of each ply, separated by comma):')      
    
# find blocked plies and update blklayup
blklayup.append( layup(angle=rawlayup[0], nplies=0) ) # initiate blklayup
for plyangle in rawlayup:
    if (plyangle == blklayup[-1].angle):
        blklayup[-1].nplies += 1
    else:
        blklayup.append( layup(angle=plyangle, nplies=1) )


# put layup info in a string for later output purpose
layupstr=''
for lp in blklayup:
    if ('.' in str(lp.angle)):
        layupstr = layupstr+str(lp.angle)+'_dp,'
    else:
        layupstr = layupstr+str(lp.angle)+'._dp,'
        
    layupstr = layupstr+str(lp.nplies)+','




#***************************************************************
#       write input_nodes_module.f90 common codes
#***************************************************************    
input_nodes.write('subroutine input_nodes()          \n')
input_nodes.write('                                  \n')
input_nodes.write('  integer :: nnode=0              \n')   
input_nodes.write('  integer :: i=0                  \n')
input_nodes.write('\n')

#***************************************************************
#       write input_edges_module.f90 common codes
#***************************************************************
input_edges.write('subroutine input_edges()          \n')
input_edges.write('                                  \n')
input_edges.write('  integer :: nedge=0              \n')
input_edges.write('  integer :: i=0                  \n')
input_edges.write('\n')


#***************************************************************
#       write input_elems_module.f90 common codes
#***************************************************************    
input_elems.write('subroutine input_elems()          \n')
input_elems.write('                                  \n') 
input_elems.write('  integer :: nelem=0, nxlam=0     \n')
input_elems.write('  integer :: i=0                  \n')
input_elems.write('\n')



#***************************************************************
#       Read input file
#       store nodes, elements, find edges; store bcd nodes
#***************************************************************
lines = abqinp.readlines()

lenlines = len(lines)

# copy header to fnm input file
for lheader, line in enumerate(lines):

    if(len(line)>=8 and line[0:8]=='*Heading'):
        fnminp.write(line)
    elif(line[0:2]=='**'):
        fnminp.write(line)
    elif(len(line)>=9 and line[0:9]=='*Preprint'):
        fnminp.write(line)
        break

# read parts and store in python lists        
for ln in range(lheader+1:lenlines):
    line=lines[ln]
    
    if(line[0:2]=='**'):
        fnminp.write(line)
        
    elif(len(line)>=5 and line[0:5]=='*Part'):
        pname=line(12:)
        parts.append(part(name=pname, nodes=[], elems=[]))
        fnminp.write(line)
        
        isnode=False
        iselem=False
        cycle=False
        
        for ip in range(ln+1:lenlines):
            
            pline=lines[ip]
            
            cycle=False
            
            if
        

#isheader=True
#isnode=False
#iselem=False
#cycle=False
#
#for ln, line in enumerate(lines):
#    cycle=False
#    
#    if(len(line)>=5 and line[0:5]=='*Node'):
#        isheader=False
#        isnode=True
#        cycle=True
#
#    elif(len(line)>=8 and line[0:8]=='*Element'):
#        isnode=False
#        iselem=True
#        cycle=True
#        
#        if ('C3D8' in line) or ('SC8' in line):
#            eltype.append('xlam')
#            elcount.append(0)
#        else:
#            print 'warning: abaqus element type in line:', ln, 'is not supported in fnm.'
#            print 'supported abaqus element types are: C3D8(R), SC8(R)'
#
#    elif(len(line)>=9 and line[0:9]=='*End Part'):
#        iselem=False
#      
#    if not cycle:
#        if isheader:
#            fnminp.write(line)
#        elif isnode:
#        # store nodes first
#            l = []
#            for t in line.split(','):
#                try:
#                    l.append(float(t))
#                except ValueError:
#                    pass
#
#            if ndim==2:
#                # nodes.append(node(x=l[1], y=l[2], z=0.0))
#                print 'error: ndim = 2 not supported!' 
#            else:
#                nodes.append(node(x=l[1], y=l[2], z=l[3]))  
#
#            # append the ndedg row to the correct length (nnode)
#            ndedg.append([]) 
#            
#            # append the ndedg col to the correct length (nnode)
#            for i in range(len(ndedg)):
#                ndedg[i]=[0]*len(ndedg)
#            
#        elif iselem:
#        # build edges, store elements' nodecnc and edgecnc
#            
#            # increase the elem count of the latest elem type
#            elcount[-1]+=1
#        
#            l = []
#            for t in line.split(','):
#                try:
#                    l.append(int(t))
#                except ValueError:
#                    pass
# 
#            id=l[0]
#            nds=l[1:]
#            
#            elems.append(element(index=id,nodes=nds,edges=[]))
#            
#            if ndim==2:
#                print 'error: ndim = 2 not supported!'             
#            else:
#            # in 3D FNM for composites, only edges parallel to shell plane are breakable
#                for j in range(2):
#                    for i in range(len(nds)/2):
#                        row=nds[i+j*len(nds)/2]-1
#                        
#                        if i==len(nds)/2-1:
#                            col=nds[j*len(nds)/2]-1
#                        else:
#                            col=nds[i+1+j*len(nds)/2]-1
#                            
#                        if ndedg[row][col]==0:
#                        # this pair of nodes hasn't been assigned to an edge
#                            fn1=len(nodes)                          # indices of 2 fl. nodes on this new edge
#                            fn2=fn1+1
#                            for nf in range(2):                     # create floating nodes on this edge
#                                nodes.append(node(0.0,0.0,0.0))
#                                
#                            edges.append(edge(nodes=[row+1,col+1,fn1+1,fn2+1]))     # create this edge
#                            ndedg[row][col]=len(edges)              # fill the new edge index in the ndedg matrix
#                            ndedg[col][row]=-(len(edges))           # nodes in rev. order makes the same edge in rev. dir.
#                            elems[-1].edges.append(ndedg[row][col]) # append this edge no. in this elem
#                        else:
#                        # this pair of nodes has been assigned to an edge
#                            elems[-1].edges.append(ndedg[row][col]) # append this edge no. in this elem            
#            
#        else:
#            break
#            #fnminp.write(line)


#***************************************************************        
#       Write Nodes
#***************************************************************
nplyblk=len(blklayup)
plynnd=len(nodes)
plynedge=len(edges)

fnminp.write('*Node\n')
input_nodes.write('        nnode='+str(nplyblk*plynnd)+'   \n')
input_nodes.write('        allocate(input_nodes(nnode))   \n')


for iply in range(nplyblk):
    
    for cntr0, nd in enumerate(nodes):
        cntr=cntr0+1+iply*plynnd
        if ndim==2:
            #fnminp.write(str(cntr)+', '+str(nd.x)+', '+str(nd.y)+'\n')
            #input_nodes.write('        call update(input_nodes('+str(cntr)+'),x=['+str(nd.x)+'_dp,'+str(nd.y)+'_dp],u=[zero,zero])\n')
            print 'error: ndim = 2 not supported!' 
        else:
            fnminp.write(str(cntr)+', '+str(nd.x)+', '+str(nd.y)+', '+str(nd.z)+'\n')
            input_nodes.write('        call update(input_nodes('+str(cntr)+'),x=['+str(nd.x)+'_dp,'+str(nd.y)+'_dp,'+str(nd.z)+'_dp],u=[zero,zero,zero])\n')


#***************************************************************        
#       Write Edges
#***************************************************************
input_edges.write('        nedge='+str(nplyblk*plynedge)+'   \n')
input_edges.write('        allocate(input_edges(nedge))   \n')
input_edges.write('        input_edges=0                  \n')

    
#***************************************************************        
#       Write Elements
#***************************************************************
# get total no. of elements in the mesh
nelem=sum(elcount)
# allocate input_elems accordingly
input_elems.write('        nelem='+str(nelem)+'            \n')
input_elems.write('        allocate(input_elems(nelem))       \n')


# determine no. of nodes and uel element integer type
# uel element integer type:
#   1st digit: ndim
#   2nd digit: 0: bulk elem; 1: coh elem; 2: x(brick/wedge) elem; 3: xlam elem
#   3rd onwards: no. of real nodes
if ndim == 2:
        print 'error: fnm 2D elem type not supported for use!'
elif ndim == 3:      
    if elt == 'xlam':
        nndrl=8*nplyblk         # no. of (real)nodes in the element. (here, 8 for linear brick elmt)
        nedge=8*nplyblk         # no. of breakable edges in this element
        nndfl=2*nedge
        nnode=nndrl+nndfl
        eltuel=338
    else:
        print 'error: fnm 3D elem type not supported for use!'

# write user element definition in abaqus fnm input file
fnminp.write('*USER ELEMENT, TYPE=U'+str(eltuel)+', NODES='+str(nnode)+', COORDINATES='+str(ndim)+
    ', PROPERTIES='+str(nprops)+', VARIABLES='+str(nsvars)+'\n')
if ndim==2:
    #fnminp.write('1,2\n')
    print 'error: ndim = 2 not supported!' 
else:
    fnminp.write('1,2,3\n')
  
# allocate input_elems in fnm input_elems module
input_elems.write('        n'+elt+'='+str(elcount[eli])+'     \n')
input_elems.write('        allocate(lib_'+elt+'(n'+elt+'))   \n')

# write elem nodal connec in abaqus fnm input file
fnminp.write('*Element, type=U'+str(eltuel)+', elset=uel_'+elt+'\n')

# find start and end elem indices for elements of the same type
if eli == 0:
    elstart=0
else:
    elstart=sum(elcount[:eli])
    
elend=elstart+elcount[eli]


# write elem nodal connec of elements of the same type
for cntr0, el in enumerate(elems[elstart:elend]):
    cntr = cntr0+1    # index of this elem in its own type library lib_elt
    
    # elem line to be printed in abaqus fnm input file
    eline = str(el.index)+','
    
    # elem node and edge cnc to be printed in fnm input_elems module (in string format)
    nodecnc = ''
    edgecnc = ''
    
    # elem ply node and edge cnc stored in lists
    plynodelist = []
    plyedgelist = []
    
    # include real nodes in lists
    for k in el.nodes:
        eline   = eline+str(k)+','
        nodecnc = nodecnc+str(k)+','
        plynodelist.append(k)
        
    # include edge fl. nodes for x-version elements
    if (elt[0] == 'x'):    
        for m in el.edges:
            # update edge cnc and plyedgelist
            edgecnc = edgecnc+str(abs(m))+','
            plyedgelist.append(abs(m))
            # update edge fl. nodes in eline, nodecnc and plynodelist
            if m>0:     # elem edge's endnodes are in the same order as the definition of the edge
                eline   = eline+str(edges[m-1].nodes[2])+','+str(edges[m-1].nodes[3])+','
                nodecnc = nodecnc+str(edges[m-1].nodes[2])+','+str(edges[m-1].nodes[3])+','
                plynodelist.append(edges[m-1].nodes[2])
                plynodelist.append(edges[m-1].nodes[3])
            elif m<0:   # elem edge's endnodes are in the reverse order as the definition of the edge
                eline   = eline+str(edges[-m-1].nodes[3])+','+str(edges[-m-1].nodes[2])+','
                nodecnc = nodecnc+str(edges[-m-1].nodes[3])+','+str(edges[-m-1].nodes[2])+','
                plynodelist.append(edges[-m-1].nodes[3])
                plynodelist.append(edges[-m-1].nodes[2])
            else:
                print 'warning: edges in elem ',cntr,' are empty'
    
    # duplicate elem nodes and edges according to blklayup
    for ip in range(1,nplyblk):
        eline = eline+'\n'    # print in a newline to avoid line being too long
        for nd in plynodelist:
            # add nodes for the ip_th ply into eline and nodecnc strings
            eline   = eline+str(nd+ip*plynnd)+','
            nodecnc = nodecnc+str(nd+ip*plynnd)+','
        for eg in plyedgelist:
            # add edges for the ip_th ply into edgecnc string
            edgecnc = edgecnc+str(eg+ip*plynedge)+','
        

    
    # wrap up eline (omit the last comma and append \n) and print in abaqus fnm input file
    fnminp.write(eline[:-1]+'\n')
    
    
    # write node cnc and edge cnc (for x version elems only) to the respective type arrays in fnm input_elems module
    if (elt == 'xlam'):
        # write elem type and keys (key to its type library lib_xlam is cntr; key to its glb library input_elems is el.index)
        input_elems.write('        call prepare(lib_'+elt+'('+str(cntr)+'),key='+str(el.index)+', & \n')
        # write elem associated material keys (three matkeys are needed, for ply, matrix crack and delamination respectively)
        input_elems.write('& bulkmat='+str(ply_mkey)+', cohmat='+str(mcrack_mkey)+', interfmat='+str(delam_mkey)+', & \n')    # update matkeys later accord. to mat section assignment
        # write elem nodal connectivity
        if len(nodecnc) <= 100:
            input_elems.write('& nodecnc=['+nodecnc[:-1]+'], & \n')
        else:
            iclist=[0]
            ic0=-1
            for ic, c in enumerate(nodecnc):
                ic0=ic0+1
                if ic0 >= 80 and c==',':
                    iclist.append(ic)
                    ic0=-1
            input_elems.write('& nodecnc=[ & \n')
            for i, ic in enumerate(iclist[:-1]):
                istart=ic
                iend=iclist[i+1]
                input_elems.write('& '+nodecnc[istart:iend]+' & \n')
            input_elems.write('& '+nodecnc[iend:-1]+'], & \n')
            
        # write elem edge connectivity
        if len(edgecnc) <= 100:
            input_elems.write('& edgecnc=['+edgecnc[:-1]+'], & \n')
        else:
            iclist=[0]
            ic0=-1
            for ic, c in enumerate(edgecnc):
                ic0=ic0+1
                if ic0 >= 80 and c==',':
                    iclist.append(ic)
                    ic0=-1
            input_elems.write('& edgecnc=[ & \n')
            for i, ic in enumerate(iclist[:-1]):
                istart=ic
                iend=iclist[i+1]
                input_elems.write('& '+edgecnc[istart:iend]+' & \n')
            input_elems.write('& '+edgecnc[iend:-1]+'], & \n')
            
        # write elem layup
        input_elems.write('& layup=reshape(['+layupstr[:-1]+'],[2,'+str(nplyblk)+']) )\n') # need to write layup into strings of corresponding format
        
    else:
        print 'unsupported fnm elem type:', elt, 'for preprocessing!'
        print 'currently supported fnm elem types: xlam'
        #input_elems.write('        call prepare(lib_'+elt+'('+str(cntr)+'),key='+str(el.index)+', & \n')
        #input_elems.write('& connec=['+nodecnc[:-1]+'], & \n')
        #input_elems.write('& matkey=1 ) \n')   # update matkey later accord. to mat section assignment
    
    # write elem type and typekey (elem index in its own type array) 
    input_elems.write('        call update(input_elems('+str(el.index)+'),elname="'+elt+'",eltype="'+elt+'",typekey='+str(cntr)+') \n')
    input_elems.write('\n')

# print the mandatory uel property line (not needed for calculation)
fnminp.write('*UEL PROPERTY, elset=uel_'+elt+'\n')
fnminp.write('1\n')


#***************************************************************        
#       Write the rest
#***************************************************************

for line in lines[ln:]:
    fnminp.write(line)





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
