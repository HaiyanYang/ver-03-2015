################################################################
####### Preprocessing for writing FNM lib_mat module ###########
################################################################

from fnm_classes import* # glb objects defined for FNM
import math


#***************************************************************
#   define material information
#***************************************************************

# create a library of materials
list_mat=[]         # library of all types of materials
list_iso=[]         # library of isotropic materials
list_lamina=[]      # library of lamina materials
list_interface=[]   # library of interface-type materials (cohesive, could be for matrix crack/delamination)

# indices of UD ply bulk material, matrix crack cohesive material and delamination cohesive material
# in the list_mat arrays
ply_mkey=1      # first material in the list should be for ud ply bulk (first of type lamina)
mcrack_mkey=2    # second should be for matrix cracking (first of type interface)
delam_mkey=3   # third should be for delamination (second of type interface)


# *** change/add below for other lamina/cohesive properties ***

# -----------------------------------------------------------
# update material definitions here for the current analysis
# -----------------------------------------------------------

# IM7/8552 UD carbon epoxy lamina material properties
laminate_Ply=lamina( lamina_modulus(\
                                        E1=161000., 
                                        E2=11400., 
                                        G12=5170., 
                                        G23=3980., 
                                        nu12=0.34, 
                                        nu23=0.43   
                               ),
                lamina_strength(\
                                        Xt=2806., 
                                        Xc=1400., 
                                        Yt=60., 
                                        Yc=185., 
                                        Sl=90., 
                                        St=90.     
                                ),
                lamina_matrixtoughness(\
                                        GmcI=0.2, 
                                        GmcII=1., 
                                        eta=1.
                                       ),
                lamina_fibretoughness(\
                                        GfcT=100., 
                                        GfcC=100.   
                                     )  )

# IM7/8552 matrix crack cohesive material properties
laminate_Matrix=interface(   interface_modulus(\
                                        Dnn=1000000., 
                                        Dtt=1000000., 
                                        Dll=1000000.
                                     ),
                    interface_strength(\
                                        tau_nc=60., 
                                        tau_tc=90., 
                                        tau_lc=90.
                                      ),
                    interface_toughness(\
                                        Gnc=0.293, 
                                        Gtc=0.631, 
                                        Glc=0.631, 
                                        eta=1.
                                       )    )


# IM7/8552 delamination cohesive material properties
# (assumed to be the same for matrix crack, except that strengths are reduced to artificially enlarge cohesive zone)
laminate_Delamination=interface(   interface_modulus(\
                                        Dnn=1000000., 
                                        Dtt=1000000., 
                                        Dll=1000000.
                                     ),
                    interface_strength(\
                                        tau_nc=60., 
                                        tau_tc=90., 
                                        tau_lc=90.
                                      ),
                    interface_toughness(\
                                        Gnc=0.293, 
                                        Gtc=0.631, 
                                        Glc=0.631, 
                                        eta=1.
                                       )    )


# -----------------------------------------------------------------------
# append current analysis materials in the respective material-type list
# -----------------------------------------------------------------------
list_lamina.append(laminate_Ply)
list_interface.append(laminate_Matrix)
list_interface.append(laminate_Delamination)


# -----------------------------------------------------------------------
# append current analysis materials in the global material list
# -----------------------------------------------------------------------
list_mat.append(material(name="'laminate Ply'",          type="'lamina'",      typekey=1))
list_mat.append(material(name="'laminate Matrix Crack'",   type="'interface'",   typekey=1))
list_mat.append(material(name="'laminate Delamination'",   type="'interface'",   typekey=2))



nmat=len(list_mat)             # no. of materials in the model, = sum of no. of the following individual materials
niso=len(list_iso)             # no. of isotropic materials
nlamina=len(list_lamina)       # no. of lamina materials
ninterface=len(list_interface) # no. of cohesive materials



#***************************************************************
#   Open Fortran modules to be written during pre-processing
#***************************************************************
lib_mat=open('init_lib_mat.f90','w')    # array of all materials


#***************************************************************
#       write lib_mat_module.f90
#***************************************************************
lib_mat.write('    subroutine initialize_lib_mat()                            \n')
lib_mat.write('                                                               \n')
lib_mat.write('        integer :: i=0, nmat=0, niso=0, nlamina=0, ninterface=0\n')

lib_mat.write('        nmat='+str(nmat)+'                                            \n')
lib_mat.write('        niso='+str(niso)+'                                            \n')
lib_mat.write('        nlamina='+str(nlamina)+'                                      \n')
lib_mat.write('        ninterface='+str(ninterface)+'                                \n')
lib_mat.write('        if(nmat>0) allocate(lib_mat(nmat))                       \n')
lib_mat.write('        if(niso>0) allocate(lib_iso(niso))                       \n')
lib_mat.write('        if(nlamina>0) allocate(lib_lamina(nlamina))              \n')
lib_mat.write('        if(ninterface>0) allocate(lib_interface(ninterface))     \n')

# write material section info
if(nmat>0):
    for i in range(nmat):
        lib_mat.write('        call update(lib_mat('+str(i+1)+'),matname='+list_mat[i].name+\
        ',mattype='+list_mat[i].type+',typekey='+str(list_mat[i].typekey)+')\n')


# write lamina material properties
if(nlamina>0):
    for i in range(nlamina):
        lib_mat.write('        call update(lib_lamina('+str(i+1)+'), & \n')
        lib_mat.write('      & lamina_modulus(& \n') 
        lib_mat.write('      & E1='                +str(list_lamina[i].modulus.E1)+               '_dp,& \n')
        lib_mat.write('      & E2='                +str(list_lamina[i].modulus.E2)+               '_dp,& \n')
        lib_mat.write('      & G12='                +str(list_lamina[i].modulus.G12)+              '_dp,& \n')
        lib_mat.write('      & G23='                +str(list_lamina[i].modulus.G23)+              '_dp,& \n')
        lib_mat.write('      & nu12='                +str(list_lamina[i].modulus.nu12)+             '_dp,& \n')
        lib_mat.write('      & nu23='                +str(list_lamina[i].modulus.nu23)+            '_dp),& \n')
        lib_mat.write('      & lamina_strength(& \n')
        lib_mat.write('      & Xt='                +str(list_lamina[i].strength.Xt)+             '_dp,& \n')
        lib_mat.write('      & Xc='                +str(list_lamina[i].strength.Xc)+              '_dp,& \n')
        lib_mat.write('      & Yt='                +str(list_lamina[i].strength.Yt)+              '_dp,& \n')
        lib_mat.write('      & Yc='                +str(list_lamina[i].strength.Yc)+              '_dp,& \n')
        lib_mat.write('      & Sl='                +str(list_lamina[i].strength.Sl)+              '_dp,& \n')
        lib_mat.write('      & St='                +str(list_lamina[i].strength.St)+             '_dp),& \n')
        lib_mat.write('      & lamina_fibretoughness(& \n')  
        lib_mat.write('      & GfcT='                         +str(list_lamina[i].fibretoughness.GfcT)+      '_dp,& \n')
        lib_mat.write('      & GfcC='                         +str(list_lamina[i].fibretoughness.GfcC)+      '_dp)& \n')
        lib_mat.write('      & ) \n')
        lib_mat.write('\n')



# write interface material properties
if(ninterface>0):
    for i in range(ninterface):
        lib_mat.write('        call update(lib_interface('+str(i+1)+'), & \n')
        lib_mat.write('      & interface_modulus(& \n')
        lib_mat.write('      & Dnn='                   +str(list_interface[i].modulus.Dnn)+          '_dp,& \n')
        lib_mat.write('      & Dtt='                   +str(list_interface[i].modulus.Dtt)+          '_dp,& \n')
        lib_mat.write('      & Dll='                   +str(list_interface[i].modulus.Dll)+         '_dp),& \n')
        lib_mat.write('      & interface_strength(& \n')
        lib_mat.write('      & tau_nc='                    +str(list_interface[i].strength.tau_nc)+     '_dp,& \n')
        lib_mat.write('      & tau_tc='                    +str(list_interface[i].strength.tau_tc)+     '_dp,& \n')
        lib_mat.write('      & tau_lc='                    +str(list_interface[i].strength.tau_lc)+    '_dp),& \n')
        lib_mat.write('      & interface_toughness(& \n')
        lib_mat.write('      & Gnc='                     +str(list_interface[i].toughness.Gnc)+      '_dp,& \n')
        lib_mat.write('      & Gtc='                     +str(list_interface[i].toughness.Gtc)+      '_dp,& \n')
        lib_mat.write('      & Glc='                     +str(list_interface[i].toughness.Glc)+      '_dp,& \n')
        lib_mat.write('      & eta='                     +str(list_interface[i].toughness.eta)+      '_dp)) \n')
        lib_mat.write('\n')



#   close lib_mat_module.f90
lib_mat.write('    end subroutine initialize_lib_mat          \n')
lib_mat.close()