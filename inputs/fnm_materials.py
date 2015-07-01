################################################################
####### Preprocessing for FNM input_materials module ###########
################################################################

from fnm_classes import lamina_modulus, lamina_strength, lamina_fibreToughness, lamina,\
                        cohesive_modulus, cohesive_strength, cohesive_toughness, cohesive
import math

# -----------------------------------------------------------
# update material definitions below for the current analysis
# -----------------------------------------------------------

# IM7/8552 UD carbon epoxy lamina material properties
UDSinglePly = lamina( lamina_modulus(\
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
                      lamina_fibreToughness(\
                                        GfcT=100., 
                                        GfcC=100.   
                                           )  )

# IM7/8552 matrix crack cohesive material properties
matrixCrack = cohesive(   cohesive_modulus(\
                                        Dnn=1000000., 
                                        Dtt=1000000., 
                                        Dll=1000000.
                                          ),
                          cohesive_strength(\
                                        tau_nc=60., 
                                        tau_tc=90., 
                                        tau_lc=90.
                                           ),
                          cohesive_toughness(\
                                        Gnc=0.293, 
                                        Gtc=0.631, 
                                        Glc=0.631, 
                                        alpha=1.
                                            )    )


# IM7/8552 delamination cohesive material properties
interface = cohesive(   cohesive_modulus(\
                                        Dnn=1000000., 
                                        Dtt=1000000., 
                                        Dll=1000000.
                                     ),
                    cohesive_strength(\
                                        tau_nc=60., 
                                        tau_tc=90., 
                                        tau_lc=90.
                                      ),
                    cohesive_toughness(\
                                        Gnc=0.293, 
                                        Gtc=0.631, 
                                        Glc=0.631, 
                                        alpha=1.
                                       )    )


#***************************************************************
#   Open Fortran modules to be written during pre-processing
#***************************************************************
input_materials = open('input_materials.f90','w')


#***************************************************************
#       write input_materials_module.f90
#***************************************************************
input_materials.write('subroutine input_materials()         \n')
input_materials.write('use material_list_module,     only: UDSinglePly_material,&\n')
input_materials.write('                                  & matrixCrack_material,&\n')
input_materials.write('                                  & interface_material\n')
input_materials.write('use lamina_material_module,   only: set, lamina_modulus, &\n')
input_materials.write('                                  & lamina_strength, lamina_fibreToughness\n')
input_materials.write('use cohesive_material_module, only: set, cohesive_modulus, &\n')
input_materials.write('                                  & cohesive_strength, cohesive_toughness\n')
input_materials.write('                                   \n')
input_materials.write('  allocate(UDSinglePly_material)   \n')
input_materials.write('  allocate(matrixCrack_material)   \n')
input_materials.write('  allocate(interface_material)     \n')
input_materials.write('\n')

# write UDSinglePly material properties
input_materials.write('  call set(UDSinglePly_material, & \n')
input_materials.write('  & lamina_modulus(& \n') 
input_materials.write('  & E1   ='+str(UDSinglePly.modulus.E1)+          '_dp,& \n')
input_materials.write('  & E2   ='+str(UDSinglePly.modulus.E2)+          '_dp,& \n')
input_materials.write('  & G12  ='+str(UDSinglePly.modulus.G12)+         '_dp,& \n')
input_materials.write('  & G23  ='+str(UDSinglePly.modulus.G23)+         '_dp,& \n')
input_materials.write('  & nu12 ='+str(UDSinglePly.modulus.nu12)+        '_dp,& \n')
input_materials.write('  & nu23 ='+str(UDSinglePly.modulus.nu23)+        '_dp),& \n')
input_materials.write('  & lamina_strength(& \n')
input_materials.write('  & Xt   ='+str(UDSinglePly.strength.Xt)+         '_dp,& \n')
input_materials.write('  & Xc   ='+str(UDSinglePly.strength.Xc)+         '_dp,& \n')
input_materials.write('  & Yt   ='+str(UDSinglePly.strength.Yt)+         '_dp,& \n')
input_materials.write('  & Yc   ='+str(UDSinglePly.strength.Yc)+         '_dp,& \n')
input_materials.write('  & Sl   ='+str(UDSinglePly.strength.Sl)+         '_dp,& \n')
input_materials.write('  & St   ='+str(UDSinglePly.strength.St)+         '_dp),& \n')
input_materials.write('  & lamina_fibreToughness(& \n')  
input_materials.write('  & GfcT ='+str(UDSinglePly.fibreToughness.GfcT)+ '_dp,& \n')
input_materials.write('  & GfcC ='+str(UDSinglePly.fibreToughness.GfcC)+ '_dp)) \n')
input_materials.write('\n')


# write matrix crack material properties
input_materials.write('  call set(matrixCrack_material, & \n')
input_materials.write('  & cohesive_modulus(& \n')
input_materials.write('  & Dnn    ='+str(matrixCrack.modulus.Dnn)+     '_dp,& \n')
input_materials.write('  & Dtt    ='+str(matrixCrack.modulus.Dtt)+     '_dp,& \n')
input_materials.write('  & Dll    ='+str(matrixCrack.modulus.Dll)+     '_dp),& \n')
input_materials.write('  & cohesive_strength(& \n')
input_materials.write('  & tau_nc ='+str(matrixCrack.strength.tau_nc)+ '_dp,& \n')
input_materials.write('  & tau_tc ='+str(matrixCrack.strength.tau_tc)+ '_dp,& \n')
input_materials.write('  & tau_lc ='+str(matrixCrack.strength.tau_lc)+ '_dp),& \n')
input_materials.write('  & cohesive_toughness(& \n')
input_materials.write('  & Gnc    ='+str(matrixCrack.toughness.Gnc)+   '_dp,& \n')
input_materials.write('  & Gtc    ='+str(matrixCrack.toughness.Gtc)+   '_dp,& \n')
input_materials.write('  & Glc    ='+str(matrixCrack.toughness.Glc)+   '_dp,& \n')
input_materials.write('  & alpha  ='+str(matrixCrack.toughness.alpha)+ '_dp)) \n')
input_materials.write('\n')


# write interface material properties
input_materials.write('  call set(interface_material, & \n')
input_materials.write('  & cohesive_modulus(& \n')
input_materials.write('  & Dnn    ='+str(interface.modulus.Dnn)+     '_dp,& \n')
input_materials.write('  & Dtt    ='+str(interface.modulus.Dtt)+     '_dp,& \n')
input_materials.write('  & Dll    ='+str(interface.modulus.Dll)+     '_dp),& \n')
input_materials.write('  & cohesive_strength(& \n')
input_materials.write('  & tau_nc ='+str(interface.strength.tau_nc)+ '_dp,& \n')
input_materials.write('  & tau_tc ='+str(interface.strength.tau_tc)+ '_dp,& \n')
input_materials.write('  & tau_lc ='+str(interface.strength.tau_lc)+ '_dp),& \n')
input_materials.write('  & cohesive_toughness(& \n')
input_materials.write('  & Gnc    ='+str(interface.toughness.Gnc)+   '_dp,& \n')
input_materials.write('  & Gtc    ='+str(interface.toughness.Gtc)+   '_dp,& \n')
input_materials.write('  & Glc    ='+str(interface.toughness.Glc)+   '_dp,& \n')
input_materials.write('  & alpha  ='+str(interface.toughness.alpha)+ '_dp)) \n')
input_materials.write('\n')



#   close input_materials_module.f90
input_materials.write('end subroutine input_materials          \n')
input_materials.close()