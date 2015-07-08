subroutine set_fnm_materials()     
use parameter_module,         only: DP 
use material_list_module,     only: UDSinglePly_material,&
                                  & matrixCrack_material,&
                                  & interface_material
use lamina_material_module,   only: set, lamina_modulus, &
                                  & lamina_strength, lamina_fibreToughness
use cohesive_material_module, only: set, cohesive_modulus, &
                                  & cohesive_strength, cohesive_toughness
                                   
  allocate(UDSinglePly_material)   
  allocate(matrixCrack_material)   
  allocate(interface_material)     

  call set(UDSinglePly_material, & 
  & lamina_modulus(& 
  & E1   =161000.0_dp,& 
  & E2   =11400.0_dp,& 
  & G12  =5170.0_dp,& 
  & G23  =3980.0_dp,& 
  & nu12 =0.0_dp,& 
  & nu23 =0.0_dp),& 
  & lamina_strength(& 
  & Xt   =2806.0_dp,& 
  & Xc   =1400.0_dp,& 
  & Yt   =60.0_dp,& 
  & Yc   =185.0_dp,& 
  & Sl   =90.0_dp,& 
  & St   =90.0_dp),& 
  & lamina_fibreToughness(& 
  & GfcT =100.0_dp,& 
  & GfcC =100.0_dp)) 

  call set(matrixCrack_material, & 
  & cohesive_modulus(& 
  & Dnn    =1000000.0_dp,& 
  & Dtt    =1000000.0_dp,& 
  & Dll    =1000000.0_dp),& 
  & cohesive_strength(& 
  & tau_nc =60.0_dp,& 
  & tau_tc =90.0_dp,& 
  & tau_lc =90.0_dp),& 
  & cohesive_toughness(& 
  & Gnc    =0.293_dp,& 
  & Gtc    =0.631_dp,& 
  & Glc    =0.631_dp,& 
  & alpha  =1.0_dp)) 

  call set(interface_material, & 
  & cohesive_modulus(& 
  & Dnn    =1000000.0_dp,& 
  & Dtt    =1000000.0_dp,& 
  & Dll    =1000000.0_dp),& 
  & cohesive_strength(& 
  & tau_nc =60.0_dp,& 
  & tau_tc =90.0_dp,& 
  & tau_lc =90.0_dp),& 
  & cohesive_toughness(& 
  & Gnc    =0.293_dp,& 
  & Gtc    =0.631_dp,& 
  & Glc    =0.631_dp,& 
  & alpha  =1.0_dp)) 

end subroutine set_fnm_materials          
