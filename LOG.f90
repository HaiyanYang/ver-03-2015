!******************************************************************************
!
!                       PROGRAMMING LOG OF FNM3D
!
!******************************************************************************
!
!
!
!====================== DATE: 19/04/2015 ======================================
! 
! :: Notes on ig_points ::
!
!   <lamina_ig_points>
!     used by brick/wedge element integrate procedure, stores updated x, u, stress
!     strain and iterating and converged lamina sdv values at an integration point
!
!   <cohesive_ig_points>
!     used by coh3d6/8 element integrate procedure, only stores updated
!     iterating and converged cohesive sdv values at an integration point. The x,
!     u, traction and separation values are omitted for update
!
! :: Suggestions for future development ::
!   make the output-only components allocatable (x, u, stress/traction, strain/
!   separation) -> DONE!
!
!
!
!
!
!
!====================== DATE: 18/04/2015 ======================================
! 
! START THIS LOG, NOTE DOWN DETAILS OF IMPLEMENTATIONS AND CHANGES
! 
! Summarize all the modules, the derived types and their associated procedures:
!
!
!---- MATERIAL MODULES: ----
!
! :: LAMINA_MATERIAL_MODULE ::
!
!   <auxiliary objects> 
!     lamina_modulus, lamina_strength, lamina_toughness
!
!   <main objects>
!     lamina_sdv 
!       -> meaning    : solution-dependent variables used for the 
!                      consitutive law of the lamina material
!       -> procedures : extract, display
!       -> note       : external module cannot define or modify a lamina sdv;
!                       all definition and calculation of lamina sdv is done
!                       in this module only
!     lamina_material
!       -> meaning    : define a lamina material object using the aux. objects
!       -> procedures : empty, set, display, ddsdde (with failure and intact)
!       -> note       : once the material is set, no update or extraction is 
!                       allowed. It is meant to be used only for integration
!     lamina_ig_point
!       -> meaning    : define an integration point object specific for lamina
!                       material type, which stores the information corrsponding
!                       to the lamina constitutive law (cohesive law)
!       -> procedures : empty, display, update, extract
!       -> note       : this object is used frequently for extract and update, 
!                       mainly in the element integration subroutines
!   
! :: COHESIVE_MATERIAL_MODULE ::
!
!   <auxiliary objects> 
!     cohesive_modulus, cohesive_strength, cohesive_toughness
!
!   <main objects>
!     cohesive_sdv 
!       -> meaning    : solution-dependent variables used for the 
!                      consitutive law of the cohesive material
!       -> procedures : extract, display
!       -> note       : external module cannot define or modify a cohesive sdv;
!                       all definition and calculation of cohesive sdv is done
!                       in this module only
!     cohesive_material
!       -> meaning    : define a cohesive material object using the aux. objects
!       -> procedures : empty, set, display, ddsdde (with failure and intact)
!       -> note       : once the material is set, no update or extraction is 
!                       allowed. It is meant to be used only for integration
!     cohesive_ig_point
!       -> meaning    : define an integration point object specific for lamina
!                       material type, which stores the information corrsponding
!                       to the cohesive constitutive law (cohesive law)
!       -> procedures : empty, display, update, extract
!       -> note       : this object is used frequently for extract and update, 
!                       mainly in the element integration subroutines
!
!-------------------------------------------------------------------------------
!
!
!---- ELEMENT MODULES: ----
!
! :: BRICK/WEDGE_ELEMENT_MODULE ::
!   <main objects>
!     brick/wedge_element
!       -> meaning    : as name suggests
!       -> procedures : empty, set, integrate, extract
!       -> note       : once the element components are set, it cannot be 
!                       updated except through the integrate procedure
!
! :: COH3D6/COH3D8_ELEMENT_MODULE ::
!   <main objects>
!     coh3d6/coh3d8_element
!       -> meaning    : as name suggests
!       -> procedures : empty, set, integrate, extract
!       -> note       : once the element components are set, it cannot be 
!                       updated except through the integrate procedure
!
! :: basePLY_ELEMENT_MODULE ::
!   <main objects>
!     basePly_element
!       -> meaning    : an abstract 'ply' element, interfacing with wedge/brick
!                       element for actual procedures
!       -> procedures : empty, set, integrate, extract
!       -> note       : it has all the procedures of wedge/brick, and depending
!                       on a eltype parameter, it selects the procedure of the
!                       wedge or the brick element for actual calculations
!
! :: baseCOH_ELEMENT_MODULE ::
!   <main objects>
!     baseCoh_element
!       -> meaning    : an abstract cohesive element, interfacing with
!                       coh3d6/coh3d8 element for actual procedures
!       -> procedures : empty, set, integrate, extract
!       -> note       : it has all the procedures of coh3d6/coh3d8, and depending
!                       on a eltype parameter, it selects the procedure of the
!                       coh3d6 or the coh3d8 element for actual calculations
!
!-------------------------------------------------------------------------------
!
!
!---- GLOBAL MODULES: ----
!
! :: XNODE_MODULE ::
!   <main objects>
!     xnode
!       -> meaning    : an extended node object, with x, u, du, v, a, dof, ddof
!                       allocatable components for all possible analysis
!       -> procedures : empty, update, extract, display, +, -, *
!       -> note       : no set procedure because such a node cannot be set in
!                       one go. It is updated gradually during analysis (e.g., u)
!                       + and - take two xnodes to +/-, they're elemental func.
!                       * takes one xnode and one real number to multiply all
!                       components of xnode by this number
!
! :: GLOBAL_NODE_LIST_MODULE ::
!    no object defined. This module stores and saves a list of all the xnodes
!
! :: GLOBAL_MATERIAL_LIST_MODULE ::
!    no object defined. This module stores a list of all the materials defined
!
! :: GLOBAL_TOOLKIT_MODULE ::
!    no object defined. This module stores a set of tools for calculations in 
!    other modules, such as the integrate procedures of the element modules.
!
!-------------------------------------------------------------------------------
