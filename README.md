# OOFEM Enhanced features

This is an Objective Oriented Finite Element (OOFEM) package (Version 2.4) enhanced with 
-  Two nonlocal damage constitutive model for brittle material
-  Intrinsic and Extrinsic PPR cohesive model
-  Anisotropic nonlocal regularization
-  Multiscale fracture propagation by coupling nonlocal damage with cohesive law
-  Two global solution techniques including  cylindrical arc length control and local normal plane arc length control

The readers are referred to check readme-original.md file for installation or to check the website www.oofem.org for OOFEM documentations. Above enhancements are published in the following papers:

- Jin, Wencheng, and Chloé Arson. "Anisotropic nonlocal damage model for materials with intrinsic transverse isotropy." International Journal of Solids and Structures 139 (2018): 29-42.
- Jin, Wencheng, and Chloé Arson. "Nonlocal enrichment of a micromechanical damage model with tensile softening: Advantages and limitations." Computers and Geotechnics 94 (2018): 196-206.
- Jin, Wencheng, and Chloé Arson. "XFEM to couple nonlocal micromechanics damage with discrete mode I cohesive fracture." Computational Methods in Applied Mechanics and Engineering, https://doi.org/10.1016/j.cma.2019.112617.

**Please kindly cite above papers if you used any of the functions or algorithms listed in this Github repository**

## Added source code files

```
adm1.c
adm1.h
adm2.c
adm2.h
admnl1.c
admnl1.h
admnl2.c
admnl2.h
anisodamagemodel.h
idmnl1.c
intmatpprexcz.c
intmatpprexcz.h
intmatpprincz.c
intmatpprincz.h
nonlocalmaterialext.c
nonlocalmaterialext.h
plnonlocaldamage.c
plnonlocaldamage.h
plnonlocalstress.c
plnonlocalstress.h
staticstructuralarc.c
staticstructuralarc.h
structuralmaterial.c
structuralmaterial.h
```

## Added input file folders

```
1 tests/nonlocal/
2 tests/multiscal/
```

