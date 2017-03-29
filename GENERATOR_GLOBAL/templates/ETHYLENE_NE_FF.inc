    &NONBONDED
     &LENNARD-JONES
      ATOMS                     CP    CP
      EPSILON   [kcalmol]       0.08600D0
      SIGMA     [angstrom]      3.399670D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     CP   H
      EPSILON   [kcalmol]       0.03675D0
      SIGMA     [angstrom]      3.02460D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     CP  CN
      EPSILON   [kcalmol]       0.08600D0
      SIGMA     [angstrom]      3.399670D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     CP   Ne
      EPSILON   [kcalmol]       0.07821D0
      SIGMA     [angstrom]      3.09084D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     H   H
      EPSILON   [kcalmol]       0.01570D0
      SIGMA     [angstrom]      2.64953D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     H   CN
      EPSILON   [kcalmol]       0.03655D0
      SIGMA     [angstrom]      3.02460D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     H   Ne
      EPSILON   [kcalmol]       0.03342D0
      SIGMA     [angstrom]      2.71577D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     CN  CN
      EPSILON   [kcalmol]       0.08600D0
      SIGMA     [angstrom]      3.399670D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     CN  Ne
      EPSILON   [kcalmol]       0.07821D0
      SIGMA     [angstrom]      3.09084D0
     &END LENNARD-JONES
     &LENNARD-JONES
      ATOMS                     Ne  Ne
      EPSILON   [kcalmol]       0.07112D0
      SIGMA     [angstrom]      2.782D0
     &END LENNARD-JONES
    &END NONBONDED 


#
    &BOND
     ATOMS                      CP CP
     KIND                       AMBER
#     K   [kcalmol*angstrom^-2]  589.70D0
#UNITS hartree/bohr^2
     K                          0.263099D0
#     R0         [bohr]          2.6475D0
     R0         [angstrom]       sedCC_CHARGED
    &END BOND
    &BOND
     ATOMS                      CP H
     KIND                       AMBER
#     K   [kcalmol*angstrom^-2]  344.30D0
#UNITS hartree/bohr^2
K                          0.15361D0
#     R0         [bohr]          2.0540D0
     R0         [angstrom]       1.0868D0
    &END BOND
    &BOND
     ATOMS                      CN CN
     KIND                       AMBER
#     K   [kcalmol*angstrom^-2]  589.70D0
#UNITS hartree/bohr^2
     K                          0.263099D0
#     R0         [bohr]          2.5020D0
     R0         [angstrom]       1.3239D0
    &END BOND
    &BOND
     ATOMS                      CN H
     KIND                       AMBER
#     K   [kcalmol*angstrom^-2]  344.30D0
#UNITS hartree/bohr^2
K                          0.15361D0
#     R0         [bohr]          2.0540D0
     R0         [angstrom]       1.0868D0
    &END BOND
#


    &BEND
     ATOMS                      H  CN  H
     KIND                       AMBER
#     K                          37.810D0
#UNITS hartree/rad^2
     K                          0.065026914D0
     THETA0     [rad]           2.075546440D0
    &END BEND
    &BEND
     ATOMS                      H  CN  CN
     KIND                       AMBER
#     K                          50.30D0
#UNITS hartree/rad^2
     K                          0.0801782D0
     THETA0     [rad]           2.08916001D0
    &END BEND
    &BEND
     ATOMS                      H  CP  H
     KIND                       AMBER
#     K                          37.810D0
#UNITS hartree/rad^2
     K                          0.065026914D0
     THETA0     [rad]           2.075546440D0
    &END BEND
    &BEND
     ATOMS                      H  CP  CP
     KIND                       AMBER
#     K                          50.30D0
#UNITS hartree/rad^2
     K                          0.0801782D0
     THETA0     [rad]           2.08916001D0
    &END BEND

#
    &TORSION
     ATOMS                      H CN CN H
     KIND                       AMBER
     K       [kcalmol]          6.650D0
     M                          2
     PHI0    [rad]              3.1415940D0
    &END TORSION
    &TORSION
     ATOMS                      H CP CP H
     KIND                       AMBER
     K       [kcalmol]          6.650D0
     M                          2
     PHI0    [rad]              3.1415940D0
    &END TORSION

