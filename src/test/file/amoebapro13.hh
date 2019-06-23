namespace tinker_6fe8e913fe4da3d46849d10248ad2a4872b4da93 {
const char* amoebapro13_prm =
R"**(

      ##############################
      ##                          ##
      ##  Force Field Definition  ##
      ##                          ##
      ##############################


forcefield              AMOEBA-PROTEIN-2013

bond-cubic              -2.55
bond-quartic            3.793125
angle-trigonal          IN-PLANE
angle-cubic             -0.014
angle-quartic           0.000056
angle-pentic            -0.0000007
angle-sextic            0.000000022
opbendtype              ALLINGER
opbend-cubic            -0.014
opbend-quartic          0.000056
opbend-pentic           -0.0000007
opbend-sextic           0.000000022
torsionunit             0.5
vdwtype                 BUFFERED-14-7
radiusrule              CUBIC-MEAN
radiustype              R-MIN
radiussize              DIAMETER
epsilonrule             HHG
dielectric              1.0
polarization            MUTUAL
vdw-12-scale            0.0
vdw-13-scale            0.0
vdw-14-scale            1.0
vdw-15-scale            1.0
mpole-12-scale          0.0
mpole-13-scale          0.0
mpole-14-scale          0.4
mpole-15-scale          0.8
polar-12-scale          0.0
polar-13-scale          0.0
polar-14-scale          1.0
polar-15-scale          1.0
polar-12-intra          0.0
polar-13-intra          0.0
polar-14-intra          0.5
polar-15-intra          1.0
direct-11-scale         0.0
direct-12-scale         1.0
direct-13-scale         1.0
direct-14-scale         1.0
mutual-11-scale         1.0
mutual-12-scale         1.0
mutual-13-scale         1.0
mutual-14-scale         1.0


      #############################
      ##                         ##
      ##  Literature References  ##
      ##                         ##
      #############################


This file contains the 2013 version of the AMOEBA protein force field as
per the J. Chem. Theory Comput. article listed immediately below

Y. Shi, Z. Xia, J. Zhang, R. Best, J. W. Ponder and P. Ren, "Polarizable
Atomic Multipole-Based AMOEBA Force Field for Proteins", J. Chem. Theory
Comput., 9, 4046-4063 (2013)

P. Ren, C. Wu and J. W. Ponder, "Polarizable Atomic Multipole-Based
Molecular Mechanics for Organic Molecules", J. Chem. Theory Comput.,
7, 3143-3161 (2011)

J. C. Wu, J.-P. Piquemal, R. Chaudret, P. Reinhardt and P. Ren,
"Polarizable Molecular Dynamics Simulation of Zn(II) in Water Using
the AMOEBA Force Field", J. Chem. Theory Comput., 6, 2059-2070 (2010)

A. Grossfield, P. Ren, J. W. Ponder, "Ion Solvation Thermodynamics from
Simulation with a Polarizable Force Field", J. Am. Chem. Soc., 125,
15671-15682 (2003)

P. Ren and J. W. Ponder, "Polarizable Atomic Multipole Water Model for
Molecular Mechanics Simulation", J. Phys. Chem. B, 107, 5933-5947 (2003)

Monovalent ion parameters taken from Zhi Wang, Ph.D. thesis, Department
of Chemistry, Washington University in St. Louis, May 2018; available
from https://dasher.wustl.edu/ponder/


   ###############################################
   ##                                           ##
   ##  AMOEBA Protein Force Field Atom Classes  ##
   ##                                           ##
   ##   1  Backbone Amide Nitrogen              ##
   ##   2  Glycine Alpha Carbon                 ##
   ##   3  Backbone Carbonyl Carbon             ##
   ##   4  Amide or Guanidinium Hydrogen        ##
   ##   5  Amide Carbonyl Oxygen                ##
   ##   6  Methine Hydrogen                     ##
   ##   7  Methine Carbon                       ##
   ##   8  Methyl or Methylene Carbon           ##
   ##   9  Methyl or Methylene Hydrogen         ##
   ##  10  Hydroxyl Oxygen                      ##
   ##  11  Hydroxyl Hydrogen                    ##
   ##  12  Sulfide or Disulfide Sulfur          ##
   ##  13  Sulfhydryl Hydrogen                  ##
   ##  14  Thiolate Sulfur                      ##
   ##  15  Proline Backbone Nitrogen            ##
   ##  16  Proline Ring Methylene Carbon        ##
   ##  17  Phenyl Carbon                        ##
   ##  18  Phenyl Hydrogen                      ##
   ##  19  Phenolic Oxygen                      ##
   ##  20  Phenolic Hydrogen                    ##
   ##  21  Phenoxide Oxygen                     ##
   ##  22  Indole Carbon                        ##
   ##  23  Indole CH Hydrogen                   ##
   ##  24  Imidazole or Indole NH Nitrogen      ##
   ##  25  Imidazole or Indole NH Hydrogen      ##
   ##  26  Imidazole C=C Carbon                 ##
   ##  27  Imidazole CH Hydrogen                ##
   ##  28  Imidazole N=C-N Carbon               ##
   ##  29  Imidazole C=N Nitrogen               ##
   ##  30  Carboxylate Carbon                   ##
   ##  31  Carboxylate Oxygen                   ##
   ##  32  Carboxylic Acid Carbonyl Carbon      ##
   ##  33  Carboxylic Acid Carbonyl Oxygen      ##
   ##  34  Carboxylic Acid Hydroxyl Oxygen      ##
   ##  35  Carboxylic Acid Hydrogen             ##
   ##  36  Lysine/Ornithine Gamma Carbon        ##
   ##  37  Ammonium Nitrogen                    ##
   ##  38  Ammonium Hydrogen                    ##
   ##  39  Guanidinium Carbon                   ##
   ##  40  Acetyl or NMe Methyl Carbon          ##
   ##  41  N-Terminal Ammonium Nitrogen         ##
   ##  42  N-Terminal Ammonium Hydrogen         ##
   ##                                           ##
   ###############################################


      #############################
      ##                         ##
      ##  Atom Type Definitions  ##
      ##                         ##
      #############################


atom          1    1    N     "Glycine N"                    7    14.007    3
atom          2    2    CA    "Glycine CA"                   6    12.011    4
atom          3    3    C     "Glycine C"                    6    12.011    3
atom          4    4    HN    "Glycine HN"                   1     1.008    1
atom          5    5    O     "Glycine O"                    8    15.999    1
atom          6    6    H     "Glycine HA"                   1     1.008    1
atom          7    1    N     "Alanine N"                    7    14.007    3
atom          8    7    CA    "Alanine CA"                   6    12.011    4
atom          9    3    C     "Alanine C"                    6    12.011    3
atom         10    4    HN    "Alanine HN"                   1     1.008    1
atom         11    5    O     "Alanine O"                    8    15.999    1
atom         12    6    H     "Alanine HA"                   1     1.008    1
atom         13    8    C     "Alanine CB"                   6    12.011    4
atom         14    9    H     "Alanine HB"                   1     1.008    1
atom         15    7    C     "Valine CB"                    6    12.011    4
atom         16    6    H     "Valine HB"                    1     1.008    1
atom         17    8    C     "Valine CG"                    6    12.011    4
atom         18    9    H     "Valine HG"                    1     1.008    1
atom         19    8    C     "Leucine CB"                   6    12.011    4
atom         20    9    H     "Leucine HB"                   1     1.008    1
atom         21    7    C     "Leucine CG"                   6    12.011    4
atom         22    6    H     "Leucine HG"                   1     1.008    1
atom         23    8    C     "Leucine CD"                   6    12.011    4
atom         24    9    H     "Leucine HD"                   1     1.008    1
atom         25    7    C     "Isoleucine CB"                6    12.011    4
atom         26    6    H     "Isoleucine HB"                1     1.008    1
atom         27    8    C     "Isoleucine CG1"               6    12.011    4
atom         28    9    H     "Isoleucine HG1"               1     1.008    1
atom         29    8    C     "Isoleucine CG2"               6    12.011    4
atom         30    9    H     "Isoleucine HG2"               1     1.008    1
atom         31    8    C     "Isoleucine CD"                6    12.011    4
atom         32    9    H     "Isoleucine HD"                1     1.008    1
atom         33    8    C     "Serine CB"                    6    12.011    4
atom         34    9    H     "Serine HB"                    1     1.008    1
atom         35   10    OH    "Serine OG"                    8    15.999    2
atom         36   11    HO    "Serine HG"                    1     1.008    1
atom         37    7    C     "Threonine CB"                 6    12.011    4
atom         38    6    H     "Threonine HB"                 1     1.008    1
atom         39   10    OH    "Threonine OG1"                8    15.999    2
atom         40   11    HO    "Threonine HG1"                1     1.008    1
atom         41    8    C     "Threonine CG2"                6    12.011    4
atom         42    9    H     "Threonine HG2"                1     1.008    1
atom         43    8    C     "Cysteine CB"                  6    12.011    4
atom         44    9    H     "Cysteine HB"                  1     1.008    1
atom         45   12    SH    "Cysteine SG"                 16    32.066    2
atom         46   13    HS    "Cysteine HG"                  1     1.008    1
atom         47   12    SS    "Cystine SG"                  16    32.066    2
atom         48    7    CA    "Cysteine Anion CA"            6    12.011    4
atom         49   14    S     "Cysteine Anion S-"           16    32.066    1
atom         50   15    N     "Proline N"                    7    14.007    3
atom         51    7    CA    "Proline CA"                   6    12.011    4
atom         52    3    C     "Proline C"                    6    12.011    3
atom         53    5    O     "Proline O"                    8    15.999    1
atom         54    6    H     "Proline HA"                   1     1.008    1
atom         55   16    C     "Proline CB"                   6    12.011    4
atom         56    9    H     "Proline HB"                   1     1.008    1
atom         57   16    C     "Proline CG"                   6    12.011    4
atom         58    9    H     "Proline HG"                   1     1.008    1
atom         59   16    C     "Proline CD"                   6    12.011    4
atom         60    6    H     "Proline HD"                   1     1.008    1
atom         61    8    C     "Phenylalanine CB"             6    12.011    4
atom         62    9    H     "Phenylalanine HB"             1     1.008    1
atom         63   17    C     "Phenylalanine CG"             6    12.011    3
atom         64   17    C     "Phenylalanine CD"             6    12.011    3
atom         65   18    H     "Phenylalanine HD"             1     1.008    1
atom         66   17    C     "Phenylalanine CE"             6    12.011    3
atom         67   18    H     "Phenylalanine HE"             1     1.008    1
atom         68   17    C     "Phenylalanine CZ"             6    12.011    3
atom         69   18    H     "Phenylalanine HZ"             1     1.008    1
atom         70    8    C     "Tyrosine CB"                  6    12.011    4
atom         71    9    H     "Tyrosine HB"                  1     1.008    1
atom         72   17    C     "Tyrosine CG"                  6    12.011    3
atom         73   17    C     "Tyrosine CD"                  6    12.011    3
atom         74   18    H     "Tyrosine HD"                  1     1.008    1
atom         75   17    C     "Tyrosine CE"                  6    12.011    3
atom         76   18    H     "Tyrosine HE"                  1     1.008    1
atom         77   17    C     "Tyrosine CZ"                  6    12.011    3
atom         78   19    OH    "Tyrosine OH"                  8    15.999    2
atom         79   20    HO    "Tyrosine HH"                  1     1.008    1
atom         80    8    C     "Tyrosine Anion CB"            6    12.011    4
atom         81    9    H     "Tyrosine Anion HB"            1     1.008    1
atom         82   17    C     "Tyrosine Anion CG"            6    12.011    3
atom         83   17    C     "Tyrosine Anion CD"            6    12.011    3
atom         84   18    H     "Tyrosine Anion HD"            1     1.008    1
atom         85   17    C     "Tyrosine Anion CE"            6    12.011    3
atom         86   18    H     "Tyrosine Anion HE"            1     1.008    1
atom         87   17    C     "Tyrosine Anion CZ"            6    12.011    3
atom         88   21    O-    "Tyrosine Anion O-"            8    15.999    1
atom         89    8    C     "Tryptophan CB"                6    12.011    4
atom         90    9    H     "Tryptophan HB"                1     1.008    1
atom         91   22    C     "Tryptophan CG"                6    12.011    3
atom         92   22    C     "Tryptophan CD1"               6    12.011    3
atom         93   23    H     "Tryptophan HD1"               1     1.008    1
atom         94   22    C     "Tryptophan CD2"               6    12.011    3
atom         95   24    N     "Tryptophan NE1"               7    14.007    3
atom         96   25    HN    "Tryptophan HE1"               1     1.008    1
atom         97   22    C     "Tryptophan CE2"               6    12.011    3
atom         98   22    C     "Tryptophan CE3"               6    12.011    3
atom         99   23    H     "Tryptophan HE3"               1     1.008    1
atom        100   22    C     "Tryptophan CZ2"               6    12.011    3
atom        101   23    H     "Tryptophan HZ2"               1     1.008    1
atom        102   22    C     "Tryptophan CZ3"               6    12.011    3
atom        103   23    H     "Tryptophan HZ3"               1     1.008    1
atom        104   22    C     "Tryptophan CH2"               6    12.011    3
atom        105   23    H     "Tryptophan HH2"               1     1.008    1
atom        106    8    C     "Histidine (+) CB"             6    12.011    4
atom        107    9    H     "Histidine (+) HB"             1     1.008    1
atom        108   26    C     "Histidine (+) CG"             6    12.011    3
atom        109   24    N     "Histidine (+) ND1"            7    14.007    3
atom        110   25    HN    "Histidine (+) HD1"            1     1.008    1
atom        111   26    C     "Histidine (+) CD2"            6    12.011    3
atom        112   27    H     "Histidine (+) HD2"            1     1.008    1
atom        113   28    C     "Histidine (+) CE1"            6    12.011    3
atom        114   27    H     "Histidine (+) HE1"            1     1.008    1
atom        115   24    N     "Histidine (+) NE2"            7    14.007    3
atom        116   25    HN    "Histidine (+) HE2"            1     1.008    1
atom        117    8    C     "Histidine (HD) CB"            6    12.011    4
atom        118    9    H     "Histidine (HD) HB"            1     1.008    1
atom        119   26    C     "Histidine (HD) CG"            6    12.011    3
atom        120   24    N     "Histidine (HD) ND1"           7    14.007    3
atom        121   25    HN    "Histidine (HD) HD1"           1     1.008    1
atom        122   26    C     "Histidine (HD) CD2"           6    12.011    3
atom        123   27    H     "Histidine (HD) HD2"           1     1.008    1
atom        124   28    C     "Histidine (HD) CE1"           6    12.011    3
atom        125   27    H     "Histidine (HD) HE1"           1     1.008    1
atom        126   29    N     "Histidine (HD) NE2"           7    14.007    2
atom        127    8    C     "Histidine (HE) CB"            6    12.011    4
atom        128    9    H     "Histidine (HE) HB"            1     1.008    1
atom        129   26    C     "Histidine (HE) CG"            6    12.011    3
atom        130   29    N     "Histidine (HE) ND1"           7    14.007    2
atom        131   26    C     "Histidine (HE) CD2"           6    12.011    3
atom        132   27    H     "Histidine (HE) HD2"           1     1.008    1
atom        133   28    C     "Histidine (HE) CE1"           6    12.011    3
atom        134   27    H     "Histidine (HE) HE1"           1     1.008    1
atom        135   24    N     "Histidine (HE) NE2"           7    14.007    3
atom        136   25    HN    "Histidine (HE) HE2"           1     1.008    1
atom        137    8    C     "Aspartate CB"                 6    12.011    4
atom        138    9    H     "Aspartate HB"                 1     1.008    1
atom        139   30    C     "Aspartate CG"                 6    12.011    3
atom        140   31    O     "Aspartate OD"                 8    15.999    1
atom        141    8    C     "Aspartic Acid CB"             6    12.011    4
atom        142    9    H     "Aspartic Acid HB"             1     1.008    1
atom        143   32    C     "Aspartic Acid CG"             6    12.011    3
atom        144   33    O     "Aspartic Acid OD1"            8    15.999    1
atom        145   34    OH    "Aspartic Acid OD2"            8    15.999    2
atom        146   35    HO    "Aspartic Acid HD2"            1     1.008    1
atom        147    8    C     "Asparagine CB"                6    12.011    4
atom        148    9    H     "Asparagine HB"                1     1.008    1
atom        149    3    C     "Asparagine CG"                6    12.011    3
atom        150    5    O     "Asparagine OD1"               8    15.999    1
atom        151    1    N     "Asparagine ND2"               7    14.007    3
atom        152    4    HN    "Asparagine HD2"               1     1.008    1
atom        153    8    C     "Glutamate CB"                 6    12.011    4
atom        154    9    H     "Glutamate HB"                 1     1.008    1
atom        155    8    C     "Glutamate CG"                 6    12.011    4
atom        156    9    H     "Glutamate HG"                 1     1.008    1
atom        157   30    C     "Glutamate CD"                 6    12.011    3
atom        158   31    O     "Glutamate OE"                 8    15.999    1
atom        159    8    C     "Glutamic Acid CB"             6    12.011    4
atom        160    9    H     "Glutamic Acid HB"             1     1.008    1
atom        161    8    C     "Glutamic Acid CG"             6    12.011    4
atom        162    9    H     "Glutamic Acid HG"             1     1.008    1
atom        163   32    C     "Glutamic Acid CD"             6    12.011    3
atom        164   33    O     "Glutamic Acid OE1"            8    15.999    1
atom        165   34    O     "Glutamic Acid OE2"            8    15.999    2
atom        166   35    H     "Glutamic Acid HE2"            8     1.008    1
atom        167    8    C     "Glutamine CB"                 6    12.011    4
atom        168    9    H     "Glutamine HB"                 1     1.008    1
atom        169    8    C     "Glutamine CG"                 6    12.011    4
atom        170    9    H     "Glutamine HG"                 1     1.008    1
atom        171    3    C     "Glutamine CD"                 6    12.011    3
atom        172    5    O     "Glutamine OE1"                8    15.999    1
atom        173    1    N     "Glutamine NE2"                7    14.007    3
atom        174    4    HN    "Glutamine HE2"                1     1.008    1
atom        175    8    C     "Methionine CB"                6    12.011    4
atom        176    9    H     "Methionine HB"                1     1.008    1
atom        177    8    C     "Methionine CG"                6    12.011    4
atom        178    9    H     "Methionine HG"                1     1.008    1
atom        179   12    S     "Methionine SD"               16    32.066    2
atom        180    8    C     "Methionine CE"                6    12.011    4
atom        181    9    H     "Methionine HE"                1     1.008    1
atom        182    8    C     "Lysine CB"                    6    12.011    4
atom        183    9    H     "Lysine HB"                    1     1.008    1
atom        184   36    C     "Lysine CG"                    6    12.011    4
atom        185    9    H     "Lysine HG"                    1     1.008    1
atom        186    8    C     "Lysine CD"                    6    12.011    4
atom        187    9    H     "Lysine HD"                    1     1.008    1
atom        188    8    C     "Lysine CE"                    6    12.011    4
atom        189    9    H     "Lysine HE"                    1     1.008    1
atom        190   37    N     "Lysine NZ"                    7    14.007    4
atom        191   38    HN    "Lysine HN"                    1     1.008    1
atom        192    8    C     "Lysine (Neutral) CB"          6    12.011    4
atom        193    9    H     "Lysine (Neutral) HB"          1     1.008    1
atom        194   36    C     "Lysine (Neutral) CG"          6    12.011    4
atom        195    9    H     "Lysine (Neutral) HG"          1     1.008    1
atom        196    8    C     "Lysine (Neutral) CD"          6    12.011    4
atom        197    9    H     "Lysine (Neutral) HD"          1     1.008    1
atom        198    8    C     "Lysine (Neutral) CE"          6    12.011    4
atom        199    9    H     "Lysine (Neutral) HE"          1     1.008    1
atom        200   37    N     "Lysine (Neutral) NZ"          7    14.007    3
atom        201   38    HN    "Lysine (Neutral) HN"          1     1.008    1
atom        202    8    C     "Arginine CB"                  6    12.011    4
atom        203    9    H     "Arginine HB"                  1     1.008    1
atom        204    8    C     "Arginine CG"                  6    12.011    4
atom        205    9    H     "Arginine HG"                  1     1.008    1
atom        206    8    C     "Arginine CD"                  6    12.011    4
atom        207    9    H     "Arginine HD"                  1     1.008    1
atom        208    1    N     "Arginine NE"                  7    14.007    3
atom        209    4    HN    "Arginine HE"                  1     1.008    1
atom        210   39    C     "Arginine CZ"                  6    12.011    3
atom        211    1    N     "Arginine NH"                  7    14.007    3
atom        212    4    HN    "Arginine HH"                  1     1.008    1
atom        213    8    C     "Ornithine CB"                 6    12.011    4
atom        214    9    H     "Ornithine HB"                 1     1.008    1
atom        215   36    C     "Ornithine CG"                 6    12.011    4
atom        216    9    H     "Ornithine HG"                 1     1.008    1
atom        217    8    C     "Ornithine CD"                 6    12.011    4
atom        218    9    H     "Ornithine HD"                 1     1.008    1
atom        219   37    N     "Ornithine NE"                 7    14.007    4
atom        220   38    HN    "Ornithine HE"                 1     1.008    1
atom        221   40    C     "Acetyl Cap CH3"               6    12.011    4
atom        222    6    H     "Acetyl Cap H3C"               1     1.008    1
atom        223    3    C     "Acetyl Cap C"                 6    12.011    3
atom        224    5    O     "Acetyl Cap O"                 8    15.999    1
atom        225    1    N     "Amide Cap NH2"                7    14.007    3
atom        226    4    HN    "Amide Cap H2N"                1     1.008    1
atom        227    1    N     "N-MeAmide Cap N"              7    14.007    3
atom        228    4    HN    "N-MeAmide Cap HN"             1     1.008    1
atom        229   40    C     "N-MeAmide Cap CH3"            6    12.011    4
atom        230    6    H     "N-MeAmide Cap H3C"            1     1.008    1
atom        231   41    N     "N-Terminal NH3+"              7    14.007    4
atom        232   42    H     "N-Terminal H3N+"              1     1.008    1
atom        233   30    C     "C-Terminal COO-"              6    12.011    3
atom        234   31    O     "C-Terminal COO-"              8    15.999    1
atom        235   32    C     "C-Terminal COOH C=O"          6    12.011    3
atom        236   33    O     "C-Terminal COOH O=C"          8    15.999    1
atom        237   34    OH    "C-Terminal COOH OH"           8    15.999    2
atom        238   35    HO    "C-Terminal COOH HO"           1     1.008    1
atom        239   41    N     "N-Terminal PRO NH2+"          7    14.007    4
atom        240   42    HN    "N-Terminal PRO H2N+"          1     1.008    1
atom        241    7    CA    "N-Terminal PRO CA"            6    12.011    4
atom        242    3    C     "N-Terminal PRO C"             6    12.011    3
atom        243    5    O     "N-Terminal PRO O"             8    15.999    1
atom        244    6    H     "N-Terminal PRO HA"            1     1.008    1
atom        245   16    C     "N-Terminal PRO CD"            6    12.011    4
atom        246    6    H     "N-Terminal PRO HD"            1     1.008    1
atom        247   43    O     "AMOEBA Water O"               8    15.999    2
atom        248   44    H     "AMOEBA Water H"               1     1.008    1
atom        249   45    Li+   "Lithium Ion Li+"              3     6.941    0
atom        250   46    Na+   "Sodium Ion Na+"              11    22.990    0
atom        251   47    K+    "Potassium Ion K+"            19    39.098    0
atom        252   48    Rb+   "Rubidium Ion Rb+"            37    85.468    0
atom        253   49    Cs+   "Cesium Ion Cs+"              55   132.905    0
atom        254   50    Be+   "Beryllium Ion Be+2"           4     9.012    0
atom        255   51    Mg+   "Magnesium Ion Mg+2"          12    24.305    0
atom        256   52    Ca+   "Calcium Ion Ca+2"            20    40.078    0
atom        257   53    Zn+   "Zinc Ion Zn+2"               30    65.380    0
atom        258   54    F-    "Fluoride Ion F-"              9    18.998    0
atom        259   55    Cl-   "Chloride Ion Cl-"            17    35.453    0
atom        260   56    Br-   "Bromide Ion Br-"             35    79.904    0
atom        261   57    I-    "Iodide Ion I-"               53   126.904    0


      ################################
      ##                            ##
      ##  Van der Waals Parameters  ##
      ##                            ##
      ################################


vdw           1               3.7100     0.1100
vdw           2               3.8200     0.1010
vdw           3               3.8200     0.1060
vdw           4               2.5900     0.0220      0.900
vdw           5               3.3000     0.1120
vdw           6               2.9400     0.0260      0.910
vdw           7               3.6500     0.1010
vdw           8               3.8200     0.1010
vdw           9               2.9800     0.0240      0.920
vdw          10               3.4050     0.1100
vdw          11               2.6550     0.0135      0.910
vdw          12               4.0050     0.3550
vdw          14               4.2000     0.3550
vdw          13               3.0000     0.0265      0.980
vdw          15               3.7100     0.1100
vdw          16               3.8200     0.1010
vdw          17               3.8000     0.0890
vdw          18               2.9800     0.0260      0.920
vdw          19               3.4050     0.1100
vdw          20               2.6550     0.0135      0.910
vdw          21               3.3200     0.1120
vdw          22               3.8000     0.1010
vdw          23               2.9800     0.0260      0.920
vdw          24               3.7100     0.1100
vdw          25               2.5900     0.0220      0.900
vdw          26               3.8000     0.1010
vdw          27               2.9800     0.0260      0.920
vdw          28               3.8000     0.1010
vdw          29               3.7100     0.1100
vdw          30               3.8200     0.1060
vdw          31               3.5500     0.0950
vdw          32               3.8200     0.1060
vdw          33               3.3000     0.1120
vdw          34               3.4050     0.1100
vdw          35               2.6550     0.0150      0.910
vdw          36               3.8200     0.1010
vdw          37               3.8100     0.1050
vdw          38               2.4800     0.0130      0.910
vdw          39               3.6500     0.1010
vdw          40               3.8200     0.1010
vdw          41               3.7600     0.1050
vdw          42               2.7000     0.0200      0.910
vdw          43               3.4050     0.1100
vdw          44               2.6550     0.0135      0.910
vdw          45               2.2000     0.0660
vdw          46               2.9550     0.2800
vdw          47               3.6800     0.3500
vdw          48               3.9000     0.3800
vdw          49               4.1400     0.4200
vdw          50               1.8800     0.0910
vdw          51               2.9000     0.2800
vdw          52               3.5900     0.3500
vdw          53               2.6800     0.2220
vdw          54               3.4300     0.2500
vdw          55               4.1200     0.3400
vdw          56               4.3200     0.4300
vdw          57               4.6100     0.5200


      #####################################
      ##                                 ##
      ##  Van der Waals Pair Parameters  ##
      ##                                 ##
      #####################################


vdwpr         4   31          3.1000     0.0400
vdwpr         5   51          3.0000     0.1530
vdwpr         5   52          3.2700     0.1750
vdwpr         5   53          2.8900     0.1750
vdwpr        47   55          4.2360     0.1512
vdwpr        47   56          4.3790     0.1664
vdwpr        47   57          4.6360     0.1720
vdwpr        48   55          4.3150     0.1859
vdwpr        48   56          4.4480     0.2068
vdwpr        48   57          4.6900     0.2145
vdwpr        49   55          4.3450     0.2894
vdwpr        49   56          4.4750     0.3307
vdwpr        49   57          4.7110     0.3466


      ##################################
      ##                              ##
      ##  Bond Stretching Parameters  ##
      ##                              ##
      ##################################


bond          1    2          375.00     1.4370
bond          1    3          482.00     1.3450
bond          1    4          487.00     1.0280
bond          1    7          375.00     1.4370
bond          1    8          374.80     1.4460
bond          1   39          491.40     1.3250
bond          1   40          375.00     1.4370
bond          2    3          345.00     1.5090
bond          2    6          341.00     1.1120
bond          2   30          345.00     1.5090
bond          2   32          345.00     1.5090
bond          2   41          381.30     1.4480
bond          3    5          662.00     1.2255
bond          3    7          345.00     1.5090
bond          3    8          345.00     1.5090
bond          3   15          482.00     1.3450
bond          3   40          345.00     1.5090
bond          6    7          341.00     1.1120
bond          6   16          341.00     1.1120
bond          6   40          341.00     1.1120
bond          7    7          323.00     1.5250
bond          7    8          323.00     1.5250
bond          7   10          410.00     1.4130
bond          7   15          375.00     1.4370
bond          7   16          323.00     1.5250
bond          7   30          345.00     1.5090
bond          7   32          345.00     1.5090
bond          7   41          381.30     1.4480
bond          8    8          323.00     1.5250
bond          8    9          341.00     1.1120
bond          8   10          410.00     1.4130
bond          8   12          215.80     1.8050
bond          8   14          215.80     1.8130
bond          8   17          453.20     1.4990
bond          8   22          345.00     1.5090
bond          8   26          453.20     1.4930
bond          8   30          345.00     1.5090
bond          8   32          345.00     1.5090
bond          8   36          323.00     1.5250
bond          8   37          381.30     1.4480
bond          9   16          341.00     1.1120
bond          9   36          341.00     1.1120
bond         10   11          548.90     0.9470
bond         12   12          188.50     2.0190
bond         12   13          278.40     1.3420
bond         15   16          375.00     1.4370
bond         16   16          323.00     1.5247
bond         16   41          381.30     1.4480
bond         17   17          471.90     1.3887
bond         17   18          370.50     1.1000
bond         17   19          431.60     1.3550
bond         17   21          680.00     1.2747
bond         19   20          548.90     0.9470
bond         22   22          471.90     1.3887
bond         22   23          370.50     1.1010
bond         22   24          653.90     1.3550
bond         24   25          467.60     1.0300
bond         24   26          653.90     1.3730
bond         24   28          653.90     1.3520
bond         26   26          539.60     1.3710
bond         26   27          370.50     1.0810
bond         26   29          670.00     1.3740
bond         27   28          370.50     1.0810
bond         28   29          670.00     1.3270
bond         30   31          705.00     1.2553
bond         32   33          705.00     1.2255
bond         32   34          431.60     1.3498
bond         34   35          514.40     0.9737
bond         37   38          461.90     1.0150
bond         41   42          461.90     1.0150
bond         43   44          556.85     0.9572


      ################################
      ##                            ##
      ##  Angle Bending Parameters  ##
      ##                            ##
      ################################


angle         2    1    3      50.00     122.00
angle         2    1    4      32.00     117.00
angle         3    1    4      36.00     121.00
angle         3    1    7      50.00     122.00
angle         3    1   40      50.00     121.00
angle         4    1    4      29.50     123.00
angle         4    1    7      32.00     117.00
angle         4    1    8      13.70     122.40
angle         4    1   39      41.70     120.50
angle         4    1   40      32.00     117.00
angle         8    1   39      52.50     120.40
angle         1    2    3      55.00     109.80
angle         1    2    6      55.00     111.00
angle         1    2   30      61.00     109.60
angle         1    2   32      61.00     109.60
angle         3    2    6      39.00     109.50
angle         3    2   41      75.20     110.70
angle         6    2    6      40.00     107.60
angle         6    2   30      39.00     109.50
angle         6    2   32      39.00     109.50
angle         6    2   41      59.00     109.30
angle         1    3    2      70.00     113.40
angle         1    3    5      77.00     124.20
angle         1    3    7      70.00     113.40
angle         1    3    8      70.00     113.40
angle         1    3   40      41.00     113.40
angle         2    3    5      80.00     122.40
angle         2    3   15      70.00     113.40
angle         5    3    7      80.00     122.40
angle         5    3    8      80.00     122.40
angle         5    3   15      77.00     124.20
angle         5    3   40      80.00     122.40
angle         7    3   15      70.00     113.40
angle        15    3   40      41.00     113.40
angle         1    7    3      61.00     109.60
angle         1    7    6      55.00     111.00
angle         1    7    7      54.00     111.30
angle         1    7    8      54.00     111.30
angle         1    7   30      61.00     109.60
angle         1    7   32      61.00     109.60
angle         3    7    6      39.00     109.50
angle         3    7    7      58.00     110.60
angle         3    7    8      58.00     110.60
angle         3    7   15      61.00     109.60
angle         3    7   16      58.00     110.60
angle         3    7   41      75.20     110.70
angle         6    7    7      42.00     110.70
angle         6    7    8      42.00     112.80
angle         6    7   10      59.00     110.00     108.90     108.70
angle         6    7   15      55.00     111.00
angle         6    7   16      42.00     112.80
angle         6    7   30      39.00     109.50
angle         6    7   32      39.00     109.50
angle         6    7   41      59.00     109.30
angle         7    7    8      48.20     109.50     110.20     111.00
angle         7    7   10      59.70     107.50
angle         7    7   30      58.00     110.60
angle         7    7   41      56.10     108.50
angle         8    7    8      48.20     109.50     110.20     111.00
angle         8    7   10      59.70     107.50
angle         8    7   30      58.00     110.60
angle         8    7   32      58.00     110.60
angle         8    7   41      75.20     110.70
angle        15    7   16      54.00     104.00
angle        15    7   30      61.00     109.60
angle        16    7   30      58.00     110.60
angle        16    7   41      54.00     104.00
angle         1    8    8      56.10     109.50
angle         1    8    9      54.60     111.00
angle         3    8    7      47.60     110.60
angle         3    8    8      47.60     110.60
angle         3    8    9      39.00     109.50
angle         7    8    7      48.20     109.50     110.20     111.00
angle         7    8    8      48.20     109.50     110.20     111.00
angle         7    8    9      42.00     110.70
angle         7    8   10      59.70     107.50
angle         7    8   12      53.20     108.00     109.50     110.10
angle         7    8   14      53.20     108.00     109.50     110.10
angle         7    8   17      38.90     110.60
angle         7    8   22      48.20     109.50     110.20     111.00
angle         7    8   26      38.80     112.70
angle         7    8   30      47.60     110.60
angle         7    8   32      47.60     110.60
angle         7    8   36      48.20     109.50     110.20     111.00
angle         8    8    8      48.20     109.50     110.20     111.00
angle         8    8    9      42.00     110.70
angle         8    8   12      53.20     108.00     109.50     110.10
angle         8    8   30      47.60     110.60
angle         8    8   32      47.60     110.60
angle         8    8   36      48.20     109.50     110.20     111.00
angle         8    8   37      56.10     109.50
angle         9    8    9      40.00     107.80
angle         9    8   10      59.00     110.00     108.90     108.70
angle         9    8   12      53.20     110.80     110.80     108.00
angle         9    8   14      53.20     110.80     110.80     108.00
angle         9    8   17      39.60     109.50     109.30     110.40
angle         9    8   22      61.20     109.80     109.30     110.70
angle         9    8   26      39.60     109.50
angle         9    8   30      39.00     109.50
angle         9    8   32      39.00     109.50
angle         9    8   36      42.00     110.70
angle         9    8   37      59.00     109.30
angle        36    8   37      56.10     109.50
angle         7   10   11      54.00     106.80
angle         8   10   11      54.00     106.80
angle         8   12    8      60.40      95.90
angle         8   12   12      71.90     101.80
angle         8   12   13      46.80      96.00
angle         3   15    7      50.00     122.00
angle         3   15   16      50.00     122.00
angle         7   15   16      54.70     112.50
angle         6   16    6      40.00     107.80
angle         6   16   15      55.00     111.00
angle         6   16   16      42.00     110.70
angle         6   16   41      55.00     111.00
angle         7   16    9      42.00     110.70
angle         7   16   16      48.20     104.00
angle         9   16    9      40.00     107.80
angle         9   16   16      42.00     110.70
angle        15   16   16      54.00     104.00
angle        16   16   16      48.20     104.00
angle        16   16   41      54.00     104.00
angle         8   17   17      33.80     122.30
angle        17   17   17      54.70     121.70
angle        17   17   18      35.30     120.00     120.50       0.00
angle        17   17   19      43.20     120.00
angle        17   17   21      60.00     123.57
angle        17   19   20      25.90     109.00
angle         8   22   22      61.90     127.00
angle        22   22   22      63.30     120.00
angle        22   22   23      35.30     128.00
angle        22   22   24      47.50     109.00
angle        23   22   24      86.30     122.50
angle        22   24   22      86.30     109.00
angle        22   24   25      35.30     125.50
angle        25   24   26      35.30     125.50
angle        25   24   28      35.30     125.50
angle        26   24   28      86.30     110.80
angle         8   26   24      36.00     122.00
angle         8   26   26      36.00     131.00
angle         8   26   29      36.00     122.00
angle        24   26   26      47.50     104.70
angle        24   26   27      38.10     122.50
angle        24   26   29      28.80     111.50
angle        26   26   27      35.30     128.00
angle        26   26   29      47.50     110.50
angle        27   26   29      38.10     122.50
angle        24   28   24      28.80     110.30
angle        24   28   27      38.10     122.50
angle        24   28   29      28.80     112.20
angle        27   28   29      38.10     122.50
angle        26   29   28      86.30     104.30
angle         2   30   31      80.00     122.40
angle         7   30   31      80.00     122.40
angle         8   30   31      80.00     122.40
angle         8   30   33      80.00     122.40
angle         8   30   34     111.50     110.30
angle        31   30   31      57.60     134.00
angle         2   32   33      80.00     122.40
angle         2   32   34      80.00     122.40
angle         7   32   33      80.00     122.40
angle         7   32   34      80.00     122.40
angle         8   32   33      80.00     122.40
angle         8   32   34     111.50     110.30
angle        33   32   34     122.30     121.50
angle        30   34   35      49.60     108.70
angle        32   34   35      49.60     108.70
angle         8   36    8      48.20     109.50     110.20     111.00
angle         8   36    9      42.00     110.70
angle         9   36    9      40.00     107.80
angle         8   37   38      43.20     110.90
angle        38   37   38      43.50     107.00
angle         1   39    1      28.80     120.00
angle         1   40    6      55.00     111.00
angle         3   40    6      39.00     109.50
angle         6   40    6      40.00     107.80
angle         2   41   42      43.20     108.50
angle         7   41   16      54.70     112.50
angle         7   41   42      43.20     108.50
angle        16   41   42      43.20     110.90
angle        42   41   42      43.20     105.40
angle        44   43   44      48.70     108.50


      ###############################
      ##                           ##
      ##  Stretch-Bend Parameters  ##
      ##                           ##
      ###############################


strbnd        2    1    3       7.20       7.20
strbnd        2    1    4       4.30       4.30
strbnd        3    1    4       4.30       4.30
strbnd        3    1    7       7.20       7.20
strbnd        3    1   40       7.20       7.20
strbnd        4    1    7       4.30       4.30
strbnd        4    1    8       4.30       4.30
strbnd        4    1   39       4.30       4.30
strbnd        4    1   40       4.30       4.30
strbnd        8    1   39       7.20       7.20
strbnd        1    2    3      18.70      18.70
strbnd        1    2    6      11.50      11.50
strbnd        1    2   30      18.70      18.70
strbnd        1    2   32      18.70      18.70
strbnd        3    2    6      11.50      11.50
strbnd        3    2   41      18.70      18.70
strbnd        6    2   30      11.50      11.50
strbnd        6    2   32      11.50      11.50
strbnd        6    2   41      11.50      11.50
strbnd        1    3    2      18.70      18.70
strbnd        1    3    5      18.70      18.70
strbnd        1    3    7      18.70      18.70
strbnd        1    3    8      18.70      18.70
strbnd        1    3   40      18.70      18.70
strbnd        2    3    5      18.70      18.70
strbnd        2    3   15      18.70      18.70
strbnd        5    3    7      18.70      18.70
strbnd        5    3    8      18.70      18.70
strbnd        5    3   15      18.70      18.70
strbnd        5    3   40      18.70      18.70
strbnd        7    3   15      18.70      18.70
strbnd       15    3   40      18.70      18.70
strbnd        1    7    3      18.70      18.70
strbnd        1    7    6      11.50      11.50
strbnd        1    7    7      18.70      18.70
strbnd        1    7    8      18.70      18.70
strbnd        1    7   30      18.70      18.70
strbnd        1    7   32      18.70      18.70
strbnd        3    7    6      11.50      11.50
strbnd        3    7    7      18.70      18.70
strbnd        3    7    8      18.70      18.70
strbnd        3    7   15      18.70      18.70
strbnd        3    7   16      18.70      18.70
strbnd        3    7   41      18.70      18.70
strbnd        6    7    7      11.50      11.50
strbnd        6    7    8      11.50      11.50
strbnd        6    7   10      11.50      11.50
strbnd        6    7   15      11.50      11.50
strbnd        6    7   16      11.50      11.50
strbnd        6    7   30      11.50      11.50
strbnd        6    7   32      11.50      11.50
strbnd        6    7   41      11.50      11.50
strbnd        7    7    8      18.70      18.70
strbnd        7    7   10      18.70      18.70
strbnd        7    7   30      18.70      18.70
strbnd        7    7   41      18.70      18.70
strbnd        8    7    8      18.70      18.70
strbnd        8    7   10      18.70      18.70
strbnd        8    7   30      18.70      18.70
strbnd        8    7   32      18.70      18.70
strbnd        8    7   41      18.70      18.70
strbnd       15    7   16      18.70      18.70
strbnd       15    7   30      18.70      18.70
strbnd       16    7   30      18.70      18.70
strbnd       16    7   41      18.70      18.70
strbnd        1    8    8      18.70      18.70
strbnd        1    8    9      11.50      11.50
strbnd        3    8    7      18.70      18.70
strbnd        3    8    8      18.70      18.70
strbnd        3    8    9      11.50      11.50
strbnd        7    8    7      18.70      18.70
strbnd        7    8    8      18.70      18.70
strbnd        7    8    9      11.50      11.50
strbnd        7    8   10      18.70      18.70
strbnd        7    8   12      18.70      18.70
strbnd        7    8   14      18.70      18.70
strbnd        7    8   17      18.70      18.70
strbnd        7    8   22      18.70      18.70
strbnd        7    8   26      18.70      18.70
strbnd        7    8   30      18.70      18.70
strbnd        7    8   32      18.70      18.70
strbnd        7    8   36      18.70      18.70
strbnd        8    8    8      18.70      18.70
strbnd        8    8    9      11.50      11.50
strbnd        8    8   12      18.70      18.70
strbnd        8    8   30      18.70      18.70
strbnd        8    8   32      18.70      18.70
strbnd        8    8   36      18.70      18.70
strbnd        8    8   37      18.70      18.70
strbnd        9    8   10      11.50      11.50
strbnd        9    8   12      11.50      11.50
strbnd        9    8   14      11.50      11.50
strbnd        9    8   17      11.50      11.50
strbnd        9    8   22      11.50      11.50
strbnd        9    8   26      11.50      11.50
strbnd        9    8   30      11.50      11.50
strbnd        9    8   32      11.50      11.50
strbnd        9    8   36      11.50      11.50
strbnd        9    8   37      11.50      11.50
strbnd       36    8   37      18.70      18.70
strbnd        7   10   11      12.95      12.95
strbnd        8   10   11      12.95      12.95
strbnd        8   12    8      -5.75      -5.75
strbnd        8   12   12      -5.75      -5.75
strbnd        8   12   13       1.45       1.45
strbnd        3   15    7       7.20       7.20
strbnd        3   15   16       7.20       7.20
strbnd        7   15   16       7.20       7.20
strbnd        6   16   15      11.50      11.50
strbnd        6   16   16      11.50      11.50
strbnd        6   16   41      11.50      11.50
strbnd        7   16    9      11.50      11.50
strbnd        7   16   16      18.70      18.70
strbnd        9   16   16      11.50      11.50
strbnd       15   16   16      18.70      18.70
strbnd       16   16   16      18.70      18.70
strbnd       16   16   41      18.70      18.70
strbnd        8   17   17      18.70      18.70
strbnd       17   17   17      18.70      18.70
strbnd       17   17   18      11.50      11.50
strbnd       17   17   19      18.70      18.70
strbnd       17   17   21      18.70      18.70
strbnd       17   19   20      12.95      12.95
strbnd        8   22   22      18.70      18.70
strbnd       22   22   22      18.70      18.70
strbnd       22   22   23      11.50      11.50
strbnd       22   22   24      18.70      18.70
strbnd       23   22   24      11.50      11.50
strbnd       22   24   22      14.40      14.40
strbnd       22   24   25       4.30       4.30
strbnd       25   24   26       4.30       4.30
strbnd       25   24   28       4.30       4.30
strbnd       26   24   28      14.40      14.40
strbnd        8   26   24      18.70      18.70
strbnd        8   26   26      18.70      18.70
strbnd        8   26   29      18.70      18.70
strbnd       24   26   26      18.70      18.70
strbnd       24   26   27      11.50      11.50
strbnd       24   26   29      18.70      18.70
strbnd       26   26   27      11.50      11.50
strbnd       26   26   29      18.70      18.70
strbnd       27   26   29      11.50      11.50
strbnd       24   28   24      18.70      18.70
strbnd       24   28   27      11.50      11.50
strbnd       24   28   29      18.70      18.70
strbnd       27   28   29      11.50      11.50
strbnd       26   29   28      14.40      14.40
strbnd        2   30   31      18.70      18.70
strbnd        7   30   31      18.70      18.70
strbnd        8   30   31      18.70      18.70
strbnd        8   30   33      18.70      18.70
strbnd        8   30   34      18.70      18.70
strbnd       31   30   31      18.70      18.70
strbnd        2   32   33      18.70      18.70
strbnd        2   32   34      18.70      18.70
strbnd        7   32   33      18.70      18.70
strbnd        7   32   34      18.70      18.70
strbnd        8   32   33      18.70      18.70
strbnd        8   32   34      18.70      18.70
strbnd       33   32   34      18.70      18.70
strbnd        8   36    8      18.70      18.70
strbnd        8   36    9      11.50      11.50
strbnd        9   36    9      11.50      11.50
strbnd        8   37   38       4.30       4.30
strbnd        1   39    1      18.70      18.70
strbnd        1   40    6      11.50      11.50
strbnd        3   40    6      11.50      11.50
strbnd        2   41   42       4.30       4.30
strbnd        7   41   16       7.20       7.20
strbnd        7   41   42       4.30       4.30
strbnd       16   41   42       4.30       4.30


      ###############################
      ##                           ##
      ##  Urey-Bradley Parameters  ##
      ##                           ##
      ###############################


ureybrad     44   43   44      -7.60     1.5537


      ####################################
      ##                                ##
      ##  Out-of-Plane Bend Parameters  ##
      ##                                ##
      ####################################


opbend        2    1    0    0            41.70
opbend        3    1    0    0            70.50
opbend        4    1    0    0            12.90
opbend        7    1    0    0            41.70
opbend        8    1    0    0             0.70
opbend       39    1    0    0             3.60
opbend       40    1    0    0            41.70
opbend        1    3    0    0           107.90
opbend        2    3    0    0            49.60
opbend        5    3    0    0            54.00
opbend        7    3    0    0            49.60
opbend        8    3    0    0            20.90
opbend       15    3    0    0            49.60
opbend       40    3    0    0            49.60
opbend        3   15    0    0            14.40
opbend        7   15    0    0            14.40
opbend       16   15    0    0            14.40
opbend       40   15    0    0            14.40
opbend        8   17    0    0            14.40
opbend       17   17    0    0            14.40
opbend       18   17    0    0            15.10
opbend       19   17    0    0            14.40
opbend       21   17    0    0            14.40
opbend        8   22    0    0            14.40
opbend       22   22    0    0            14.40
opbend       23   22    0    0            15.10
opbend       24   22    0    0            14.40
opbend       22   24    0    0            34.50
opbend       25   24    0    0            12.90
opbend       26   24    0    0            15.10
opbend       28   24    0    0            15.10
opbend        8   26    0    0            14.40
opbend       24   26    0    0            14.40
opbend       26   26    0    0            14.40
opbend       27   26    0    0            14.40
opbend       29   26    0    0            14.40
opbend       24   28    0    0            14.40
opbend       27   28    0    0            14.40
opbend       29   28    0    0            14.40
opbend        2   30    0    0           121.60
opbend        7   30    0    0           121.60
opbend        8   30    0    0           121.60
opbend       31   30    0    0           118.70
opbend        2   32    0    0           121.60
opbend        7   32    0    0           121.60
opbend        8   32    0    0           121.60
opbend       33   32    0    0           118.70
opbend       34   32    0    0           133.10
opbend        1   39    0    0            14.40


      ############################
      ##                        ##
      ##  Torsional Parameters  ##
      ##                        ##
      ############################


torsion       3    1    2    3     -3.805 0.0 1   1.646 180.0 2   1.239 0.0 3
torsion       3    1    2    6      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       3    1    2   30     -0.856 0.0 1   0.048 180.0 2  -1.842 0.0 3
torsion       3    1    2   32     -6.037 0.0 1  -0.396 180.0 2   2.176 0.0 3
torsion       4    1    2    3      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       4    1    2    6      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       4    1    2   30      0.000 0.0 1   8.000 180.0 2   0.000 0.0 3
torsion       4    1    2   32      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       2    1    3    2     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       2    1    3    5      1.000 0.0 1   2.250 180.0 2  -2.250 0.0 3
torsion       2    1    3    7     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       2    1    3   40     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       4    1    3    2      0.000 0.0 1   1.000 180.0 2   0.800 0.0 3
torsion       4    1    3    5      0.000 0.0 1   1.000 180.0 2  -0.550 0.0 3
torsion       4    1    3    7      0.000 0.0 1   1.000 180.0 2   0.800 0.0 3
torsion       4    1    3    8      0.000 0.0 1   1.000 180.0 2   0.800 0.0 3
torsion       4    1    3   40      0.000 0.0 1   1.000 180.0 2   0.800 0.0 3
torsion       7    1    3    2     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       7    1    3    5      1.000 0.0 1   2.250 180.0 2  -2.250 0.0 3
torsion       7    1    3    7     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       7    1    3   40     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion      40    1    3    2     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion      40    1    3    5      1.000 0.0 1   2.250 180.0 2  -2.250 0.0 3
torsion      40    1    3    7     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       3    1    7    3     -1.407 0.0 1  -0.068 180.0 2   0.000 0.0 3
torsion       3    1    7    6      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion       3    1    7    7      2.576 0.0 1   1.011 180.0 2   0.825 0.0 3
torsion       3    1    7    8      2.576 0.0 1   1.011 180.0 2   0.825 0.0 3
torsion       3    1    7   30     -0.856 0.0 1   0.048 180.0 2  -1.842 0.0 3
torsion       3    1    7   32     -2.296 0.0 1  -0.249 180.0 2  -0.097 0.0 3
torsion       4    1    7    3      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       4    1    7    6      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       4    1    7    7      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       4    1    7    8      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       4    1    7   30      0.000 0.0 1   8.000 180.0 2   0.000 0.0 3
torsion       4    1    7   32      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       4    1    8    8      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       4    1    8    9      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      39    1    8    8     -0.150 0.0 1   0.550 180.0 2  -0.450 0.0 3
torsion      39    1    8    9      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion       4    1   39    1      0.000 0.0 1   4.000 180.0 2   0.000 0.0 3
torsion       8    1   39    1      0.000 0.0 1   4.500 180.0 2   0.000 0.0 3
torsion       3    1   40    6      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion       4    1   40    6      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    2    3    1     -2.411 0.0 1  -0.587 180.0 2   0.493 0.0 3
torsion       1    2    3    5      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    2    3   15      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       6    2    3    1      0.000 0.0 1   0.000 180.0 2  -0.010 0.0 3
torsion       6    2    3    5      0.000 0.0 1   0.000 180.0 2   0.235 0.0 3
torsion       6    2    3   15      0.000 0.0 1   0.000 180.0 2  -0.010 0.0 3
torsion      41    2    3    1     -0.397 0.0 1   2.196 180.0 2  -0.371 0.0 3
torsion      41    2    3    5      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      41    2    3   15     -0.397 0.0 1   2.196 180.0 2  -0.371 0.0 3
torsion       1    2   30   31     -0.059 0.0 1   3.400 180.0 2  -0.045 0.0 3
torsion       6    2   30   31     -0.154 0.0 1   0.044 180.0 2  -0.086 0.0 3
torsion       1    2   32   33      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    2   32   34     -0.523 0.0 1   3.984 180.0 2   0.545 0.0 3
torsion       6    2   32   33      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       6    2   32   34      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       3    2   41   42      1.239 0.0 1  -0.405 180.0 2   0.035 0.0 3
torsion       6    2   41   42      0.000 0.0 1  -0.081 180.0 2   0.370 0.0 3
torsion       1    3    7    1     -0.495 0.0 1   3.000 180.0 2  -0.792 0.0 3
torsion       1    3    7    6      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    3    7    7      0.929 0.0 1   0.328 180.0 2   0.000 0.0 3
torsion       1    3    7    8      0.929 0.0 1   0.328 180.0 2   0.000 0.0 3
torsion       1    3    7   15      0.476 0.0 1   1.100 180.0 2  -0.160 0.0 3
torsion       1    3    7   16      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    3    7   41     -0.397 0.0 1   2.196 180.0 2  -0.371 0.0 3
torsion       5    3    7    1      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       5    3    7    6      0.000 0.0 1   0.000 180.0 2   0.235 0.0 3
torsion       5    3    7    7      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       5    3    7    8      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       5    3    7   15      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       5    3    7   16      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       5    3    7   41      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      15    3    7    1     -1.804 0.0 1   3.612 180.0 2  -0.196 0.0 3
torsion      15    3    7    6      0.000 0.0 1   0.000 180.0 2  -0.010 0.0 3
torsion      15    3    7    7      0.355 0.0 1   1.478 180.0 2  -0.896 0.0 3
torsion      15    3    7    8      0.355 0.0 1   1.478 180.0 2  -0.896 0.0 3
torsion      15    3    7   15     -2.118 0.0 1   1.173 180.0 2   0.000 0.0 3
torsion      15    3    7   16      2.225 0.0 1  -1.252 180.0 2  -1.067 0.0 3
torsion      15    3    7   41     -0.397 0.0 1   2.196 180.0 2  -0.371 0.0 3
torsion       1    3    8    7      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    3    8    8      1.780 0.0 1   1.100 180.0 2  -1.270 0.0 3
torsion       1    3    8    9      0.000 0.0 1   0.000 180.0 2   0.230 0.0 3
torsion       5    3    8    7      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       5    3    8    8      0.454 0.0 1   0.080 180.0 2  -0.310 0.0 3
torsion       5    3    8    9     -0.154 0.0 1   0.044 180.0 2  -0.086 0.0 3
torsion       2    3   15    7     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       2    3   15   16     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       5    3   15    7      1.000 0.0 1   2.250 180.0 2  -2.250 0.0 3
torsion       5    3   15   16      1.000 0.0 1   2.250 180.0 2  -2.250 0.0 3
torsion       7    3   15    7     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       7    3   15   16     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion      40    3   15    7     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion      40    3   15   16     -1.000 0.0 1   2.000 180.0 2   2.000 0.0 3
torsion       1    3   40    6      0.000 0.0 1   0.000 180.0 2  -0.010 0.0 3
torsion       5    3   40    6      0.000 0.0 1   0.000 180.0 2   0.235 0.0 3
torsion      15    3   40    6      0.000 0.0 1   0.000 180.0 2  -0.010 0.0 3
torsion       1    7    7    6      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    7    7    8      0.800 0.0 1   1.000 180.0 2   0.500 0.0 3
torsion       1    7    7   10     -1.500 0.0 1   0.720 180.0 2   1.260 0.0 3
torsion       3    7    7    6      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       3    7    7    8      0.650 0.0 1   1.025 180.0 2   0.120 0.0 3
torsion       3    7    7   10     -3.490 0.0 1   2.340 180.0 2   0.000 0.0 3
torsion       6    7    7    6      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       6    7    7    8      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       6    7    7   10      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       6    7    7   30      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       6    7    7   41      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       8    7    7   30      0.650 0.0 1   1.025 180.0 2   0.120 0.0 3
torsion       8    7    7   41      0.800 0.0 1   1.000 180.0 2   0.500 0.0 3
torsion      10    7    7   29     -3.490 0.0 1   2.340 180.0 2   0.000 0.0 3
torsion      10    7    7   41     -1.500 0.0 1   0.720 180.0 2   1.260 0.0 3
torsion       1    7    8    3     -0.933 0.0 1  -0.929 180.0 2   0.153 0.0 3
torsion       1    7    8    7      0.902 0.0 1   0.520 180.0 2   1.000 0.0 3
torsion       1    7    8    8     -2.280 0.0 1   0.970 180.0 2   3.700 0.0 3
torsion       1    7    8    9      0.000 0.0 1   0.000 180.0 2   0.500 0.0 3
torsion       1    7    8   10      0.630 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    7    8   12     -0.160 0.0 1   1.080 180.0 2  -1.520 0.0 3
torsion       1    7    8   14     -6.170 0.0 1   1.440 180.0 2   0.760 0.0 3
torsion       1    7    8   17     -1.470 0.0 1   1.530 180.0 2   0.260 0.0 3
torsion       1    7    8   22     -1.000 0.0 1   0.100 180.0 2   1.150 0.0 3
torsion       1    7    8   26      0.730 0.0 1   1.780 180.0 2   0.000 0.0 3
torsion       1    7    8   30     -2.900 0.0 1   1.800 180.0 2   0.000 0.0 3
torsion       1    7    8   32     -4.150 0.0 1   2.060 180.0 2   0.000 0.0 3
torsion       1    7    8   36     -0.500 0.0 1   0.800 180.0 2  -1.110 0.0 3
torsion       3    7    8    3     -0.911 0.0 1   3.500 180.0 2  -0.198 0.0 3
torsion       3    7    8    7      0.740 0.0 1   0.270 180.0 2  -1.316 0.0 3
torsion       3    7    8    8      0.160 0.0 1   1.655 180.0 2  -2.520 0.0 3
torsion       3    7    8    9      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion       3    7    8   10      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       3    7    8   12     -1.010 0.0 1   1.230 180.0 2   1.000 0.0 3
torsion       3    7    8   14     -5.680 0.0 1   3.000 180.0 2   1.000 0.0 3
torsion       3    7    8   17      1.020 0.0 1   1.150 180.0 2  -0.300 0.0 3
torsion       3    7    8   22      0.700 0.0 1   1.800 180.0 2   0.000 0.0 3
torsion       3    7    8   26     -0.660 0.0 1   0.200 180.0 2   0.000 0.0 3
torsion       3    7    8   30     -6.950 0.0 1  -1.150 180.0 2   0.000 0.0 3
torsion       3    7    8   32      1.000 0.0 1  -0.130 180.0 2   0.000 0.0 3
torsion       3    7    8   36     -1.810 0.0 1   0.150 180.0 2   1.000 0.0 3
torsion       6    7    8    3      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion       6    7    8    7      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       6    7    8    8      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       6    7    8    9      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       6    7    8    9      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       6    7    8    9      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       6    7    8   10      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       6    7    8   12      0.000 0.0 1   0.000 180.0 2   0.475 0.0 3
torsion       6    7    8   14      0.000 0.0 1   0.000 180.0 2   0.475 0.0 3
torsion       6    7    8   17      0.000 0.0 1   0.000 180.0 2   0.500 0.0 3
torsion       6    7    8   22      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       6    7    8   26      0.000 0.0 1   0.000 180.0 2   0.500 0.0 3
torsion       6    7    8   30      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion       6    7    8   32      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion       6    7    8   36      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       7    7    8    8     -0.640 0.0 1   0.360 180.0 2  -0.089 0.0 3
torsion       7    7    8    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       8    7    8    7     -0.083 0.0 1  -0.045 180.0 2  -0.120 0.0 3
torsion       8    7    8    8     -0.640 0.0 1   0.360 180.0 2  -0.089 0.0 3
torsion       8    7    8    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       8    7    8    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       8    7    8    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion      10    7    8    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion      30    7    8    3      0.000 0.0 1   0.000 180.0 2   0.100 0.0 3
torsion      30    7    8    7      0.740 0.0 1   0.270 180.0 2  -1.316 0.0 3
torsion      30    7    8    8      1.190 0.0 1   0.520 180.0 2   0.520 0.0 3
torsion      30    7    8    8      1.190 0.0 1   0.520 180.0 2   0.520 0.0 3
torsion      30    7    8    9      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion      30    7    8    9      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion      30    7    8   10      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      30    7    8   12     -1.010 0.0 1   1.230 180.0 2   1.000 0.0 3
torsion      30    7    8   14     -5.680 0.0 1   3.000 180.0 2   1.000 0.0 3
torsion      30    7    8   17      1.020 0.0 1   1.150 180.0 2  -0.300 0.0 3
torsion      30    7    8   22      0.700 0.0 1   1.800 180.0 2   0.000 0.0 3
torsion      30    7    8   26     -0.660 0.0 1   0.200 180.0 2   0.000 0.0 3
torsion      30    7    8   30     -2.450 0.0 1   1.688 180.0 2   0.000 0.0 3
torsion      30    7    8   32      1.000 0.0 1  -0.130 180.0 2   0.000 0.0 3
torsion      30    7    8   36     -1.810 0.0 1   0.150 180.0 2   1.000 0.0 3
torsion      32    7    8    9      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion      41    7    8    3     -3.000 0.0 1   0.500 180.0 2   0.100 0.0 3
torsion      41    7    8    7      0.902 0.0 1   0.520 180.0 2   1.000 0.0 3
torsion      41    7    8    8     -0.513 0.0 1   0.568 180.0 2  -0.750 0.0 3
torsion      41    7    8    9      0.000 0.0 1   0.000 180.0 2   0.500 0.0 3
torsion      41    7    8   10      0.630 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      41    7    8   12     -0.160 0.0 1   1.080 180.0 2  -1.520 0.0 3
torsion      41    7    8   14     -6.170 0.0 1   1.440 180.0 2   0.760 0.0 3
torsion      41    7    8   17     -1.470 0.0 1   1.530 180.0 2   0.260 0.0 3
torsion      41    7    8   22     -1.000 0.0 1   0.100 180.0 2   1.150 0.0 3
torsion      41    7    8   26      0.730 0.0 1   1.780 180.0 2   0.000 0.0 3
torsion      41    7    8   30     -4.150 0.0 1   2.060 180.0 2   0.000 0.0 3
torsion      41    7    8   32     -3.400 0.0 1   0.900 180.0 2   0.000 0.0 3
torsion      41    7    8   36     -0.500 0.0 1   0.800 180.0 2  -1.110 0.0 3
torsion       6    7   10   11      0.000 0.0 1   0.000 180.0 2   0.250 0.0 3
torsion       7    7   10   11     -0.266 0.0 1  -0.910 180.0 2   0.322 0.0 3
torsion       8    7   10   11     -0.266 0.0 1  -0.910 180.0 2   0.322 0.0 3
torsion       3    7   15    3      0.000 0.0 1   0.000 180.0 2   0.030 0.0 3
torsion       3    7   15   16      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       6    7   15    3      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion       6    7   15   16      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion      16    7   15    3      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      16    7   15   16      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      30    7   15    3      0.000 0.0 1   0.000 180.0 2   0.030 0.0 3
torsion      30    7   15   16      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       3    7   16    9      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion       3    7   16   16      1.660 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       6    7   16    9      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       6    7   16   16      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion      15    7   16    9      0.000 0.0 1   0.000 180.0 2   0.500 0.0 3
torsion      15    7   16   16      0.201 0.0 1   0.687 180.0 2   0.000 0.0 3
torsion      30    7   16    9      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion      30    7   16   16      1.660 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      41    7   16    9      0.000 0.0 1   0.000 180.0 2   0.180 0.0 3
torsion      41    7   16   16      1.660 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    7   30   31     -0.059 0.0 1   1.900 180.0 2  -0.045 0.0 3
torsion       6    7   30   31     -0.154 0.0 1   0.044 180.0 2  -0.086 0.0 3
torsion       7    7   30   31      0.000 0.0 1  -0.489 180.0 2   0.000 0.0 3
torsion       8    7   30   31      0.000 0.0 1  -0.489 180.0 2   0.000 0.0 3
torsion      15    7   30   31      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      16    7   30   31      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       1    7   32   33      0.000 0.0 1  -0.136 180.0 2   0.000 0.0 3
torsion       1    7   32   34      0.000 0.0 1   1.645 180.0 2   0.000 0.0 3
torsion       6    7   32   33      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       6    7   32   34      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       8    7   32   33      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       8    7   32   34      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       3    7   41   16      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       3    7   41   42      1.239 0.0 1  -0.405 180.0 2   0.035 0.0 3
torsion       6    7   41   16      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion       6    7   41   42      0.000 0.0 1  -0.081 180.0 2   0.370 0.0 3
torsion       7    7   41   42      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       8    7   41   42      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion      16    7   41   16      2.216 0.0 1   0.215 180.0 2   0.810 0.0 3
torsion      16    7   41   42     -0.939 0.0 1  -0.159 180.0 2   0.246 0.0 3
torsion       1    8    8    8      0.000 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion       1    8    8    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       3    8    8    7      0.750 0.0 1  -1.240 180.0 2  -0.190 0.0 3
torsion       3    8    8    9      0.000 0.0 1   0.000 180.0 2   0.380 0.0 3
torsion       7    8    8    8     -1.100 0.0 1  -1.150 180.0 2   0.000 0.0 3
torsion       7    8    8    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       7    8    8   12     -0.500 0.0 1  -0.310 180.0 2  -1.290 0.0 3
torsion       7    8    8   30     -2.500 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       7    8    8   32      1.000 0.0 1   1.000 180.0 2   0.500 0.0 3
torsion       8    8    8    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       9    8    8    9      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       9    8    8   12      0.000 0.0 1   0.000 180.0 2   0.475 0.0 3
torsion       9    8    8   30      0.000 0.0 1   0.000 180.0 2   1.300 0.0 3
torsion       9    8    8   32      0.000 0.0 1   0.000 180.0 2   0.300 0.0 3
torsion       9    8    8   36      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       9    8    8   37      0.000 0.0 1   0.000 180.0 2   0.374 0.0 3
torsion      36    8    8   37     -0.800 0.0 1   1.064 180.0 2   0.000 0.0 3
torsion       7    8   10   11      3.000 0.0 1  -1.000 180.0 2   0.000 0.0 3
torsion       9    8   10   11      0.000 0.0 1   0.000 180.0 2   0.250 0.0 3
torsion       7    8   12   12      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       7    8   12   13      0.042 0.0 1  -0.091 180.0 2   0.068 0.0 3
torsion       8    8   12    8     -1.250 0.0 1   0.170 180.0 2   0.509 0.0 3
torsion       9    8   12    8      0.000 0.0 1   0.000 180.0 2   0.660 0.0 3
torsion       9    8   12   12      0.300 0.0 1   0.000 180.0 2   0.600 0.0 3
torsion       9    8   12   13      0.000 0.0 1   0.000 180.0 2   0.343 0.0 3
torsion       7    8   17   17      0.000 0.0 1  -0.574 180.0 2   0.000 0.0 3
torsion       9    8   17   17      0.000 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       7    8   22   22      0.000 0.0 1  -0.450 180.0 2   0.000 0.0 3
torsion       9    8   22   22      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       7    8   26   24     -2.000 0.0 1   0.000 180.0 2   1.000 0.0 3
torsion       7    8   26   26      0.350 0.0 1   0.280 180.0 2   0.000 0.0 3
torsion       7    8   26   29     -1.520 0.0 1   0.000 180.0 2   0.000 0.0 3
torsion       9    8   26   24      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       9    8   26   26      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       9    8   26   29      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       7    8   30   31      0.000 0.0 1   1.700 180.0 2   0.000 0.0 3
torsion       8    8   30   31      0.000 0.0 1   1.460 180.0 2   0.000 0.0 3
torsion       9    8   30   31     -0.154 0.0 1   0.044 180.0 2  -0.086 0.0 3
torsion       7    8   32   33     -0.092 0.0 1   1.124 180.0 2   0.435 0.0 3
torsion       7    8   32   34      0.677 0.0 1   1.020 180.0 2   0.163 0.0 3
torsion       8    8   32   33      0.000 0.0 1   0.800 180.0 2   0.000 0.0 3
torsion       8    8   32   34      0.000 0.0 1   0.800 180.0 2   0.000 0.0 3
torsion       9    8   32   33     -0.154 0.0 1   0.044 180.0 2  -0.086 0.0 3
torsion       9    8   32   34      0.250 0.0 1   0.850 180.0 2   0.000 0.0 3
torsion       7    8   36    8      1.000 0.0 1   0.070 180.0 2  -1.220 0.0 3
torsion       7    8   36    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       8    8   36    8      0.000 0.0 1   0.000 180.0 2  -1.000 0.0 3
torsion       8    8   36    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       9    8   36    8      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       9    8   36    9      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion      37    8   36    8     -0.800 0.0 1   1.064 180.0 2   0.000 0.0 3
torsion      37    8   36    9      0.000 0.0 1   0.000 180.0 2   0.374 0.0 3
torsion       8    8   37   38      0.000 0.0 1   0.000 180.0 2  -0.110 0.0 3
torsion       9    8   37   38      0.000 0.0 1  -0.081 180.0 2   0.370 0.0 3
torsion      36    8   37   38      0.000 0.0 1   0.000 180.0 2  -0.110 0.0 3
torsion       8   12   12    8     -0.940 0.0 1  -6.900 180.0 2   0.300 0.0 3
torsion       3   15   16    6      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion       3   15   16   16     -2.216 0.0 1   0.215 180.0 2  -0.810 0.0 3
torsion       7   15   16    6      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion       7   15   16   16     -2.216 0.0 1   0.215 180.0 2  -0.810 0.0 3
torsion       6   16   16    9      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       6   16   16   16      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       7   16   16    9      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       7   16   16   16     -2.000 0.0 1  -0.500 180.0 2   0.000 0.0 3
torsion       9   16   16    9      0.000 0.0 1   0.000 180.0 2   0.299 0.0 3
torsion       9   16   16   15      0.000 0.0 1   0.000 180.0 2   0.500 0.0 3
torsion       9   16   16   16      0.000 0.0 1   0.000 180.0 2   0.341 0.0 3
torsion       9   16   16   41      0.000 0.0 1   0.000 180.0 2   0.500 0.0 3
torsion      15   16   16   16     -2.201 0.0 1   0.687 180.0 2   0.000 0.0 3
torsion      16   16   16   41     -0.201 0.0 1   0.687 180.0 2   0.000 0.0 3
torsion       6   16   41    7      0.000 0.0 1   0.000 180.0 2  -0.126 0.0 3
torsion       6   16   41   42      0.000 0.0 1  -0.014 180.0 2   0.295 0.0 3
torsion      16   16   41    7      2.216 0.0 1   0.215 180.0 2   0.810 0.0 3
torsion      16   16   41   42     -0.939 0.0 1  -0.159 180.0 2   0.246 0.0 3
torsion       8   17   17   17     -0.610 0.0 1   4.212 180.0 2   0.000 0.0 3
torsion       8   17   17   18      0.000 0.0 1   4.104 180.0 2   0.000 0.0 3
torsion      17   17   17   17     -0.670 0.0 1   4.004 180.0 2   0.000 0.0 3
torsion      17   17   17   18      0.550 0.0 1   4.534 180.0 2  -0.550 0.0 3
torsion      17   17   17   19      0.000 0.0 1   4.470 180.0 2   0.000 0.0 3
torsion      17   17   17   21      0.000 0.0 1   4.470 180.0 2   0.000 0.0 3
torsion      18   17   17   18      0.000 0.0 1   4.072 180.0 2   0.000 0.0 3
torsion      18   17   17   19      0.000 0.0 1   4.470 180.0 2   0.000 0.0 3
torsion      18   17   17   21      0.000 0.0 1   4.470 180.0 2   0.000 0.0 3
torsion      17   17   19   20      0.000 0.0 1   2.300 180.0 2   0.000 0.0 3
torsion       8   22   22   22     -0.610 0.0 1   4.212 180.0 2   0.000 0.0 3
torsion       8   22   22   23      0.000 0.0 1   3.104 180.0 2   0.000 0.0 3
torsion       8   22   22   24      0.000 0.0 1   4.470 180.0 2   0.000 0.0 3
torsion      22   22   22   22     -0.670 0.0 1   4.304 180.0 2   0.000 0.0 3
torsion      22   22   22   23      0.250 0.0 1   4.534 180.0 2  -0.550 0.0 3
torsion      22   22   22   24      0.000 0.0 1   4.470 180.0 2   0.000 0.0 3
torsion      23   22   22   23      0.000 0.0 1   3.072 180.0 2   0.000 0.0 3
torsion      23   22   22   24     -3.150 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      22   22   24   22      0.000 0.0 1  10.000 180.0 2   0.000 0.0 3
torsion      22   22   24   25     -3.150 0.0 1   5.000 180.0 2   0.000 0.0 3
torsion      23   22   24   22     -6.650 0.0 1  10.000 180.0 2   0.000 0.0 3
torsion      23   22   24   25     -0.530 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      25   24   26    8      0.000 0.0 1   4.104 180.0 2   0.000 0.0 3
torsion      25   24   26   26     -3.150 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      25   24   26   27     -0.530 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      28   24   26    8      0.000 0.0 1   4.212 180.0 2   0.000 0.0 3
torsion      28   24   26   26      0.000 0.0 1   6.000 180.0 2   0.000 0.0 3
torsion      28   24   26   27      0.000 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      25   24   28   24     -2.744 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      25   24   28   27     -0.530 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      25   24   28   29     -2.744 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      26   24   28   24      0.000 0.0 1   8.000 180.0 2   0.000 0.0 3
torsion      26   24   28   27      0.000 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      26   24   28   29      0.000 0.0 1   8.000 180.0 2   0.000 0.0 3
torsion       8   26   26   24      0.000 0.0 1   4.212 180.0 2   0.000 0.0 3
torsion       8   26   26   27      0.000 0.0 1   4.102 180.0 2   0.000 0.0 3
torsion       8   26   26   29      0.000 0.0 1   4.212 180.0 2   0.000 0.0 3
torsion      24   26   26   24      0.900 0.0 1   8.000 180.0 2   0.000 0.0 3
torsion      24   26   26   27      0.755 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      24   26   26   29      0.900 0.0 1   8.000 180.0 2   0.000 0.0 3
torsion      27   26   26   29      0.755 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion       8   26   29   28      0.000 0.0 1   4.212 180.0 2   0.000 0.0 3
torsion      26   26   29   28      0.000 0.0 1   6.000 180.0 2   0.000 0.0 3
torsion      27   26   29   28      0.000 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion      24   28   29   26      0.000 0.0 1   8.000 180.0 2   0.000 0.0 3
torsion      27   28   29   26      0.000 0.0 1   3.000 180.0 2   0.000 0.0 3
torsion       2   32   34   35      0.000 0.0 1   5.390 180.0 2   1.230 0.0 3
torsion       7   32   34   35      0.000 0.0 1   5.390 180.0 2   1.230 0.0 3
torsion       8   32   34   35      0.000 0.0 1   5.390 180.0 2   1.230 0.0 3
torsion      33   32   34   35     -1.200 0.0 1   5.390 180.0 2   0.400 0.0 3


      #############################
      ##                         ##
      ##  Pi-Torsion Parameters  ##
      ##                         ##
      #############################


pitors        1    3            6.85
pitors        3   15            6.85
pitors       17   17            6.85
pitors       22   22            6.85
pitors       22   24            6.85
pitors       24   26            6.85
pitors       24   28            6.85
pitors       26   26            6.85


      ##################################
      ##                              ##
      ##  Torsion-Torsion Parameters  ##
      ##                              ##
      ##################################


tortors       3    1    2    3    1        25   25
  -180.0   -180.0    -0.41030
  -165.0   -180.0    -0.67838
  -150.0   -180.0    -1.32284
  -135.0   -180.0    -2.35817
  -120.0   -180.0    -3.54199
  -105.0   -180.0    -4.74648
   -90.0   -180.0    -5.73614
   -75.0   -180.0    -6.35918
   -60.0   -180.0    -6.35328
   -45.0   -180.0    -5.25542
   -30.0   -180.0    -3.16646
   -15.0   -180.0    -1.09098
     0.0   -180.0     1.38787
    15.0   -180.0    -1.09098
    30.0   -180.0    -3.16646
    45.0   -180.0    -5.25542
    60.0   -180.0    -6.35328
    75.0   -180.0    -6.36063
    90.0   -180.0    -5.73614
   105.0   -180.0    -4.74648
   120.0   -180.0    -3.54199
   135.0   -180.0    -2.35817
   150.0   -180.0    -1.32284
   165.0   -180.0    -0.67848
   180.0   -180.0    -0.41030
  -180.0   -165.0    -0.41413
  -165.0   -165.0    -0.90052
  -150.0   -165.0    -1.68618
  -135.0   -165.0    -2.64483
  -120.0   -165.0    -3.73864
  -105.0   -165.0    -4.87747
   -90.0   -165.0    -5.89392
   -75.0   -165.0    -6.55928
   -60.0   -165.0    -6.49751
   -45.0   -165.0    -5.22067
   -30.0   -165.0    -2.94954
   -15.0   -165.0    -0.08376
     0.0   -165.0     0.32603
    15.0   -165.0    -1.02648
    30.0   -165.0    -3.08291
    45.0   -165.0    -4.96362
    60.0   -165.0    -5.89147
    75.0   -165.0    -5.87160
    90.0   -165.0    -5.26258
   105.0   -165.0    -4.24252
   120.0   -165.0    -3.01099
   135.0   -165.0    -1.80556
   150.0   -165.0    -0.85965
   165.0   -165.0    -0.38834
   180.0   -165.0    -0.41413
  -180.0   -150.0    -0.34547
  -165.0   -150.0    -0.91119
  -150.0   -150.0    -1.61128
  -135.0   -150.0    -2.39523
  -120.0   -150.0    -3.36374
  -105.0   -150.0    -4.51571
   -90.0   -150.0    -5.59097
   -75.0   -150.0    -6.28648
   -60.0   -150.0    -4.83329
   -45.0   -150.0    -4.62413
   -30.0   -150.0    -1.82424
   -15.0   -150.0     1.19531
     0.0   -150.0     0.65267
    15.0   -150.0    -0.72633
    30.0   -150.0    -2.63732
    45.0   -150.0    -4.32002
    60.0   -150.0    -5.24397
    75.0   -150.0    -5.11358
    90.0   -150.0    -4.47212
   105.0   -150.0    -3.43923
   120.0   -150.0    -2.21872
   135.0   -150.0    -1.09788
   150.0   -150.0    -0.34336
   165.0   -150.0    -0.09173
   180.0   -150.0    -0.34547
  -180.0   -135.0     0.11457
  -165.0   -135.0    -0.35767
  -150.0   -135.0    -0.86749
  -135.0   -135.0    -1.46100
  -120.0   -135.0    -2.41372
  -105.0   -135.0    -3.62454
   -90.0   -135.0    -4.73673
   -75.0   -135.0    -5.35985
   -60.0   -135.0    -4.95241
   -45.0   -135.0    -2.94717
   -30.0   -135.0     1.07212
   -15.0   -135.0     1.77725
     0.0   -135.0     1.29409
    15.0   -135.0    -0.05629
    30.0   -135.0    -1.78936
    45.0   -135.0    -3.31494
    60.0   -135.0    -4.15255
    75.0   -135.0    -4.14746
    90.0   -135.0    -3.48022
   105.0   -135.0    -2.57323
   120.0   -135.0    -1.38411
   135.0   -135.0    -0.26179
   150.0   -135.0     0.30336
   165.0   -135.0     0.39970
   180.0   -135.0     0.11457
  -180.0   -120.0     0.94012
  -165.0   -120.0     0.68619
  -150.0   -120.0     0.39856
  -135.0   -120.0    -0.12378
  -120.0   -120.0    -1.19823
  -105.0   -120.0    -2.58368
   -90.0   -120.0    -3.50628
   -75.0   -120.0    -3.94292
   -60.0   -120.0    -3.11264
   -45.0   -120.0    -0.39145
   -30.0   -120.0     2.74127
   -15.0   -120.0     2.72433
     0.0   -120.0     2.15560
    15.0   -120.0     1.05732
    30.0   -120.0    -0.33944
    45.0   -120.0    -1.83317
    60.0   -120.0    -2.90987
    75.0   -120.0    -3.08145
    90.0   -120.0    -2.45277
   105.0   -120.0    -1.58770
   120.0   -120.0    -0.45731
   135.0   -120.0     0.42439
   150.0   -120.0     0.94220
   165.0   -120.0     1.14392
   180.0   -120.0     0.94012
  -180.0   -105.0     1.95487
  -165.0   -105.0     1.89156
  -150.0   -105.0     1.74844
  -135.0   -105.0     1.09827
  -120.0   -105.0    -0.23826
  -105.0   -105.0    -1.44395
   -90.0   -105.0    -2.34446
   -75.0   -105.0    -2.46809
   -60.0   -105.0    -1.05326
   -45.0   -105.0     1.58629
   -30.0   -105.0     3.35587
   -15.0   -105.0     3.45811
     0.0   -105.0     3.26537
    15.0   -105.0     2.74774
    30.0   -105.0     1.68837
    45.0   -105.0     0.01966
    60.0   -105.0    -1.39335
    75.0   -105.0    -2.00163
    90.0   -105.0    -1.75770
   105.0   -105.0    -0.79043
   120.0   -105.0     0.25289
   135.0   -105.0     1.13229
   150.0   -105.0     1.65296
   165.0   -105.0     1.97933
   180.0   -105.0     1.95487
  -180.0    -90.0     2.74331
  -165.0    -90.0     2.74719
  -150.0    -90.0     2.65974
  -135.0    -90.0     1.59315
  -120.0    -90.0     0.44790
  -105.0    -90.0    -0.60009
   -90.0    -90.0    -1.24738
   -75.0    -90.0    -1.04633
   -60.0    -90.0     0.58385
   -45.0    -90.0     2.52552
   -30.0    -90.0     3.38653
   -15.0    -90.0     3.96283
     0.0    -90.0     4.49297
    15.0    -90.0     4.70534
    30.0    -90.0     3.46788
    45.0    -90.0     1.48393
    60.0    -90.0    -0.28428
    75.0    -90.0    -1.08330
    90.0    -90.0    -1.11458
   105.0    -90.0    -0.40231
   120.0    -90.0     0.62829
   135.0    -90.0     1.55961
   150.0    -90.0     2.27845
   165.0    -90.0     2.58188
   180.0    -90.0     2.74331
  -180.0    -75.0     2.77751
  -165.0    -75.0     2.93260
  -150.0    -75.0     2.55866
  -135.0    -75.0     1.67411
  -120.0    -75.0     0.59458
  -105.0    -75.0    -0.17181
   -90.0    -75.0    -0.48079
   -75.0    -75.0     0.13385
   -60.0    -75.0     0.58100
   -45.0    -75.0     1.59048
   -30.0    -75.0     2.74299
   -15.0    -75.0     4.08945
     0.0    -75.0     5.44049
    15.0    -75.0     5.94623
    30.0    -75.0     4.24483
    45.0    -75.0     1.98183
    60.0    -75.0     0.19866
    75.0    -75.0    -0.95978
    90.0    -75.0    -1.06248
   105.0    -75.0    -0.43657
   120.0    -75.0     0.70386
   135.0    -75.0     1.60144
   150.0    -75.0     2.18997
   165.0    -75.0     2.53115
   180.0    -75.0     2.77751
  -180.0    -60.0     1.75943
  -165.0    -60.0     1.79376
  -150.0    -60.0     1.47445
  -135.0    -60.0     0.86891
  -120.0    -60.0     0.19601
  -105.0    -60.0    -0.17426
   -90.0    -60.0    -1.15692
   -75.0    -60.0    -1.49008
   -60.0    -60.0    -1.04413
   -45.0    -60.0     0.12852
   -30.0    -60.0     1.88655
   -15.0    -60.0     3.95852
     0.0    -60.0     5.75134
    15.0    -60.0     5.31247
    30.0    -60.0     3.44060
    45.0    -60.0     1.34478
    60.0    -60.0    -0.72835
    75.0    -60.0    -1.64238
    90.0    -60.0    -1.77690
   105.0    -60.0    -1.08681
   120.0    -60.0    -0.01095
   135.0    -60.0     0.71371
   150.0    -60.0     1.18267
   165.0    -60.0     1.51071
   180.0    -60.0     1.75943
  -180.0    -45.0    -0.27385
  -165.0    -45.0    -0.26642
  -150.0    -45.0    -0.36210
  -135.0    -45.0    -0.41066
  -120.0    -45.0    -1.56700
  -105.0    -45.0    -2.58412
   -90.0    -45.0    -3.31258
   -75.0    -45.0    -3.48890
   -60.0    -45.0    -2.82064
   -45.0    -45.0    -1.15606
   -30.0    -45.0     1.17974
   -15.0    -45.0     3.56251
     0.0    -45.0     4.30355
    15.0    -45.0     3.60759
    30.0    -45.0     1.82415
    45.0    -45.0    -0.23679
    60.0    -45.0    -2.11507
    75.0    -45.0    -3.16780
    90.0    -45.0    -3.23377
   105.0    -45.0    -2.42000
   120.0    -45.0    -1.56505
   135.0    -45.0    -1.07442
   150.0    -45.0    -0.69448
   165.0    -45.0    -0.39460
   180.0    -45.0    -0.27385
  -180.0    -30.0    -2.58364
  -165.0    -30.0    -2.46195
  -150.0    -30.0    -2.03519
  -135.0    -30.0    -3.32323
  -120.0    -30.0    -4.15184
  -105.0    -30.0    -4.81605
   -90.0    -30.0    -5.26121
   -75.0    -30.0    -5.15227
   -60.0    -30.0    -4.04692
   -45.0    -30.0    -2.03450
   -30.0    -30.0     0.38805
   -15.0    -30.0     1.96687
     0.0    -30.0     2.31228
    15.0    -30.0     1.55486
    30.0    -30.0    -0.06654
    45.0    -30.0    -2.00755
    60.0    -30.0    -3.84363
    75.0    -30.0    -4.69688
    90.0    -30.0    -4.67162
   105.0    -30.0    -4.04096
   120.0    -30.0    -3.61554
   135.0    -30.0    -3.21036
   150.0    -30.0    -2.88548
   165.0    -30.0    -2.68532
   180.0    -30.0    -2.58364
  -180.0    -15.0    -4.27536
  -165.0    -15.0    -3.60137
  -150.0    -15.0    -4.68528
  -135.0    -15.0    -5.26619
  -120.0    -15.0    -5.75171
  -105.0    -15.0    -6.18609
   -90.0    -15.0    -6.40257
   -75.0    -15.0    -5.99165
   -60.0    -15.0    -4.68897
   -45.0    -15.0    -2.79424
   -30.0    -15.0    -0.88273
   -15.0    -15.0     0.42755
     0.0    -15.0     0.73235
    15.0    -15.0     0.05777
    30.0    -15.0    -1.33943
    45.0    -15.0    -3.22837
    60.0    -15.0    -4.85236
    75.0    -15.0    -5.72489
    90.0    -15.0    -5.75813
   105.0    -15.0    -5.59602
   120.0    -15.0    -5.30079
   135.0    -15.0    -4.95535
   150.0    -15.0    -4.67044
   165.0    -15.0    -4.46238
   180.0    -15.0    -4.27536
  -180.0      0.0    -4.22230
  -165.0      0.0    -5.04625
  -150.0      0.0    -5.41586
  -135.0      0.0    -5.74456
  -120.0      0.0    -6.10583
  -105.0      0.0    -6.42560
   -90.0      0.0    -6.50266
   -75.0      0.0    -6.06098
   -60.0      0.0    -5.05604
   -45.0      0.0    -3.40770
   -30.0      0.0    -1.55029
   -15.0      0.0    -0.29366
     0.0      0.0     0.07505
    15.0      0.0    -0.34966
    30.0      0.0    -1.55069
    45.0      0.0    -3.40780
    60.0      0.0    -5.05604
    75.0      0.0    -6.06088
    90.0      0.0    -6.50276
   105.0      0.0    -6.42560
   120.0      0.0    -6.10583
   135.0      0.0    -5.74466
   150.0      0.0    -5.41586
   165.0      0.0    -5.04645
   180.0      0.0    -4.22230
  -180.0     15.0    -4.27526
  -165.0     15.0    -4.46238
  -150.0     15.0    -4.67044
  -135.0     15.0    -4.95535
  -120.0     15.0    -5.30069
  -105.0     15.0    -5.59602
   -90.0     15.0    -5.75823
   -75.0     15.0    -5.72479
   -60.0     15.0    -4.85236
   -45.0     15.0    -3.22837
   -30.0     15.0    -1.33903
   -15.0     15.0     0.05777
     0.0     15.0     0.73245
    15.0     15.0     0.42655
    30.0     15.0    -0.88283
    45.0     15.0    -2.79424
    60.0     15.0    -4.68897
    75.0     15.0    -5.99165
    90.0     15.0    -6.40257
   105.0     15.0    -6.18609
   120.0     15.0    -5.75161
   135.0     15.0    -5.26619
   150.0     15.0    -4.68528
   165.0     15.0    -3.60147
   180.0     15.0    -4.27526
  -180.0     30.0    -2.58364
  -165.0     30.0    -2.68532
  -150.0     30.0    -2.88548
  -135.0     30.0    -3.21036
  -120.0     30.0    -3.61554
  -105.0     30.0    -4.04096
   -90.0     30.0    -4.67162
   -75.0     30.0    -4.69688
   -60.0     30.0    -3.84373
   -45.0     30.0    -2.00765
   -30.0     30.0    -0.06654
   -15.0     30.0     1.55486
     0.0     30.0     2.31188
    15.0     30.0     1.95817
    30.0     30.0     0.38855
    45.0     30.0    -2.03450
    60.0     30.0    -4.04682
    75.0     30.0    -5.15227
    90.0     30.0    -5.26121
   105.0     30.0    -4.81605
   120.0     30.0    -4.15184
   135.0     30.0    -3.32353
   150.0     30.0    -2.03519
   165.0     30.0    -2.46185
   180.0     30.0    -2.58364
  -180.0     45.0    -0.27385
  -165.0     45.0    -0.39460
  -150.0     45.0    -0.69448
  -135.0     45.0    -1.07442
  -120.0     45.0    -1.56505
  -105.0     45.0    -2.42000
   -90.0     45.0    -3.23377
   -75.0     45.0    -3.16790
   -60.0     45.0    -2.11507
   -45.0     45.0    -0.23679
   -30.0     45.0     1.82415
   -15.0     45.0     3.60759
     0.0     45.0     4.30295
    15.0     45.0     3.56251
    30.0     45.0     1.17974
    45.0     45.0    -1.15606
    60.0     45.0    -2.82064
    75.0     45.0    -3.48890
    90.0     45.0    -3.31258
   105.0     45.0    -2.58412
   120.0     45.0    -1.56700
   135.0     45.0    -0.41066
   150.0     45.0    -0.36210
   165.0     45.0    -0.26642
   180.0     45.0    -0.27385
  -180.0     60.0     1.75943
  -165.0     60.0     1.51071
  -150.0     60.0     1.18267
  -135.0     60.0     0.71371
  -120.0     60.0    -0.01095
  -105.0     60.0    -1.08681
   -90.0     60.0    -1.77710
   -75.0     60.0    -1.64228
   -60.0     60.0    -0.72835
   -45.0     60.0     1.34478
   -30.0     60.0     3.44050
   -15.0     60.0     5.31247
     0.0     60.0     5.75134
    15.0     60.0     3.95862
    30.0     60.0     1.88655
    45.0     60.0     0.12792
    60.0     60.0    -1.04403
    75.0     60.0    -1.49008
    90.0     60.0    -1.15692
   105.0     60.0    -0.17426
   120.0     60.0     0.19601
   135.0     60.0     0.86911
   150.0     60.0     1.47445
   165.0     60.0     1.79366
   180.0     60.0     1.75943
  -180.0     75.0     2.77761
  -165.0     75.0     2.53125
  -150.0     75.0     2.18987
  -135.0     75.0     1.60144
  -120.0     75.0     0.70386
  -105.0     75.0    -0.43657
   -90.0     75.0    -1.06238
   -75.0     75.0    -0.95978
   -60.0     75.0     0.19866
   -45.0     75.0     1.98173
   -30.0     75.0     4.24483
   -15.0     75.0     5.94623
     0.0     75.0     5.44049
    15.0     75.0     4.08935
    30.0     75.0     2.74299
    45.0     75.0     1.59048
    60.0     75.0     0.58100
    75.0     75.0     0.13385
    90.0     75.0    -0.48069
   105.0     75.0    -0.17171
   120.0     75.0     0.59458
   135.0     75.0     1.67401
   150.0     75.0     2.55876
   165.0     75.0     2.93260
   180.0     75.0     2.77761
  -180.0     90.0     2.74331
  -165.0     90.0     2.58198
  -150.0     90.0     2.27845
  -135.0     90.0     1.55961
  -120.0     90.0     0.62829
  -105.0     90.0    -0.40231
   -90.0     90.0    -1.11458
   -75.0     90.0    -1.08320
   -60.0     90.0    -0.28428
   -45.0     90.0     1.48393
   -30.0     90.0     3.46788
   -15.0     90.0     4.70554
     0.0     90.0     4.49327
    15.0     90.0     3.96283
    30.0     90.0     3.38653
    45.0     90.0     2.52482
    60.0     90.0     0.58385
    75.0     90.0    -1.04633
    90.0     90.0    -1.24738
   105.0     90.0    -0.60009
   120.0     90.0     0.44790
   135.0     90.0     1.59315
   150.0     90.0     2.65974
   165.0     90.0     2.74719
   180.0     90.0     2.74331
  -180.0    105.0     1.95497
  -165.0    105.0     1.97933
  -150.0    105.0     1.65296
  -135.0    105.0     1.13219
  -120.0    105.0     0.25289
  -105.0    105.0    -0.79043
   -90.0    105.0    -1.75760
   -75.0    105.0    -2.00163
   -60.0    105.0    -1.39335
   -45.0    105.0     0.01966
   -30.0    105.0     1.68837
   -15.0    105.0     2.74754
     0.0    105.0     3.26537
    15.0    105.0     3.45811
    30.0    105.0     3.35587
    45.0    105.0     1.58629
    60.0    105.0    -1.05326
    75.0    105.0    -2.46809
    90.0    105.0    -2.34446
   105.0    105.0    -1.44395
   120.0    105.0    -0.23836
   135.0    105.0     1.09817
   150.0    105.0     1.74844
   165.0    105.0     1.89146
   180.0    105.0     1.95497
  -180.0    120.0     0.94002
  -165.0    120.0     1.14392
  -150.0    120.0     0.94220
  -135.0    120.0     0.42439
  -120.0    120.0    -0.45731
  -105.0    120.0    -1.58770
   -90.0    120.0    -2.45287
   -75.0    120.0    -3.08135
   -60.0    120.0    -2.90997
   -45.0    120.0    -1.83307
   -30.0    120.0    -0.33944
   -15.0    120.0     1.05732
     0.0    120.0     2.15560
    15.0    120.0     2.72433
    30.0    120.0     2.73567
    45.0    120.0    -0.51565
    60.0    120.0    -3.11264
    75.0    120.0    -3.94292
    90.0    120.0    -3.50628
   105.0    120.0    -2.58378
   120.0    120.0    -1.19823
   135.0    120.0    -0.12378
   150.0    120.0     0.39856
   165.0    120.0     0.68619
   180.0    120.0     0.94002
  -180.0    135.0     0.11457
  -165.0    135.0     0.39970
  -150.0    135.0     0.30326
  -135.0    135.0    -0.26169
  -120.0    135.0    -1.38411
  -105.0    135.0    -2.57323
   -90.0    135.0    -3.48022
   -75.0    135.0    -4.14746
   -60.0    135.0    -4.15255
   -45.0    135.0    -3.31494
   -30.0    135.0    -1.78936
   -15.0    135.0    -0.05629
     0.0    135.0     1.29409
    15.0    135.0     1.77725
    30.0    135.0     1.07212
    45.0    135.0    -2.94717
    60.0    135.0    -4.95241
    75.0    135.0    -5.35985
    90.0    135.0    -4.73673
   105.0    135.0    -3.62454
   120.0    135.0    -2.41382
   135.0    135.0    -1.46090
   150.0    135.0    -0.86739
   165.0    135.0    -0.35767
   180.0    135.0     0.11457
  -180.0    150.0    -0.34547
  -165.0    150.0    -0.09173
  -150.0    150.0    -0.34336
  -135.0    150.0    -1.09788
  -120.0    150.0    -2.21872
  -105.0    150.0    -3.43923
   -90.0    150.0    -4.47202
   -75.0    150.0    -5.11368
   -60.0    150.0    -5.15707
   -45.0    150.0    -4.32002
   -30.0    150.0    -2.63712
   -15.0    150.0    -0.72643
     0.0    150.0     0.43227
    15.0    150.0     1.19521
    30.0    150.0    -1.82424
    45.0    150.0    -4.62413
    60.0    150.0    -6.11149
    75.0    150.0    -6.28648
    90.0    150.0    -5.59097
   105.0    150.0    -4.51581
   120.0    150.0    -3.36364
   135.0    150.0    -2.39523
   150.0    150.0    -1.61138
   165.0    150.0    -0.91119
   180.0    150.0    -0.34547
  -180.0    165.0    -0.41413
  -165.0    165.0    -0.38834
  -150.0    165.0    -0.85965
  -135.0    165.0    -1.80556
  -120.0    165.0    -3.01099
  -105.0    165.0    -4.24252
   -90.0    165.0    -5.26258
   -75.0    165.0    -5.87160
   -60.0    165.0    -5.89147
   -45.0    165.0    -4.96362
   -30.0    165.0    -3.08291
   -15.0    165.0    -1.02648
     0.0    165.0     0.32603
    15.0    165.0    -0.08376
    30.0    165.0    -2.94954
    45.0    165.0    -5.22077
    60.0    165.0    -6.49751
    75.0    165.0    -6.55928
    90.0    165.0    -5.89392
   105.0    165.0    -4.87737
   120.0    165.0    -3.73864
   135.0    165.0    -2.64493
   150.0    165.0    -1.68618
   165.0    165.0    -0.90052
   180.0    165.0    -0.41413
  -180.0    180.0    -0.41030
  -165.0    180.0    -0.67838
  -150.0    180.0    -1.32284
  -135.0    180.0    -2.35817
  -120.0    180.0    -3.54199
  -105.0    180.0    -4.74648
   -90.0    180.0    -5.73614
   -75.0    180.0    -6.35918
   -60.0    180.0    -6.35328
   -45.0    180.0    -5.25542
   -30.0    180.0    -3.16646
   -15.0    180.0    -1.09098
     0.0    180.0     1.38787
    15.0    180.0    -1.09098
    30.0    180.0    -3.16646
    45.0    180.0    -5.25542
    60.0    180.0    -6.35328
    75.0    180.0    -6.36063
    90.0    180.0    -5.73614
   105.0    180.0    -4.74648
   120.0    180.0    -3.54199
   135.0    180.0    -2.35817
   150.0    180.0    -1.32284
   165.0    180.0    -0.67848
   180.0    180.0    -0.41030

tortors       3    1    2   32   34        25   25
  -180.0   -180.0    -0.41030
  -165.0   -180.0    -0.67838
  -150.0   -180.0    -1.32284
  -135.0   -180.0    -2.35817
  -120.0   -180.0    -3.54199
  -105.0   -180.0    -4.74648
   -90.0   -180.0    -5.73614
   -75.0   -180.0    -6.35918
   -60.0   -180.0    -6.35328
   -45.0   -180.0    -5.25542
   -30.0   -180.0    -3.16646
   -15.0   -180.0    -1.09098
     0.0   -180.0     1.38787
    15.0   -180.0    -1.09098
    30.0   -180.0    -3.16646
    45.0   -180.0    -5.25542
    60.0   -180.0    -6.35328
    75.0   -180.0    -6.36063
    90.0   -180.0    -5.73614
   105.0   -180.0    -4.74648
   120.0   -180.0    -3.54199
   135.0   -180.0    -2.35817
   150.0   -180.0    -1.32284
   165.0   -180.0    -0.67848
   180.0   -180.0    -0.41030
  -180.0   -165.0    -0.41413
  -165.0   -165.0    -0.90052
  -150.0   -165.0    -1.68618
  -135.0   -165.0    -2.64483
  -120.0   -165.0    -3.73864
  -105.0   -165.0    -4.87747
   -90.0   -165.0    -5.89392
   -75.0   -165.0    -6.55928
   -60.0   -165.0    -6.49751
   -45.0   -165.0    -5.22067
   -30.0   -165.0    -2.94954
   -15.0   -165.0    -0.08376
     0.0   -165.0     0.32603
    15.0   -165.0    -1.02648
    30.0   -165.0    -3.08291
    45.0   -165.0    -4.96362
    60.0   -165.0    -5.89147
    75.0   -165.0    -5.87160
    90.0   -165.0    -5.26258
   105.0   -165.0    -4.24252
   120.0   -165.0    -3.01099
   135.0   -165.0    -1.80556
   150.0   -165.0    -0.85965
   165.0   -165.0    -0.38834
   180.0   -165.0    -0.41413
  -180.0   -150.0    -0.34547
  -165.0   -150.0    -0.91119
  -150.0   -150.0    -1.61128
  -135.0   -150.0    -2.39523
  -120.0   -150.0    -3.36374
  -105.0   -150.0    -4.51571
   -90.0   -150.0    -5.59097
   -75.0   -150.0    -6.28648
   -60.0   -150.0    -4.83329
   -45.0   -150.0    -4.62413
   -30.0   -150.0    -1.82424
   -15.0   -150.0     1.19531
     0.0   -150.0     0.65267
    15.0   -150.0    -0.72633
    30.0   -150.0    -2.63732
    45.0   -150.0    -4.32002
    60.0   -150.0    -5.24397
    75.0   -150.0    -5.11358
    90.0   -150.0    -4.47212
   105.0   -150.0    -3.43923
   120.0   -150.0    -2.21872
   135.0   -150.0    -1.09788
   150.0   -150.0    -0.34336
   165.0   -150.0    -0.09173
   180.0   -150.0    -0.34547
  -180.0   -135.0     0.11457
  -165.0   -135.0    -0.35767
  -150.0   -135.0    -0.86749
  -135.0   -135.0    -1.46100
  -120.0   -135.0    -2.41372
  -105.0   -135.0    -3.62454
   -90.0   -135.0    -4.73673
   -75.0   -135.0    -5.35985
   -60.0   -135.0    -4.95241
   -45.0   -135.0    -2.94717
   -30.0   -135.0     1.07212
   -15.0   -135.0     1.77725
     0.0   -135.0     1.29409
    15.0   -135.0    -0.05629
    30.0   -135.0    -1.78936
    45.0   -135.0    -3.31494
    60.0   -135.0    -4.15255
    75.0   -135.0    -4.14746
    90.0   -135.0    -3.48022
   105.0   -135.0    -2.57323
   120.0   -135.0    -1.38411
   135.0   -135.0    -0.26179
   150.0   -135.0     0.30336
   165.0   -135.0     0.39970
   180.0   -135.0     0.11457
  -180.0   -120.0     0.94012
  -165.0   -120.0     0.68619
  -150.0   -120.0     0.39856
  -135.0   -120.0    -0.12378
  -120.0   -120.0    -1.19823
  -105.0   -120.0    -2.58368
   -90.0   -120.0    -3.50628
   -75.0   -120.0    -3.94292
   -60.0   -120.0    -3.11264
   -45.0   -120.0    -0.39145
   -30.0   -120.0     2.74127
   -15.0   -120.0     2.72433
     0.0   -120.0     2.15560
    15.0   -120.0     1.05732
    30.0   -120.0    -0.33944
    45.0   -120.0    -1.83317
    60.0   -120.0    -2.90987
    75.0   -120.0    -3.08145
    90.0   -120.0    -2.45277
   105.0   -120.0    -1.58770
   120.0   -120.0    -0.45731
   135.0   -120.0     0.42439
   150.0   -120.0     0.94220
   165.0   -120.0     1.14392
   180.0   -120.0     0.94012
  -180.0   -105.0     1.95487
  -165.0   -105.0     1.89156
  -150.0   -105.0     1.74844
  -135.0   -105.0     1.09827
  -120.0   -105.0    -0.23826
  -105.0   -105.0    -1.44395
   -90.0   -105.0    -2.34446
   -75.0   -105.0    -2.46809
   -60.0   -105.0    -1.05326
   -45.0   -105.0     1.58629
   -30.0   -105.0     3.35587
   -15.0   -105.0     3.45811
     0.0   -105.0     3.26537
    15.0   -105.0     2.74774
    30.0   -105.0     1.68837
    45.0   -105.0     0.01966
    60.0   -105.0    -1.39335
    75.0   -105.0    -2.00163
    90.0   -105.0    -1.75770
   105.0   -105.0    -0.79043
   120.0   -105.0     0.25289
   135.0   -105.0     1.13229
   150.0   -105.0     1.65296
   165.0   -105.0     1.97933
   180.0   -105.0     1.95487
  -180.0    -90.0     2.74331
  -165.0    -90.0     2.74719
  -150.0    -90.0     2.65974
  -135.0    -90.0     1.59315
  -120.0    -90.0     0.44790
  -105.0    -90.0    -0.60009
   -90.0    -90.0    -1.24738
   -75.0    -90.0    -1.04633
   -60.0    -90.0     0.58385
   -45.0    -90.0     2.52552
   -30.0    -90.0     3.38653
   -15.0    -90.0     3.96283
     0.0    -90.0     4.49297
    15.0    -90.0     4.70534
    30.0    -90.0     3.46788
    45.0    -90.0     1.48393
    60.0    -90.0    -0.28428
    75.0    -90.0    -1.08330
    90.0    -90.0    -1.11458
   105.0    -90.0    -0.40231
   120.0    -90.0     0.62829
   135.0    -90.0     1.55961
   150.0    -90.0     2.27845
   165.0    -90.0     2.58188
   180.0    -90.0     2.74331
  -180.0    -75.0     2.77751
  -165.0    -75.0     2.93260
  -150.0    -75.0     2.55866
  -135.0    -75.0     1.67411
  -120.0    -75.0     0.59458
  -105.0    -75.0    -0.17181
   -90.0    -75.0    -0.48079
   -75.0    -75.0     0.13385
   -60.0    -75.0     0.58100
   -45.0    -75.0     1.59048
   -30.0    -75.0     2.74299
   -15.0    -75.0     4.08945
     0.0    -75.0     5.44049
    15.0    -75.0     5.94623
    30.0    -75.0     4.24483
    45.0    -75.0     1.98183
    60.0    -75.0     0.19866
    75.0    -75.0    -0.95978
    90.0    -75.0    -1.06248
   105.0    -75.0    -0.43657
   120.0    -75.0     0.70386
   135.0    -75.0     1.60144
   150.0    -75.0     2.18997
   165.0    -75.0     2.53115
   180.0    -75.0     2.77751
  -180.0    -60.0     1.75943
  -165.0    -60.0     1.79376
  -150.0    -60.0     1.47445
  -135.0    -60.0     0.86891
  -120.0    -60.0     0.19601
  -105.0    -60.0    -0.17426
   -90.0    -60.0    -1.15692
   -75.0    -60.0    -1.49008
   -60.0    -60.0    -1.04413
   -45.0    -60.0     0.12852
   -30.0    -60.0     1.88655
   -15.0    -60.0     3.95852
     0.0    -60.0     5.75134
    15.0    -60.0     5.31247
    30.0    -60.0     3.44060
    45.0    -60.0     1.34478
    60.0    -60.0    -0.72835
    75.0    -60.0    -1.64238
    90.0    -60.0    -1.77690
   105.0    -60.0    -1.08681
   120.0    -60.0    -0.01095
   135.0    -60.0     0.71371
   150.0    -60.0     1.18267
   165.0    -60.0     1.51071
   180.0    -60.0     1.75943
  -180.0    -45.0    -0.27385
  -165.0    -45.0    -0.26642
  -150.0    -45.0    -0.36210
  -135.0    -45.0    -0.41066
  -120.0    -45.0    -1.56700
  -105.0    -45.0    -2.58412
   -90.0    -45.0    -3.31258
   -75.0    -45.0    -3.48890
   -60.0    -45.0    -2.82064
   -45.0    -45.0    -1.15606
   -30.0    -45.0     1.17974
   -15.0    -45.0     3.56251
     0.0    -45.0     4.30355
    15.0    -45.0     3.60759
    30.0    -45.0     1.82415
    45.0    -45.0    -0.23679
    60.0    -45.0    -2.11507
    75.0    -45.0    -3.16780
    90.0    -45.0    -3.23377
   105.0    -45.0    -2.42000
   120.0    -45.0    -1.56505
   135.0    -45.0    -1.07442
   150.0    -45.0    -0.69448
   165.0    -45.0    -0.39460
   180.0    -45.0    -0.27385
  -180.0    -30.0    -2.58364
  -165.0    -30.0    -2.46195
  -150.0    -30.0    -2.03519
  -135.0    -30.0    -3.32323
  -120.0    -30.0    -4.15184
  -105.0    -30.0    -4.81605
   -90.0    -30.0    -5.26121
   -75.0    -30.0    -5.15227
   -60.0    -30.0    -4.04692
   -45.0    -30.0    -2.03450
   -30.0    -30.0     0.38805
   -15.0    -30.0     1.96687
     0.0    -30.0     2.31228
    15.0    -30.0     1.55486
    30.0    -30.0    -0.06654
    45.0    -30.0    -2.00755
    60.0    -30.0    -3.84363
    75.0    -30.0    -4.69688
    90.0    -30.0    -4.67162
   105.0    -30.0    -4.04096
   120.0    -30.0    -3.61554
   135.0    -30.0    -3.21036
   150.0    -30.0    -2.88548
   165.0    -30.0    -2.68532
   180.0    -30.0    -2.58364
  -180.0    -15.0    -4.27536
  -165.0    -15.0    -3.60137
  -150.0    -15.0    -4.68528
  -135.0    -15.0    -5.26619
  -120.0    -15.0    -5.75171
  -105.0    -15.0    -6.18609
   -90.0    -15.0    -6.40257
   -75.0    -15.0    -5.99165
   -60.0    -15.0    -4.68897
   -45.0    -15.0    -2.79424
   -30.0    -15.0    -0.88273
   -15.0    -15.0     0.42755
     0.0    -15.0     0.73235
    15.0    -15.0     0.05777
    30.0    -15.0    -1.33943
    45.0    -15.0    -3.22837
    60.0    -15.0    -4.85236
    75.0    -15.0    -5.72489
    90.0    -15.0    -5.75813
   105.0    -15.0    -5.59602
   120.0    -15.0    -5.30079
   135.0    -15.0    -4.95535
   150.0    -15.0    -4.67044
   165.0    -15.0    -4.46238
   180.0    -15.0    -4.27536
  -180.0      0.0    -4.22230
  -165.0      0.0    -5.04625
  -150.0      0.0    -5.41586
  -135.0      0.0    -5.74456
  -120.0      0.0    -6.10583
  -105.0      0.0    -6.42560
   -90.0      0.0    -6.50266
   -75.0      0.0    -6.06098
   -60.0      0.0    -5.05604
   -45.0      0.0    -3.40770
   -30.0      0.0    -1.55029
   -15.0      0.0    -0.29366
     0.0      0.0     0.07505
    15.0      0.0    -0.34966
    30.0      0.0    -1.55069
    45.0      0.0    -3.40780
    60.0      0.0    -5.05604
    75.0      0.0    -6.06088
    90.0      0.0    -6.50276
   105.0      0.0    -6.42560
   120.0      0.0    -6.10583
   135.0      0.0    -5.74466
   150.0      0.0    -5.41586
   165.0      0.0    -5.04645
   180.0      0.0    -4.22230
  -180.0     15.0    -4.27526
  -165.0     15.0    -4.46238
  -150.0     15.0    -4.67044
  -135.0     15.0    -4.95535
  -120.0     15.0    -5.30069
  -105.0     15.0    -5.59602
   -90.0     15.0    -5.75823
   -75.0     15.0    -5.72479
   -60.0     15.0    -4.85236
   -45.0     15.0    -3.22837
   -30.0     15.0    -1.33903
   -15.0     15.0     0.05777
     0.0     15.0     0.73245
    15.0     15.0     0.42655
    30.0     15.0    -0.88283
    45.0     15.0    -2.79424
    60.0     15.0    -4.68897
    75.0     15.0    -5.99165
    90.0     15.0    -6.40257
   105.0     15.0    -6.18609
   120.0     15.0    -5.75161
   135.0     15.0    -5.26619
   150.0     15.0    -4.68528
   165.0     15.0    -3.60147
   180.0     15.0    -4.27526
  -180.0     30.0    -2.58364
  -165.0     30.0    -2.68532
  -150.0     30.0    -2.88548
  -135.0     30.0    -3.21036
  -120.0     30.0    -3.61554
  -105.0     30.0    -4.04096
   -90.0     30.0    -4.67162
   -75.0     30.0    -4.69688
   -60.0     30.0    -3.84373
   -45.0     30.0    -2.00765
   -30.0     30.0    -0.06654
   -15.0     30.0     1.55486
     0.0     30.0     2.31188
    15.0     30.0     1.95817
    30.0     30.0     0.38855
    45.0     30.0    -2.03450
    60.0     30.0    -4.04682
    75.0     30.0    -5.15227
    90.0     30.0    -5.26121
   105.0     30.0    -4.81605
   120.0     30.0    -4.15184
   135.0     30.0    -3.32353
   150.0     30.0    -2.03519
   165.0     30.0    -2.46185
   180.0     30.0    -2.58364
  -180.0     45.0    -0.27385
  -165.0     45.0    -0.39460
  -150.0     45.0    -0.69448
  -135.0     45.0    -1.07442
  -120.0     45.0    -1.56505
  -105.0     45.0    -2.42000
   -90.0     45.0    -3.23377
   -75.0     45.0    -3.16790
   -60.0     45.0    -2.11507
   -45.0     45.0    -0.23679
   -30.0     45.0     1.82415
   -15.0     45.0     3.60759
     0.0     45.0     4.30295
    15.0     45.0     3.56251
    30.0     45.0     1.17974
    45.0     45.0    -1.15606
    60.0     45.0    -2.82064
    75.0     45.0    -3.48890
    90.0     45.0    -3.31258
   105.0     45.0    -2.58412
   120.0     45.0    -1.56700
   135.0     45.0    -0.41066
   150.0     45.0    -0.36210
   165.0     45.0    -0.26642
   180.0     45.0    -0.27385
  -180.0     60.0     1.75943
  -165.0     60.0     1.51071
  -150.0     60.0     1.18267
  -135.0     60.0     0.71371
  -120.0     60.0    -0.01095
  -105.0     60.0    -1.08681
   -90.0     60.0    -1.77710
   -75.0     60.0    -1.64228
   -60.0     60.0    -0.72835
   -45.0     60.0     1.34478
   -30.0     60.0     3.44050
   -15.0     60.0     5.31247
     0.0     60.0     5.75134
    15.0     60.0     3.95862
    30.0     60.0     1.88655
    45.0     60.0     0.12792
    60.0     60.0    -1.04403
    75.0     60.0    -1.49008
    90.0     60.0    -1.15692
   105.0     60.0    -0.17426
   120.0     60.0     0.19601
   135.0     60.0     0.86911
   150.0     60.0     1.47445
   165.0     60.0     1.79366
   180.0     60.0     1.75943
  -180.0     75.0     2.77761
  -165.0     75.0     2.53125
  -150.0     75.0     2.18987
  -135.0     75.0     1.60144
  -120.0     75.0     0.70386
  -105.0     75.0    -0.43657
   -90.0     75.0    -1.06238
   -75.0     75.0    -0.95978
   -60.0     75.0     0.19866
   -45.0     75.0     1.98173
   -30.0     75.0     4.24483
   -15.0     75.0     5.94623
     0.0     75.0     5.44049
    15.0     75.0     4.08935
    30.0     75.0     2.74299
    45.0     75.0     1.59048
    60.0     75.0     0.58100
    75.0     75.0     0.13385
    90.0     75.0    -0.48069
   105.0     75.0    -0.17171
   120.0     75.0     0.59458
   135.0     75.0     1.67401
   150.0     75.0     2.55876
   165.0     75.0     2.93260
   180.0     75.0     2.77761
  -180.0     90.0     2.74331
  -165.0     90.0     2.58198
  -150.0     90.0     2.27845
  -135.0     90.0     1.55961
  -120.0     90.0     0.62829
  -105.0     90.0    -0.40231
   -90.0     90.0    -1.11458
   -75.0     90.0    -1.08320
   -60.0     90.0    -0.28428
   -45.0     90.0     1.48393
   -30.0     90.0     3.46788
   -15.0     90.0     4.70554
     0.0     90.0     4.49327
    15.0     90.0     3.96283
    30.0     90.0     3.38653
    45.0     90.0     2.52482
    60.0     90.0     0.58385
    75.0     90.0    -1.04633
    90.0     90.0    -1.24738
   105.0     90.0    -0.60009
   120.0     90.0     0.44790
   135.0     90.0     1.59315
   150.0     90.0     2.65974
   165.0     90.0     2.74719
   180.0     90.0     2.74331
  -180.0    105.0     1.95497
  -165.0    105.0     1.97933
  -150.0    105.0     1.65296
  -135.0    105.0     1.13219
  -120.0    105.0     0.25289
  -105.0    105.0    -0.79043
   -90.0    105.0    -1.75760
   -75.0    105.0    -2.00163
   -60.0    105.0    -1.39335
   -45.0    105.0     0.01966
   -30.0    105.0     1.68837
   -15.0    105.0     2.74754
     0.0    105.0     3.26537
    15.0    105.0     3.45811
    30.0    105.0     3.35587
    45.0    105.0     1.58629
    60.0    105.0    -1.05326
    75.0    105.0    -2.46809
    90.0    105.0    -2.34446
   105.0    105.0    -1.44395
   120.0    105.0    -0.23836
   135.0    105.0     1.09817
   150.0    105.0     1.74844
   165.0    105.0     1.89146
   180.0    105.0     1.95497
  -180.0    120.0     0.94002
  -165.0    120.0     1.14392
  -150.0    120.0     0.94220
  -135.0    120.0     0.42439
  -120.0    120.0    -0.45731
  -105.0    120.0    -1.58770
   -90.0    120.0    -2.45287
   -75.0    120.0    -3.08135
   -60.0    120.0    -2.90997
   -45.0    120.0    -1.83307
   -30.0    120.0    -0.33944
   -15.0    120.0     1.05732
     0.0    120.0     2.15560
    15.0    120.0     2.72433
    30.0    120.0     2.73567
    45.0    120.0    -0.51565
    60.0    120.0    -3.11264
    75.0    120.0    -3.94292
    90.0    120.0    -3.50628
   105.0    120.0    -2.58378
   120.0    120.0    -1.19823
   135.0    120.0    -0.12378
   150.0    120.0     0.39856
   165.0    120.0     0.68619
   180.0    120.0     0.94002
  -180.0    135.0     0.11457
  -165.0    135.0     0.39970
  -150.0    135.0     0.30326
  -135.0    135.0    -0.26169
  -120.0    135.0    -1.38411
  -105.0    135.0    -2.57323
   -90.0    135.0    -3.48022
   -75.0    135.0    -4.14746
   -60.0    135.0    -4.15255
   -45.0    135.0    -3.31494
   -30.0    135.0    -1.78936
   -15.0    135.0    -0.05629
     0.0    135.0     1.29409
    15.0    135.0     1.77725
    30.0    135.0     1.07212
    45.0    135.0    -2.94717
    60.0    135.0    -4.95241
    75.0    135.0    -5.35985
    90.0    135.0    -4.73673
   105.0    135.0    -3.62454
   120.0    135.0    -2.41382
   135.0    135.0    -1.46090
   150.0    135.0    -0.86739
   165.0    135.0    -0.35767
   180.0    135.0     0.11457
  -180.0    150.0    -0.34547
  -165.0    150.0    -0.09173
  -150.0    150.0    -0.34336
  -135.0    150.0    -1.09788
  -120.0    150.0    -2.21872
  -105.0    150.0    -3.43923
   -90.0    150.0    -4.47202
   -75.0    150.0    -5.11368
   -60.0    150.0    -5.15707
   -45.0    150.0    -4.32002
   -30.0    150.0    -2.63712
   -15.0    150.0    -0.72643
     0.0    150.0     0.43227
    15.0    150.0     1.19521
    30.0    150.0    -1.82424
    45.0    150.0    -4.62413
    60.0    150.0    -6.11149
    75.0    150.0    -6.28648
    90.0    150.0    -5.59097
   105.0    150.0    -4.51581
   120.0    150.0    -3.36364
   135.0    150.0    -2.39523
   150.0    150.0    -1.61138
   165.0    150.0    -0.91119
   180.0    150.0    -0.34547
  -180.0    165.0    -0.41413
  -165.0    165.0    -0.38834
  -150.0    165.0    -0.85965
  -135.0    165.0    -1.80556
  -120.0    165.0    -3.01099
  -105.0    165.0    -4.24252
   -90.0    165.0    -5.26258
   -75.0    165.0    -5.87160
   -60.0    165.0    -5.89147
   -45.0    165.0    -4.96362
   -30.0    165.0    -3.08291
   -15.0    165.0    -1.02648
     0.0    165.0     0.32603
    15.0    165.0    -0.08376
    30.0    165.0    -2.94954
    45.0    165.0    -5.22077
    60.0    165.0    -6.49751
    75.0    165.0    -6.55928
    90.0    165.0    -5.89392
   105.0    165.0    -4.87737
   120.0    165.0    -3.73864
   135.0    165.0    -2.64493
   150.0    165.0    -1.68618
   165.0    165.0    -0.90052
   180.0    165.0    -0.41413
  -180.0    180.0    -0.41030
  -165.0    180.0    -0.67838
  -150.0    180.0    -1.32284
  -135.0    180.0    -2.35817
  -120.0    180.0    -3.54199
  -105.0    180.0    -4.74648
   -90.0    180.0    -5.73614
   -75.0    180.0    -6.35918
   -60.0    180.0    -6.35328
   -45.0    180.0    -5.25542
   -30.0    180.0    -3.16646
   -15.0    180.0    -1.09098
     0.0    180.0     1.38787
    15.0    180.0    -1.09098
    30.0    180.0    -3.16646
    45.0    180.0    -5.25542
    60.0    180.0    -6.35328
    75.0    180.0    -6.36063
    90.0    180.0    -5.73614
   105.0    180.0    -4.74648
   120.0    180.0    -3.54199
   135.0    180.0    -2.35817
   150.0    180.0    -1.32284
   165.0    180.0    -0.67848
   180.0    180.0    -0.41030

tortors       3    1    2    3   15        25   25
  -180.0   -180.0    -0.41030
  -165.0   -180.0    -0.67838
  -150.0   -180.0    -1.32284
  -135.0   -180.0    -2.35817
  -120.0   -180.0    -3.54199
  -105.0   -180.0    -4.74648
   -90.0   -180.0    -5.73614
   -75.0   -180.0    -6.35918
   -60.0   -180.0    -6.35328
   -45.0   -180.0    -5.25542
   -30.0   -180.0    -3.16646
   -15.0   -180.0    -1.09098
     0.0   -180.0     1.38787
    15.0   -180.0    -1.09098
    30.0   -180.0    -3.16646
    45.0   -180.0    -5.25542
    60.0   -180.0    -6.35328
    75.0   -180.0    -6.36063
    90.0   -180.0    -5.73614
   105.0   -180.0    -4.74648
   120.0   -180.0    -3.54199
   135.0   -180.0    -2.35817
   150.0   -180.0    -1.32284
   165.0   -180.0    -0.67848
   180.0   -180.0    -0.41030
  -180.0   -165.0    -0.41413
  -165.0   -165.0    -0.90052
  -150.0   -165.0    -1.68618
  -135.0   -165.0    -2.64483
  -120.0   -165.0    -3.73864
  -105.0   -165.0    -4.87747
   -90.0   -165.0    -5.89392
   -75.0   -165.0    -6.55928
   -60.0   -165.0    -6.49751
   -45.0   -165.0    -5.22067
   -30.0   -165.0    -2.94954
   -15.0   -165.0    -0.08376
     0.0   -165.0     0.32603
    15.0   -165.0    -1.02648
    30.0   -165.0    -3.08291
    45.0   -165.0    -4.96362
    60.0   -165.0    -5.89147
    75.0   -165.0    -5.87160
    90.0   -165.0    -5.26258
   105.0   -165.0    -4.24252
   120.0   -165.0    -3.01099
   135.0   -165.0    -1.80556
   150.0   -165.0    -0.85965
   165.0   -165.0    -0.38834
   180.0   -165.0    -0.41413
  -180.0   -150.0    -0.34547
  -165.0   -150.0    -0.91119
  -150.0   -150.0    -1.61128
  -135.0   -150.0    -2.39523
  -120.0   -150.0    -3.36374
  -105.0   -150.0    -4.51571
   -90.0   -150.0    -5.59097
   -75.0   -150.0    -6.28648
   -60.0   -150.0    -4.83329
   -45.0   -150.0    -4.62413
   -30.0   -150.0    -1.82424
   -15.0   -150.0     1.19531
     0.0   -150.0     0.65267
    15.0   -150.0    -0.72633
    30.0   -150.0    -2.63732
    45.0   -150.0    -4.32002
    60.0   -150.0    -5.24397
    75.0   -150.0    -5.11358
    90.0   -150.0    -4.47212
   105.0   -150.0    -3.43923
   120.0   -150.0    -2.21872
   135.0   -150.0    -1.09788
   150.0   -150.0    -0.34336
   165.0   -150.0    -0.09173
   180.0   -150.0    -0.34547
  -180.0   -135.0     0.11457
  -165.0   -135.0    -0.35767
  -150.0   -135.0    -0.86749
  -135.0   -135.0    -1.46100
  -120.0   -135.0    -2.41372
  -105.0   -135.0    -3.62454
   -90.0   -135.0    -4.73673
   -75.0   -135.0    -5.35985
   -60.0   -135.0    -4.95241
   -45.0   -135.0    -2.94717
   -30.0   -135.0     1.07212
   -15.0   -135.0     1.77725
     0.0   -135.0     1.29409
    15.0   -135.0    -0.05629
    30.0   -135.0    -1.78936
    45.0   -135.0    -3.31494
    60.0   -135.0    -4.15255
    75.0   -135.0    -4.14746
    90.0   -135.0    -3.48022
   105.0   -135.0    -2.57323
   120.0   -135.0    -1.38411
   135.0   -135.0    -0.26179
   150.0   -135.0     0.30336
   165.0   -135.0     0.39970
   180.0   -135.0     0.11457
  -180.0   -120.0     0.94012
  -165.0   -120.0     0.68619
  -150.0   -120.0     0.39856
  -135.0   -120.0    -0.12378
  -120.0   -120.0    -1.19823
  -105.0   -120.0    -2.58368
   -90.0   -120.0    -3.50628
   -75.0   -120.0    -3.94292
   -60.0   -120.0    -3.11264
   -45.0   -120.0    -0.39145
   -30.0   -120.0     2.74127
   -15.0   -120.0     2.72433
     0.0   -120.0     2.15560
    15.0   -120.0     1.05732
    30.0   -120.0    -0.33944
    45.0   -120.0    -1.83317
    60.0   -120.0    -2.90987
    75.0   -120.0    -3.08145
    90.0   -120.0    -2.45277
   105.0   -120.0    -1.58770
   120.0   -120.0    -0.45731
   135.0   -120.0     0.42439
   150.0   -120.0     0.94220
   165.0   -120.0     1.14392
   180.0   -120.0     0.94012
  -180.0   -105.0     1.95487
  -165.0   -105.0     1.89156
  -150.0   -105.0     1.74844
  -135.0   -105.0     1.09827
  -120.0   -105.0    -0.23826
  -105.0   -105.0    -1.44395
   -90.0   -105.0    -2.34446
   -75.0   -105.0    -2.46809
   -60.0   -105.0    -1.05326
   -45.0   -105.0     1.58629
   -30.0   -105.0     3.35587
   -15.0   -105.0     3.45811
     0.0   -105.0     3.26537
    15.0   -105.0     2.74774
    30.0   -105.0     1.68837
    45.0   -105.0     0.01966
    60.0   -105.0    -1.39335
    75.0   -105.0    -2.00163
    90.0   -105.0    -1.75770
   105.0   -105.0    -0.79043
   120.0   -105.0     0.25289
   135.0   -105.0     1.13229
   150.0   -105.0     1.65296
   165.0   -105.0     1.97933
   180.0   -105.0     1.95487
  -180.0    -90.0     2.74331
  -165.0    -90.0     2.74719
  -150.0    -90.0     2.65974
  -135.0    -90.0     1.59315
  -120.0    -90.0     0.44790
  -105.0    -90.0    -0.60009
   -90.0    -90.0    -1.24738
   -75.0    -90.0    -1.04633
   -60.0    -90.0     0.58385
   -45.0    -90.0     2.52552
   -30.0    -90.0     3.38653
   -15.0    -90.0     3.96283
     0.0    -90.0     4.49297
    15.0    -90.0     4.70534
    30.0    -90.0     3.46788
    45.0    -90.0     1.48393
    60.0    -90.0    -0.28428
    75.0    -90.0    -1.08330
    90.0    -90.0    -1.11458
   105.0    -90.0    -0.40231
   120.0    -90.0     0.62829
   135.0    -90.0     1.55961
   150.0    -90.0     2.27845
   165.0    -90.0     2.58188
   180.0    -90.0     2.74331
  -180.0    -75.0     2.77751
  -165.0    -75.0     2.93260
  -150.0    -75.0     2.55866
  -135.0    -75.0     1.67411
  -120.0    -75.0     0.59458
  -105.0    -75.0    -0.17181
   -90.0    -75.0    -0.48079
   -75.0    -75.0     0.13385
   -60.0    -75.0     0.58100
   -45.0    -75.0     1.59048
   -30.0    -75.0     2.74299
   -15.0    -75.0     4.08945
     0.0    -75.0     5.44049
    15.0    -75.0     5.94623
    30.0    -75.0     4.24483
    45.0    -75.0     1.98183
    60.0    -75.0     0.19866
    75.0    -75.0    -0.95978
    90.0    -75.0    -1.06248
   105.0    -75.0    -0.43657
   120.0    -75.0     0.70386
   135.0    -75.0     1.60144
   150.0    -75.0     2.18997
   165.0    -75.0     2.53115
   180.0    -75.0     2.77751
  -180.0    -60.0     1.75943
  -165.0    -60.0     1.79376
  -150.0    -60.0     1.47445
  -135.0    -60.0     0.86891
  -120.0    -60.0     0.19601
  -105.0    -60.0    -0.17426
   -90.0    -60.0    -1.15692
   -75.0    -60.0    -1.49008
   -60.0    -60.0    -1.04413
   -45.0    -60.0     0.12852
   -30.0    -60.0     1.88655
   -15.0    -60.0     3.95852
     0.0    -60.0     5.75134
    15.0    -60.0     5.31247
    30.0    -60.0     3.44060
    45.0    -60.0     1.34478
    60.0    -60.0    -0.72835
    75.0    -60.0    -1.64238
    90.0    -60.0    -1.77690
   105.0    -60.0    -1.08681
   120.0    -60.0    -0.01095
   135.0    -60.0     0.71371
   150.0    -60.0     1.18267
   165.0    -60.0     1.51071
   180.0    -60.0     1.75943
  -180.0    -45.0    -0.27385
  -165.0    -45.0    -0.26642
  -150.0    -45.0    -0.36210
  -135.0    -45.0    -0.41066
  -120.0    -45.0    -1.56700
  -105.0    -45.0    -2.58412
   -90.0    -45.0    -3.31258
   -75.0    -45.0    -3.48890
   -60.0    -45.0    -2.82064
   -45.0    -45.0    -1.15606
   -30.0    -45.0     1.17974
   -15.0    -45.0     3.56251
     0.0    -45.0     4.30355
    15.0    -45.0     3.60759
    30.0    -45.0     1.82415
    45.0    -45.0    -0.23679
    60.0    -45.0    -2.11507
    75.0    -45.0    -3.16780
    90.0    -45.0    -3.23377
   105.0    -45.0    -2.42000
   120.0    -45.0    -1.56505
   135.0    -45.0    -1.07442
   150.0    -45.0    -0.69448
   165.0    -45.0    -0.39460
   180.0    -45.0    -0.27385
  -180.0    -30.0    -2.58364
  -165.0    -30.0    -2.46195
  -150.0    -30.0    -2.03519
  -135.0    -30.0    -3.32323
  -120.0    -30.0    -4.15184
  -105.0    -30.0    -4.81605
   -90.0    -30.0    -5.26121
   -75.0    -30.0    -5.15227
   -60.0    -30.0    -4.04692
   -45.0    -30.0    -2.03450
   -30.0    -30.0     0.38805
   -15.0    -30.0     1.96687
     0.0    -30.0     2.31228
    15.0    -30.0     1.55486
    30.0    -30.0    -0.06654
    45.0    -30.0    -2.00755
    60.0    -30.0    -3.84363
    75.0    -30.0    -4.69688
    90.0    -30.0    -4.67162
   105.0    -30.0    -4.04096
   120.0    -30.0    -3.61554
   135.0    -30.0    -3.21036
   150.0    -30.0    -2.88548
   165.0    -30.0    -2.68532
   180.0    -30.0    -2.58364
  -180.0    -15.0    -4.27536
  -165.0    -15.0    -3.60137
  -150.0    -15.0    -4.68528
  -135.0    -15.0    -5.26619
  -120.0    -15.0    -5.75171
  -105.0    -15.0    -6.18609
   -90.0    -15.0    -6.40257
   -75.0    -15.0    -5.99165
   -60.0    -15.0    -4.68897
   -45.0    -15.0    -2.79424
   -30.0    -15.0    -0.88273
   -15.0    -15.0     0.42755
     0.0    -15.0     0.73235
    15.0    -15.0     0.05777
    30.0    -15.0    -1.33943
    45.0    -15.0    -3.22837
    60.0    -15.0    -4.85236
    75.0    -15.0    -5.72489
    90.0    -15.0    -5.75813
   105.0    -15.0    -5.59602
   120.0    -15.0    -5.30079
   135.0    -15.0    -4.95535
   150.0    -15.0    -4.67044
   165.0    -15.0    -4.46238
   180.0    -15.0    -4.27536
  -180.0      0.0    -4.22230
  -165.0      0.0    -5.04625
  -150.0      0.0    -5.41586
  -135.0      0.0    -5.74456
  -120.0      0.0    -6.10583
  -105.0      0.0    -6.42560
   -90.0      0.0    -6.50266
   -75.0      0.0    -6.06098
   -60.0      0.0    -5.05604
   -45.0      0.0    -3.40770
   -30.0      0.0    -1.55029
   -15.0      0.0    -0.29366
     0.0      0.0     0.07505
    15.0      0.0    -0.34966
    30.0      0.0    -1.55069
    45.0      0.0    -3.40780
    60.0      0.0    -5.05604
    75.0      0.0    -6.06088
    90.0      0.0    -6.50276
   105.0      0.0    -6.42560
   120.0      0.0    -6.10583
   135.0      0.0    -5.74466
   150.0      0.0    -5.41586
   165.0      0.0    -5.04645
   180.0      0.0    -4.22230
  -180.0     15.0    -4.27526
  -165.0     15.0    -4.46238
  -150.0     15.0    -4.67044
  -135.0     15.0    -4.95535
  -120.0     15.0    -5.30069
  -105.0     15.0    -5.59602
   -90.0     15.0    -5.75823
   -75.0     15.0    -5.72479
   -60.0     15.0    -4.85236
   -45.0     15.0    -3.22837
   -30.0     15.0    -1.33903
   -15.0     15.0     0.05777
     0.0     15.0     0.73245
    15.0     15.0     0.42655
    30.0     15.0    -0.88283
    45.0     15.0    -2.79424
    60.0     15.0    -4.68897
    75.0     15.0    -5.99165
    90.0     15.0    -6.40257
   105.0     15.0    -6.18609
   120.0     15.0    -5.75161
   135.0     15.0    -5.26619
   150.0     15.0    -4.68528
   165.0     15.0    -3.60147
   180.0     15.0    -4.27526
  -180.0     30.0    -2.58364
  -165.0     30.0    -2.68532
  -150.0     30.0    -2.88548
  -135.0     30.0    -3.21036
  -120.0     30.0    -3.61554
  -105.0     30.0    -4.04096
   -90.0     30.0    -4.67162
   -75.0     30.0    -4.69688
   -60.0     30.0    -3.84373
   -45.0     30.0    -2.00765
   -30.0     30.0    -0.06654
   -15.0     30.0     1.55486
     0.0     30.0     2.31188
    15.0     30.0     1.95817
    30.0     30.0     0.38855
    45.0     30.0    -2.03450
    60.0     30.0    -4.04682
    75.0     30.0    -5.15227
    90.0     30.0    -5.26121
   105.0     30.0    -4.81605
   120.0     30.0    -4.15184
   135.0     30.0    -3.32353
   150.0     30.0    -2.03519
   165.0     30.0    -2.46185
   180.0     30.0    -2.58364
  -180.0     45.0    -0.27385
  -165.0     45.0    -0.39460
  -150.0     45.0    -0.69448
  -135.0     45.0    -1.07442
  -120.0     45.0    -1.56505
  -105.0     45.0    -2.42000
   -90.0     45.0    -3.23377
   -75.0     45.0    -3.16790
   -60.0     45.0    -2.11507
   -45.0     45.0    -0.23679
   -30.0     45.0     1.82415
   -15.0     45.0     3.60759
     0.0     45.0     4.30295
    15.0     45.0     3.56251
    30.0     45.0     1.17974
    45.0     45.0    -1.15606
    60.0     45.0    -2.82064
    75.0     45.0    -3.48890
    90.0     45.0    -3.31258
   105.0     45.0    -2.58412
   120.0     45.0    -1.56700
   135.0     45.0    -0.41066
   150.0     45.0    -0.36210
   165.0     45.0    -0.26642
   180.0     45.0    -0.27385
  -180.0     60.0     1.75943
  -165.0     60.0     1.51071
  -150.0     60.0     1.18267
  -135.0     60.0     0.71371
  -120.0     60.0    -0.01095
  -105.0     60.0    -1.08681
   -90.0     60.0    -1.77710
   -75.0     60.0    -1.64228
   -60.0     60.0    -0.72835
   -45.0     60.0     1.34478
   -30.0     60.0     3.44050
   -15.0     60.0     5.31247
     0.0     60.0     5.75134
    15.0     60.0     3.95862
    30.0     60.0     1.88655
    45.0     60.0     0.12792
    60.0     60.0    -1.04403
    75.0     60.0    -1.49008
    90.0     60.0    -1.15692
   105.0     60.0    -0.17426
   120.0     60.0     0.19601
   135.0     60.0     0.86911
   150.0     60.0     1.47445
   165.0     60.0     1.79366
   180.0     60.0     1.75943
  -180.0     75.0     2.77761
  -165.0     75.0     2.53125
  -150.0     75.0     2.18987
  -135.0     75.0     1.60144
  -120.0     75.0     0.70386
  -105.0     75.0    -0.43657
   -90.0     75.0    -1.06238
   -75.0     75.0    -0.95978
   -60.0     75.0     0.19866
   -45.0     75.0     1.98173
   -30.0     75.0     4.24483
   -15.0     75.0     5.94623
     0.0     75.0     5.44049
    15.0     75.0     4.08935
    30.0     75.0     2.74299
    45.0     75.0     1.59048
    60.0     75.0     0.58100
    75.0     75.0     0.13385
    90.0     75.0    -0.48069
   105.0     75.0    -0.17171
   120.0     75.0     0.59458
   135.0     75.0     1.67401
   150.0     75.0     2.55876
   165.0     75.0     2.93260
   180.0     75.0     2.77761
  -180.0     90.0     2.74331
  -165.0     90.0     2.58198
  -150.0     90.0     2.27845
  -135.0     90.0     1.55961
  -120.0     90.0     0.62829
  -105.0     90.0    -0.40231
   -90.0     90.0    -1.11458
   -75.0     90.0    -1.08320
   -60.0     90.0    -0.28428
   -45.0     90.0     1.48393
   -30.0     90.0     3.46788
   -15.0     90.0     4.70554
     0.0     90.0     4.49327
    15.0     90.0     3.96283
    30.0     90.0     3.38653
    45.0     90.0     2.52482
    60.0     90.0     0.58385
    75.0     90.0    -1.04633
    90.0     90.0    -1.24738
   105.0     90.0    -0.60009
   120.0     90.0     0.44790
   135.0     90.0     1.59315
   150.0     90.0     2.65974
   165.0     90.0     2.74719
   180.0     90.0     2.74331
  -180.0    105.0     1.95497
  -165.0    105.0     1.97933
  -150.0    105.0     1.65296
  -135.0    105.0     1.13219
  -120.0    105.0     0.25289
  -105.0    105.0    -0.79043
   -90.0    105.0    -1.75760
   -75.0    105.0    -2.00163
   -60.0    105.0    -1.39335
   -45.0    105.0     0.01966
   -30.0    105.0     1.68837
   -15.0    105.0     2.74754
     0.0    105.0     3.26537
    15.0    105.0     3.45811
    30.0    105.0     3.35587
    45.0    105.0     1.58629
    60.0    105.0    -1.05326
    75.0    105.0    -2.46809
    90.0    105.0    -2.34446
   105.0    105.0    -1.44395
   120.0    105.0    -0.23836
   135.0    105.0     1.09817
   150.0    105.0     1.74844
   165.0    105.0     1.89146
   180.0    105.0     1.95497
  -180.0    120.0     0.94002
  -165.0    120.0     1.14392
  -150.0    120.0     0.94220
  -135.0    120.0     0.42439
  -120.0    120.0    -0.45731
  -105.0    120.0    -1.58770
   -90.0    120.0    -2.45287
   -75.0    120.0    -3.08135
   -60.0    120.0    -2.90997
   -45.0    120.0    -1.83307
   -30.0    120.0    -0.33944
   -15.0    120.0     1.05732
     0.0    120.0     2.15560
    15.0    120.0     2.72433
    30.0    120.0     2.73567
    45.0    120.0    -0.51565
    60.0    120.0    -3.11264
    75.0    120.0    -3.94292
    90.0    120.0    -3.50628
   105.0    120.0    -2.58378
   120.0    120.0    -1.19823
   135.0    120.0    -0.12378
   150.0    120.0     0.39856
   165.0    120.0     0.68619
   180.0    120.0     0.94002
  -180.0    135.0     0.11457
  -165.0    135.0     0.39970
  -150.0    135.0     0.30326
  -135.0    135.0    -0.26169
  -120.0    135.0    -1.38411
  -105.0    135.0    -2.57323
   -90.0    135.0    -3.48022
   -75.0    135.0    -4.14746
   -60.0    135.0    -4.15255
   -45.0    135.0    -3.31494
   -30.0    135.0    -1.78936
   -15.0    135.0    -0.05629
     0.0    135.0     1.29409
    15.0    135.0     1.77725
    30.0    135.0     1.07212
    45.0    135.0    -2.94717
    60.0    135.0    -4.95241
    75.0    135.0    -5.35985
    90.0    135.0    -4.73673
   105.0    135.0    -3.62454
   120.0    135.0    -2.41382
   135.0    135.0    -1.46090
   150.0    135.0    -0.86739
   165.0    135.0    -0.35767
   180.0    135.0     0.11457
  -180.0    150.0    -0.34547
  -165.0    150.0    -0.09173
  -150.0    150.0    -0.34336
  -135.0    150.0    -1.09788
  -120.0    150.0    -2.21872
  -105.0    150.0    -3.43923
   -90.0    150.0    -4.47202
   -75.0    150.0    -5.11368
   -60.0    150.0    -5.15707
   -45.0    150.0    -4.32002
   -30.0    150.0    -2.63712
   -15.0    150.0    -0.72643
     0.0    150.0     0.43227
    15.0    150.0     1.19521
    30.0    150.0    -1.82424
    45.0    150.0    -4.62413
    60.0    150.0    -6.11149
    75.0    150.0    -6.28648
    90.0    150.0    -5.59097
   105.0    150.0    -4.51581
   120.0    150.0    -3.36364
   135.0    150.0    -2.39523
   150.0    150.0    -1.61138
   165.0    150.0    -0.91119
   180.0    150.0    -0.34547
  -180.0    165.0    -0.41413
  -165.0    165.0    -0.38834
  -150.0    165.0    -0.85965
  -135.0    165.0    -1.80556
  -120.0    165.0    -3.01099
  -105.0    165.0    -4.24252
   -90.0    165.0    -5.26258
   -75.0    165.0    -5.87160
   -60.0    165.0    -5.89147
   -45.0    165.0    -4.96362
   -30.0    165.0    -3.08291
   -15.0    165.0    -1.02648
     0.0    165.0     0.32603
    15.0    165.0    -0.08376
    30.0    165.0    -2.94954
    45.0    165.0    -5.22077
    60.0    165.0    -6.49751
    75.0    165.0    -6.55928
    90.0    165.0    -5.89392
   105.0    165.0    -4.87737
   120.0    165.0    -3.73864
   135.0    165.0    -2.64493
   150.0    165.0    -1.68618
   165.0    165.0    -0.90052
   180.0    165.0    -0.41413
  -180.0    180.0    -0.41030
  -165.0    180.0    -0.67838
  -150.0    180.0    -1.32284
  -135.0    180.0    -2.35817
  -120.0    180.0    -3.54199
  -105.0    180.0    -4.74648
   -90.0    180.0    -5.73614
   -75.0    180.0    -6.35918
   -60.0    180.0    -6.35328
   -45.0    180.0    -5.25542
   -30.0    180.0    -3.16646
   -15.0    180.0    -1.09098
     0.0    180.0     1.38787
    15.0    180.0    -1.09098
    30.0    180.0    -3.16646
    45.0    180.0    -5.25542
    60.0    180.0    -6.35328
    75.0    180.0    -6.36063
    90.0    180.0    -5.73614
   105.0    180.0    -4.74648
   120.0    180.0    -3.54199
   135.0    180.0    -2.35817
   150.0    180.0    -1.32284
   165.0    180.0    -0.67848
   180.0    180.0    -0.41030


      ###################################
      ##                               ##
      ##  Atomic Multipole Parameters  ##
      ##                               ##
      ###################################


multipole     1    2    4              -0.22483
                                        0.12482    0.00000    0.32382
                                        0.10689
                                        0.00000   -0.84021
                                       -0.14175    0.00000    0.73332
multipole     2    1    3              -0.09722
                                        0.29356    0.00000    0.16702
                                        0.08444
                                        0.00000   -0.32288
                                       -0.24845    0.00000    0.23844
multipole     2    1  233              -0.31113
                                        0.29356    0.00000    0.16702
                                        0.08444
                                        0.00000   -0.32288
                                       -0.24845    0.00000    0.23844
multipole     2    1  235              -0.18217
                                        0.30273    0.00000    0.14966
                                       -0.00540
                                        0.00000   -0.58515
                                       -0.15734    0.00000    0.59055
multipole     2  231    3              -0.06142
                                        0.15310    0.00000    0.58485
                                       -1.10051
                                        0.00000   -0.52551
                                       -0.23601    0.00000    1.62602
multipole     3    5    2               0.79586
                                       -0.13770    0.00000    0.29742
                                        0.14120
                                        0.00000   -0.29053
                                        0.00632    0.00000    0.14933
multipole     4    1    2               0.13964
                                        0.00855    0.00000   -0.12105
                                        0.00073
                                        0.00000   -0.01705
                                       -0.03619    0.00000    0.01632
multipole     5    3    2              -0.76219
                                        0.02701    0.00000   -0.11735
                                       -0.52988
                                        0.00000    0.17528
                                        0.04451    0.00000    0.35460
multipole     6    2    6               0.07437
                                        0.01314    0.00000   -0.06475
                                        0.08205
                                        0.00000    0.01765
                                        0.00706    0.00000   -0.09970
multipole     7    8   10              -0.15418
                                        0.09038    0.00000    0.44023
                                        0.11417
                                        0.00000   -1.09602
                                       -0.18333    0.00000    0.98185
multipole     7   48   10              -0.14115
                                        0.12232    0.00000    0.52872
                                       -0.43015
                                        0.00000   -0.90054
                                       -0.20486    0.00000    1.33069
multipole     8    7    9   12         -0.17302
                                        0.20245    0.00000    0.14886
                                       -0.19394
                                        0.00000   -0.38164
                                       -0.02856    0.00000    0.57558
multipole     8    7  233   12         -0.36199
                                       -0.12137    0.00000    0.18627
                                       -0.91798
                                        0.00000    0.05376
                                        0.15506    0.00000    0.86422
multipole     8    7  235   12         -0.11441
                                        0.20245    0.00000    0.14886
                                       -0.19394
                                        0.00000   -0.38164
                                       -0.02856    0.00000    0.57558
multipole     8  231    9   12          0.04440
                                        0.09865    0.00000    0.46297
                                       -1.02943
                                        0.00000   -0.69704
                                       -0.49818    0.00000    1.72647
multipole     9   11    8               0.81010
                                       -0.03747    0.00000    0.23377
                                        0.31967
                                        0.00000   -0.32818
                                       -0.03412    0.00000    0.00851
multipole     9   11   48               0.85003
                                       -0.02491    0.00000    0.25777
                                        0.45171
                                        0.00000   -0.12926
                                       -0.10936    0.00000   -0.32245
multipole    10    7    8               0.12044
                                       -0.00442    0.00000   -0.13836
                                        0.03666
                                        0.00000    0.02185
                                       -0.10048    0.00000   -0.05851
multipole    10    7   48               0.12637
                                       -0.02506    0.00000   -0.11601
                                        0.15583
                                        0.00000    0.22096
                                        0.19740    0.00000   -0.37679
multipole    11    9    8              -0.75149
                                       -0.00495    0.00000   -0.15871
                                       -0.53066
                                        0.00000    0.27820
                                        0.10276    0.00000    0.25246
multipole    11    9   48              -0.77014
                                        0.03419    0.00000   -0.15399
                                       -0.55689
                                        0.00000   -0.04221
                                        0.10937    0.00000    0.59910
multipole    12    8   13               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   15               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   19               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   25               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   33               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   37               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   43               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   61               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   70               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   80               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8   89               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  106               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  117               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  127               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  137               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  141               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  147               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  153               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  159               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  167               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  175               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  182               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  192               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12    8  202               0.08514
                                        0.02910    0.00000    0.05267
                                        0.04149
                                        0.00000   -0.00658
                                       -0.00600    0.00000   -0.03491
multipole    12   48   43               0.07834
                                        0.07720    0.00000    0.12670
                                        0.05776
                                        0.00000   -0.05378
                                       -0.00817    0.00000   -0.00398
multipole    13    8   12              -0.15959
                                       -0.01994    0.00000    0.37655
                                       -0.25740
                                        0.00000   -0.29981
                                       -0.04929    0.00000    0.55721
multipole    14   13    8               0.07420
                                       -0.02399    0.00000   -0.11050
                                        0.00932
                                        0.00000    0.05597
                                       -0.06329    0.00000   -0.06529
multipole    15   16    8              -0.06548
                                        0.18863    0.00000   -0.25977
                                        0.38994
                                        0.00000    0.04868
                                       -0.23538    0.00000   -0.43862
multipole    16   15    8               0.07117
                                       -0.02192    0.00000   -0.09528
                                       -0.12013
                                        0.00000    0.08153
                                       -0.00394    0.00000    0.03860
multipole    17   15    8              -0.16832
                                       -0.08635    0.00000    0.33587
                                       -0.27420
                                        0.00000   -0.17464
                                       -0.25682    0.00000    0.44884
multipole    18   17   15               0.06566
                                       -0.02088    0.00000   -0.10166
                                        0.01623
                                        0.00000    0.04728
                                        0.00231    0.00000   -0.06351
multipole    19    8   21              -0.13173
                                        0.18350    0.00000    0.29386
                                        0.00999
                                        0.00000   -0.65954
                                       -0.21063    0.00000    0.64955
multipole    20   19   21               0.06916
                                        0.05664    0.00000   -0.13463
                                        0.04076
                                        0.00000    0.11032
                                        0.10592    0.00000   -0.15108
multipole    21   22   19              -0.02887
                                        0.11379    0.00000   -0.33524
                                        0.32588
                                        0.00000    0.22134
                                       -0.11449    0.00000   -0.54722
multipole    22   21   19               0.08767
                                       -0.01393    0.00000   -0.07934
                                       -0.05712
                                        0.00000    0.16068
                                        0.00191    0.00000   -0.10356
multipole    23   21   19              -0.17741
                                       -0.04820    0.00000    0.27453
                                       -0.28572
                                        0.00000   -0.21363
                                       -0.06394    0.00000    0.49935
multipole    24   23   21               0.05874
                                        0.00734    0.00000   -0.09498
                                        0.04549
                                        0.00000    0.01789
                                        0.01922    0.00000   -0.06338
multipole    25   26   27              -0.08555
                                       -0.10128    0.00000   -0.23524
                                        0.20588
                                        0.00000    0.35986
                                       -0.06061    0.00000   -0.56574
multipole    26   25   27               0.08230
                                        0.04789    0.00000   -0.08924
                                        0.01057
                                        0.00000    0.11342
                                        0.06251    0.00000   -0.12399
multipole    27   31   25              -0.11939
                                        0.32104    0.00000    0.10406
                                        0.26796
                                        0.00000   -0.55943
                                       -0.32286    0.00000    0.29147
multipole    28   27   25               0.07096
                                       -0.02472    0.00000   -0.08760
                                       -0.02354
                                        0.00000    0.07793
                                       -0.09866    0.00000   -0.05439
multipole    29   25   27              -0.17217
                                        0.01336    0.00000    0.38349
                                       -0.31978
                                        0.00000   -0.13786
                                        0.11205    0.00000    0.45764
multipole    30   29   25               0.06823
                                       -0.01947    0.00000   -0.09588
                                        0.04494
                                        0.00000    0.05530
                                        0.00419    0.00000   -0.10024
multipole    31   27   32              -0.16882
                                        0.02110    0.00000    0.24422
                                       -0.25302
                                        0.00000   -0.13559
                                        0.02868    0.00000    0.38861
multipole    32   31   27               0.06001
                                        0.00612    0.00000   -0.10697
                                        0.03438
                                        0.00000    0.06956
                                       -0.03511    0.00000   -0.10394
multipole    33    8   35               0.16229
                                        0.22880    0.00000    0.30596
                                        0.12912
                                        0.00000   -0.06006
                                       -0.66812    0.00000   -0.06906
multipole    34   33    8               0.05600
                                       -0.14483    0.00000   -0.07492
                                       -0.00102
                                        0.00000   -0.06173
                                       -0.03044    0.00000    0.06275
multipole    35   36   33              -0.44145
                                        0.23477    0.00000    0.13825
                                        0.06803
                                        0.00000   -0.49343
                                       -0.14987    0.00000    0.42541
multipole    36   35   33               0.23017
                                       -0.12968    0.00000    0.01379
                                       -0.06824
                                        0.00000   -0.11228
                                       -0.08774    0.00000    0.18052
multipole    37   39    8               0.10834
                                        0.24982    0.00000    0.27937
                                        0.08477
                                        0.00000   -0.58656
                                       -0.19826    0.00000    0.50179
multipole    38   37   39               0.06805
                                        0.02048    0.00000   -0.09220
                                       -0.00303
                                        0.00000   -0.01146
                                       -0.05463    0.00000    0.01449
multipole    39   40   37              -0.40177
                                        0.35143    0.00000    0.04331
                                        0.31002
                                        0.00000   -0.84324
                                       -0.38927    0.00000    0.53322
multipole    40   39   37               0.24568
                                       -0.10665    0.00000    0.01809
                                       -0.15372
                                        0.00000   -0.05391
                                        0.00480    0.00000    0.20763
multipole    41   37   42              -0.16081
                                       -0.01376    0.00000    0.41130
                                       -0.26001
                                        0.00000   -0.18818
                                       -0.12509    0.00000    0.44819
multipole    42   41   37               0.06784
                                       -0.05387    0.00000   -0.13358
                                        0.00765
                                        0.00000    0.10555
                                       -0.15569    0.00000   -0.11320
multipole    43    8   45              -0.30144
                                        0.07672    0.00000    0.42645
                                       -0.35442
                                        0.00000   -0.51033
                                       -0.35331    0.00000    0.86475
multipole    43    8   47              -0.30144
                                        0.07672    0.00000    0.42645
                                       -0.35442
                                        0.00000   -0.51033
                                       -0.35331    0.00000    0.86475
multipole    43   48   49              -0.08027
                                        0.13313    0.00000    0.28094
                                       -0.73267
                                        0.00000    0.44660
                                       -0.05681    0.00000    0.28607
multipole    44   43    8               0.14381
                                       -0.07645    0.00000    0.01699
                                       -0.09149
                                        0.00000   -0.03033
                                       -0.15284    0.00000    0.12182
multipole    44   43   48               0.00975
                                       -0.14579    0.00000   -0.00614
                                       -0.22868
                                        0.00000    0.13700
                                       -0.21056    0.00000    0.09168
multipole    45   46   43              -0.03264
                                        0.49460    0.00000   -0.04347
                                        1.63727
                                        0.00000   -2.40778
                                       -0.45566    0.00000    0.77051
multipole    46   45   43               0.10947
                                       -0.10504    0.00000   -0.03823
                                       -0.08656
                                        0.00000   -0.01943
                                       -0.23485    0.00000    0.10599
multipole    47   47   43               0.07683
                                        0.30088    0.00000    0.27072
                                        1.04951
                                        0.00000   -2.57126
                                       -0.50502    0.00000    1.52175
multipole    48    7    9   12         -0.23113
                                        0.27663    0.00000    0.12440
                                       -0.20625
                                        0.00000   -0.28389
                                        0.15574    0.00000    0.49014
multipole    48    7  233   12          0.60118
                                        0.27663    0.00000    0.12440
                                       -0.20625
                                        0.00000   -0.28389
                                        0.15574    0.00000    0.49014
multipole    49   43   48              -0.85155
                                        0.15871    0.00000    0.56464
                                       -1.49210
                                        0.00000   -1.58110
                                       -0.20334    0.00000    3.07320
multipole    50   51   59              -0.13558
                                        0.39735    0.00000    0.41339
                                        0.27347
                                        0.00000   -1.36200
                                       -0.22439    0.00000    1.08853
multipole    51   50   52   54         -0.26070
                                        0.08712    0.00000    0.33920
                                       -0.37488
                                        0.00000   -0.33230
                                        0.21154    0.00000    0.70718
multipole    51   50  233   54         -0.37342
                                        0.08712    0.00000    0.33920
                                       -0.37488
                                        0.00000   -0.33230
                                        0.21154    0.00000    0.70718
multipole    52   53   51               0.90470
                                        0.06004    0.00000    0.25405
                                        0.35305
                                        0.00000   -0.30405
                                        0.20110    0.00000   -0.04900
multipole    53   52   51              -0.76984
                                       -0.05467    0.00000   -0.23073
                                       -0.46593
                                        0.00000    0.27107
                                       -0.00005    0.00000    0.19486
multipole    54   51   55               0.09550
                                        0.05530    0.00000    0.17764
                                       -0.07670
                                        0.00000   -0.13879
                                       -0.00767    0.00000    0.21549
multipole    55   51   57              -0.13205
                                        0.32755    0.00000    0.46077
                                        0.17752
                                        0.00000   -0.61073
                                       -0.35102    0.00000    0.43321
multipole    55  241   57              -0.12680
                                        0.31318    0.00000    0.41341
                                        0.26389
                                        0.00000   -0.70262
                                       -0.10425    0.00000    0.43873
multipole    56   55   51               0.09686
                                       -0.08812    0.00000   -0.10597
                                       -0.00614
                                        0.00000    0.10040
                                       -0.16651    0.00000   -0.09426
multipole    56   55  241               0.11510
                                        0.00090    0.00000   -0.07548
                                       -0.00989
                                        0.00000   -0.01502
                                        0.00444    0.00000    0.02491
multipole    57   55   59              -0.16591
                                        0.43907    0.00000    0.23052
                                        0.21670
                                        0.00000   -0.57363
                                       -0.09937    0.00000    0.35693
multipole    57   55  245              -0.14420
                                        0.42355    0.00000    0.25235
                                        0.34755
                                        0.00000   -0.69761
                                       -0.08434    0.00000    0.35006
multipole    58   57   55               0.08886
                                       -0.00234    0.00000   -0.08982
                                        0.07325
                                        0.00000    0.05889
                                       -0.03278    0.00000   -0.13214
multipole    58   57  245               0.12121
                                       -0.00270    0.00000   -0.06030
                                       -0.00918
                                        0.00000   -0.01098
                                        0.00194    0.00000    0.02016
multipole    59   50   57              -0.03932
                                        0.17998    0.00000    0.43068
                                       -0.10117
                                        0.00000   -0.38138
                                        0.01443    0.00000    0.48255
multipole    60   59   50               0.06588
                                       -0.03733    0.00000   -0.07724
                                        0.02650
                                        0.00000   -0.00717
                                       -0.08958    0.00000   -0.01933
multipole    61   63    8              -0.11216
                                        0.38620    0.00000    0.08762
                                        0.47377
                                        0.00000   -0.59723
                                       -0.25956    0.00000    0.12346
multipole    62   61   63               0.09205
                                        0.04929    0.00000   -0.11825
                                        0.02631
                                        0.00000    0.10410
                                        0.09871    0.00000   -0.13041
multipole    63   61   64              -0.04246
                                        0.01409    0.00000    0.21398
                                        0.05117
                                        0.00000   -0.26849
                                       -0.07370    0.00000    0.21732
multipole    64   63   66               0.04681
                                        0.17684    0.00000    0.12228
                                        0.48702
                                        0.00000   -0.40448
                                        0.34751    0.00000   -0.08254
multipole    65   64   63               0.00023
                                       -0.01320    0.00000   -0.17110
                                        0.08705
                                        0.00000    0.10296
                                        0.00895    0.00000   -0.19001
multipole    66   64   68              -0.02656
                                       -0.02571    0.00000    0.00360
                                       -0.00110
                                        0.00000   -0.14275
                                       -0.09779    0.00000    0.14385
multipole    67   66   64               0.00719
                                        0.00137    0.00000   -0.18600
                                        0.05207
                                        0.00000   -0.00064
                                        0.00676    0.00000   -0.05143
multipole    68   66   66              -0.02866
                                       -0.02759    0.00000   -0.00699
                                       -0.01545
                                        0.00000   -0.19273
                                       -0.12394    0.00000    0.20818
multipole    69   68   66               0.00685
                                       -0.00666    0.00000   -0.18375
                                        0.05602
                                        0.00000   -0.00631
                                        0.00545    0.00000   -0.04971
multipole    70   72    8              -0.10731
                                        0.38034    0.00000    0.08372
                                        0.50724
                                        0.00000   -0.62647
                                       -0.18907    0.00000    0.11923
multipole    71   70   72               0.09011
                                        0.04522    0.00000   -0.12453
                                        0.04012
                                        0.00000    0.06958
                                        0.07019    0.00000   -0.10970
multipole    72   70   73              -0.05309
                                        0.06610    0.00000    0.18549
                                        0.04921
                                        0.00000   -0.33224
                                       -0.11333    0.00000    0.28303
multipole    73   72   75               0.07134
                                        0.24299    0.00000    0.14377
                                        0.57252
                                        0.00000   -0.47541
                                        0.42651    0.00000   -0.09711
multipole    74   73   72               0.00345
                                       -0.01772    0.00000   -0.17185
                                        0.14345
                                        0.00000    0.12181
                                        0.02031    0.00000   -0.26526
multipole    75   73   77              -0.12512
                                       -0.11733    0.00000    0.00000
                                        0.04101
                                        0.00000   -0.14140
                                       -0.08647    0.00000    0.10039
multipole    76   75   73               0.01090
                                        0.02337    0.00000   -0.19140
                                        0.03985
                                        0.00000    0.06312
                                       -0.01694    0.00000   -0.10297
multipole    77   78   75               0.35744
                                        0.02001    0.00000    0.42701
                                       -0.19251
                                        0.00000   -0.38151
                                        0.27981    0.00000    0.57402
multipole    78   77   79              -0.46438
                                        0.19713    0.00000    0.00614
                                        0.42591
                                        0.00000   -0.39511
                                       -0.34625    0.00000   -0.03080
multipole    79   78   77               0.22899
                                       -0.06109    0.00000   -0.01309
                                       -0.11816
                                        0.00000   -0.15815
                                        0.04555    0.00000    0.27631
multipole    80    8   82              -0.14909
                                        0.11344    0.00000    0.39909
                                        0.03628
                                        0.00000   -0.47377
                                        0.02261    0.00000    0.43749
multipole    81   80   81               0.07954
                                        0.01404    0.00000   -0.09182
                                        0.04405
                                        0.00000    0.01752
                                        0.03456    0.00000   -0.06157
multipole    82   80   83              -0.21590
                                        0.00000    0.00000    0.19889
                                        0.43753
                                        0.00000   -0.84736
                                        0.00000    0.00000    0.40983
multipole    83   82   85              -0.04549
                                        0.32216    0.00000    0.13839
                                        0.73596
                                        0.00000   -0.60409
                                        0.25285    0.00000   -0.13187
multipole    84   83   82               0.01759
                                       -0.02121    0.00000   -0.10542
                                        0.08078
                                        0.00000    0.10479
                                       -0.00935    0.00000   -0.18557
multipole    85   83   87              -0.15892
                                        0.03729    0.00000    0.15665
                                        0.59432
                                        0.00000   -0.78857
                                        0.29227    0.00000    0.19425
multipole    86   85   87              -0.03324
                                        0.03781    0.00000   -0.15021
                                        0.16274
                                        0.00000    0.11972
                                        0.03565    0.00000   -0.28246
multipole    87   88   85               0.62054
                                        0.00000    0.00000    0.43281
                                        0.41363
                                        0.00000   -0.34255
                                        0.00000    0.00000   -0.07108
multipole    88   87   85              -0.91150
                                        0.00000    0.00000   -0.26161
                                       -0.27372
                                        0.00000    0.34305
                                        0.00000    0.00000   -0.06933
multipole    89    8   91              -0.16339
                                       -0.12477    0.00000    0.53121
                                       -0.30582
                                        0.00000   -0.43704
                                        0.40713    0.00000    0.74286
multipole    90   89   91               0.08819
                                        0.06555    0.00000   -0.10115
                                       -0.10247
                                        0.00000    0.15243
                                        0.03296    0.00000   -0.04996
multipole    91   92   94              -0.14846
                                       -0.32645    0.00000   -0.13851
                                        0.16205
                                        0.00000   -0.20283
                                        0.53464    0.00000    0.04078
multipole    92   95   91               0.09676
                                        0.02936    0.00000    0.48291
                                       -0.05924
                                        0.00000   -0.39203
                                        0.57677    0.00000    0.45127
multipole    93   92   95               0.02989
                                       -0.00149    0.00000   -0.20147
                                       -0.02175
                                        0.00000    0.07153
                                        0.04234    0.00000   -0.04978
multipole    94   91   98               0.00548
                                        0.26321    0.00000   -0.09442
                                        0.27216
                                        0.00000   -0.00293
                                       -0.14493    0.00000   -0.26923
multipole    95   92   97              -0.03588
                                        0.13305    0.00000    0.32467
                                        0.47653
                                        0.00000   -0.93115
                                        0.99076    0.00000    0.45462
multipole    96   95   92               0.08805
                                        0.00174    0.00000   -0.22101
                                        0.07562
                                        0.00000    0.06395
                                        0.01790    0.00000   -0.13957
multipole    97   95  100               0.27382
                                       -0.01350    0.00000   -0.01650
                                       -0.51808
                                        0.00000   -0.74455
                                       -0.04923    0.00000    1.26263
multipole    98  102   94              -0.20979
                                       -0.34622    0.00000    0.10170
                                       -0.36832
                                        0.00000    0.33475
                                        0.35509    0.00000    0.03357
multipole    99   98  102               0.00102
                                        0.00505    0.00000   -0.17714
                                        0.06900
                                        0.00000    0.06904
                                       -0.03904    0.00000   -0.13804
multipole   100  104   97              -0.11296
                                       -0.13812    0.00000   -0.03163
                                       -0.06766
                                        0.00000    0.01833
                                       -0.00493    0.00000    0.04933
multipole   101  100  104               0.00446
                                        0.00790    0.00000   -0.18506
                                        0.03323
                                        0.00000    0.01707
                                       -0.02718    0.00000   -0.05030
multipole   102   98  104               0.06902
                                        0.27757    0.00000    0.10967
                                        0.56142
                                        0.00000   -0.53401
                                        0.44766    0.00000   -0.02741
multipole   103  102   98               0.00474
                                        0.00358    0.00000   -0.18227
                                        0.04939
                                        0.00000    0.01436
                                        0.01490    0.00000   -0.06375
multipole   104  102  100              -0.01987
                                        0.02835    0.00000    0.04818
                                        0.20675
                                        0.00000   -0.27604
                                        0.00990    0.00000    0.06929
multipole   105  104  102               0.00374
                                        0.00125    0.00000   -0.18390
                                        0.04446
                                        0.00000   -0.01545
                                       -0.01437    0.00000   -0.02901
multipole   106    8  108              -0.10044
                                        0.08050    0.00000    0.31541
                                       -0.04428
                                        0.00000   -0.57160
                                       -0.15478    0.00000    0.61588
multipole   107  106  108               0.11241
                                        0.00094    0.00000   -0.14233
                                        0.02742
                                        0.00000    0.10370
                                        0.00990    0.00000   -0.13112
multipole   108  106  109               0.16205
                                        0.16361    0.00000    0.27606
                                        0.18556
                                        0.00000   -0.33531
                                       -0.34553    0.00000    0.14975
multipole   109  108  113              -0.10600
                                        0.06737    0.00000    0.13621
                                        0.52452
                                        0.00000   -0.63296
                                        0.39427    0.00000    0.10844
multipole   110  109  108               0.20247
                                       -0.04159    0.00000   -0.13094
                                        0.01709
                                        0.00000   -0.01033
                                       -0.01012    0.00000   -0.00676
multipole   111  108  115               0.04476
                                        0.19634    0.00000   -0.22474
                                        0.58980
                                        0.00000   -0.17266
                                       -0.21786    0.00000   -0.41714
multipole   112  111  115               0.08897
                                        0.00602    0.00000   -0.21192
                                        0.01372
                                        0.00000    0.02533
                                        0.02467    0.00000   -0.03905
multipole   113  109  115               0.39335
                                        0.28599    0.00000    0.15648
                                        0.09076
                                        0.00000   -0.28047
                                       -0.26152    0.00000    0.18971
multipole   114  113  115               0.09278
                                       -0.00527    0.00000   -0.23564
                                        0.04222
                                        0.00000   -0.00568
                                       -0.00439    0.00000   -0.03654
multipole   115  111  113              -0.08186
                                       -0.06539    0.00000    0.06061
                                        0.34422
                                        0.00000   -0.66901
                                        0.31146    0.00000    0.32479
multipole   116  115  113               0.14211
                                       -0.01142    0.00000   -0.24975
                                        0.05866
                                        0.00000   -0.00893
                                       -0.00227    0.00000   -0.04973
multipole   117    8  119              -0.10846
                                        0.06915    0.00000    0.34129
                                       -0.13274
                                        0.00000   -0.55071
                                       -0.03813    0.00000    0.68345
multipole   118  117  119               0.08849
                                        0.06029    0.00000   -0.13194
                                        0.02277
                                        0.00000    0.08141
                                        0.08284    0.00000   -0.10418
multipole   119  117  120               0.04660
                                        0.25752    0.00000    0.19196
                                        0.36112
                                        0.00000   -0.48799
                                       -0.42867    0.00000    0.12687
multipole   120  119  124              -0.16996
                                        0.07326    0.00000    0.05464
                                        0.44029
                                        0.00000   -0.71866
                                        0.32632    0.00000    0.27837
multipole   121  120  124               0.13110
                                        0.02575    0.00000   -0.14558
                                        0.00657
                                        0.00000   -0.01572
                                        0.00411    0.00000    0.00915
multipole   122  119  126               0.08866
                                        0.18994    0.00000   -0.17003
                                        0.37335
                                        0.00000   -0.27499
                                       -0.18897    0.00000   -0.09836
multipole   123  122  126               0.02049
                                        0.03210    0.00000   -0.20183
                                        0.00721
                                        0.00000    0.04324
                                        0.00501    0.00000   -0.05045
multipole   124  120  126               0.34991
                                        0.30562    0.00000    0.16281
                                       -0.05686
                                        0.00000   -0.39010
                                       -0.16611    0.00000    0.44696
multipole   125  124  126               0.03404
                                        0.02495    0.00000   -0.20391
                                        0.00173
                                        0.00000   -0.02858
                                        0.01013    0.00000    0.02685
multipole   126  122  124              -0.50635
                                        0.18134    0.00000    0.26499
                                       -0.07393
                                        0.00000   -0.04385
                                       -0.12432    0.00000    0.11778
multipole   127    8  129              -0.11338
                                        0.22119    0.00000    0.30702
                                       -0.06116
                                        0.00000   -0.60019
                                        0.03111    0.00000    0.66135
multipole   128  127  129               0.08990
                                        0.03505    0.00000   -0.12631
                                        0.02431
                                        0.00000    0.08784
                                        0.08710    0.00000   -0.11215
multipole   129  127  130               0.20611
                                        0.29762    0.00000    0.05580
                                        0.13649
                                        0.00000   -0.54039
                                       -0.01448    0.00000    0.40390
multipole   130  129  133              -0.56975
                                        0.21589    0.00000    0.18560
                                        0.00927
                                        0.00000   -0.03572
                                       -0.10374    0.00000    0.02645
multipole   131  129  135              -0.02289
                                        0.28325    0.00000   -0.15990
                                        0.65227
                                        0.00000   -0.37726
                                       -0.23589    0.00000   -0.27501
multipole   132  131  135               0.03876
                                       -0.00320    0.00000   -0.17520
                                        0.06180
                                        0.00000    0.00375
                                       -0.07542    0.00000   -0.06555
multipole   133  130  135               0.33110
                                        0.22778    0.00000    0.23060
                                        0.19791
                                        0.00000   -0.37320
                                       -0.38619    0.00000    0.17529
multipole   134  133  135               0.03651
                                       -0.02898    0.00000   -0.20060
                                        0.03456
                                        0.00000   -0.00906
                                       -0.01436    0.00000   -0.02550
multipole   135  131  133              -0.13121
                                       -0.01623    0.00000    0.07052
                                        0.22148
                                        0.00000   -0.71860
                                        0.21693    0.00000    0.49712
multipole   136  135  133               0.10796
                                        0.00816    0.00000   -0.20977
                                        0.04667
                                        0.00000   -0.02930
                                        0.04575    0.00000   -0.01737
multipole   137    8  139              -0.33682
                                       -0.09973    0.00000    0.38664
                                       -0.56710
                                        0.00000   -0.27996
                                        0.12769    0.00000    0.84706
multipole   138  137  139               0.04965
                                        0.04752    0.00000   -0.13859
                                        0.11994
                                        0.00000    0.04796
                                        0.03771    0.00000   -0.16790
multipole   139  140 -140               1.01811
                                       -0.00488    0.00000   -0.15412
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   140  139  137              -0.85879
                                       -0.08949    0.00000   -0.07764
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   141    8  143              -0.17303
                                        0.44945    0.00000    0.30728
                                       -0.07251
                                        0.00000   -0.67903
                                        0.08175    0.00000    0.75154
multipole   142  141  143               0.10446
                                       -0.08174    0.00000   -0.18518
                                       -0.00080
                                        0.00000    0.13227
                                       -0.04979    0.00000   -0.13147
multipole   143  144  145               0.85923
                                        0.11304    0.00000    0.10040
                                        0.31104
                                        0.00000   -0.30726
                                       -0.25489    0.00000   -0.00378
multipole   144  143  145              -0.63514
                                       -0.04784    0.00000   -0.07224
                                       -0.41873
                                        0.00000    0.25275
                                       -0.15048    0.00000    0.16598
multipole   145  143  146              -0.43939
                                        0.21541    0.00000   -0.02221
                                        0.36981
                                        0.00000   -0.36139
                                       -0.26932    0.00000   -0.00842
multipole   146  145  143               0.24242
                                       -0.05297    0.00000   -0.03262
                                       -0.15414
                                        0.00000   -0.15395
                                       -0.00414    0.00000    0.30809
multipole   147    8  149              -0.15914
                                        0.20696    0.00000    0.39766
                                        0.04998
                                        0.00000   -0.56038
                                        0.16613    0.00000    0.51040
multipole   148  147  149               0.10149
                                        0.02977    0.00000   -0.15557
                                       -0.01454
                                        0.00000    0.17295
                                        0.10589    0.00000   -0.15841
multipole   149  151  150               0.76984
                                        0.28242    0.00000   -0.15875
                                        0.01831
                                        0.00000   -0.16686
                                        0.09056    0.00000    0.14855
multipole   150  149  151              -0.74313
                                       -0.04904    0.00000   -0.17638
                                       -0.37631
                                        0.00000    0.37623
                                       -0.05205    0.00000    0.00008
multipole   151  149  150              -0.28938
                                        0.12918    0.00000   -0.08416
                                        0.48804
                                        0.00000   -0.69748
                                       -0.01277    0.00000    0.20944
multipole   152  151  149               0.14092
                                       -0.00409    0.00000   -0.17075
                                        0.00327
                                        0.00000   -0.00106
                                        0.19714    0.00000   -0.00221
multipole   153    8  155              -0.10324
                                        0.09831    0.00000    0.33772
                                        0.29786
                                        0.00000   -0.73875
                                       -0.12488    0.00000    0.44089
multipole   154  153  155               0.06492
                                        0.02097    0.00000   -0.14796
                                        0.09102
                                        0.00000    0.12078
                                       -0.12342    0.00000   -0.21180
multipole   155  153  157              -0.35416
                                       -0.08035    0.00000    0.17352
                                       -0.33142
                                        0.00000   -0.02713
                                        0.00000    0.00000    0.35855
multipole   156  155  157               0.06420
                                        0.07412    0.00000   -0.09544
                                        0.40117
                                        0.00000   -0.10252
                                        0.00000    0.00000   -0.29865
multipole   157  158 -158               1.05691
                                       -0.08698    0.00000    0.00867
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   158  157  158              -0.89737
                                       -0.03977    0.00000   -0.09456
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   159    8  161              -0.10668
                                        0.29642    0.00000    0.28365
                                        0.39855
                                        0.00000   -0.41154
                                        0.04273    0.00000    0.01299
multipole   160  159    8               0.09513
                                       -0.05887    0.00000   -0.04640
                                       -0.03611
                                        0.00000   -0.02016
                                       -0.01670    0.00000    0.05627
multipole   161  159  163              -0.23637
                                        0.34729    0.00000    0.09253
                                       -0.40537
                                        0.00000   -0.00346
                                       -0.52649    0.00000    0.40883
multipole   162  161  163               0.10600
                                       -0.20910    0.00000   -0.08841
                                       -0.14892
                                        0.00000    0.17218
                                       -0.34482    0.00000   -0.02326
multipole   163  164  165               0.90086
                                        0.11874    0.00000    0.10181
                                       -0.08663
                                        0.00000    0.03649
                                       -0.12175    0.00000    0.05014
multipole   164  163  165              -0.60677
                                        0.02332    0.00000   -0.04031
                                       -0.46754
                                        0.00000    0.22363
                                        0.10574    0.00000    0.24391
multipole   165  163  164              -0.53100
                                        0.05936    0.00000   -0.00202
                                        0.03421
                                        0.00000    0.00132
                                        0.18532    0.00000   -0.03553
multipole   166  165  163               0.24071
                                       -0.19323    0.00000   -0.29128
                                        0.07422
                                        0.00000    0.09278
                                       -0.19140    0.00000   -0.16700
multipole   167    8  169              -0.07575
                                        0.20886    0.00000    0.29708
                                        0.20862
                                        0.00000   -0.53751
                                       -0.38031    0.00000    0.32889
multipole   168  167    8               0.09458
                                       -0.04304    0.00000   -0.08810
                                       -0.17118
                                        0.00000    0.05903
                                        0.10731    0.00000    0.11215
multipole   169  167  171              -0.27255
                                        0.04245    0.00000    0.29740
                                        0.02354
                                        0.00000   -0.03872
                                        0.02994    0.00000    0.01518
multipole   170  169  171               0.10149
                                        0.02977    0.00000   -0.15557
                                       -0.01454
                                        0.00000    0.17295
                                        0.10589    0.00000   -0.15841
multipole   171  173  172               0.76984
                                        0.28242    0.00000   -0.15875
                                        0.01831
                                        0.00000   -0.16686
                                        0.09056    0.00000    0.14855
multipole   172  171  173              -0.74313
                                       -0.04904    0.00000   -0.17638
                                       -0.37631
                                        0.00000    0.37623
                                       -0.05205    0.00000    0.00008
multipole   173  171  172              -0.28938
                                        0.12918    0.00000   -0.08416
                                        0.48804
                                        0.00000   -0.69748
                                       -0.01277    0.00000    0.20944
multipole   174  173  171               0.14092
                                       -0.00409    0.00000   -0.17075
                                        0.00327
                                        0.00000   -0.00106
                                        0.19714    0.00000   -0.00221
multipole   175    8  177              -0.14215
                                        0.12943    0.00000    0.31867
                                        0.18924
                                        0.00000   -0.71964
                                       -0.21203    0.00000    0.53040
multipole   176  175  177               0.08383
                                        0.08229    0.00000   -0.07905
                                       -0.01759
                                        0.00000    0.08087
                                        0.02604    0.00000   -0.06328
multipole   177  179  175              -0.22038
                                        0.23726    0.00000   -0.16299
                                        0.36732
                                        0.00000   -0.44121
                                       -0.37426    0.00000    0.07389
multipole   178  177  179               0.07884
                                        0.03596    0.00000   -0.06087
                                       -0.03114
                                        0.00000    0.01997
                                        0.06174    0.00000    0.01117
multipole   179  177  180               0.14780
                                        0.54807    0.00000    0.51535
                                        1.30094
                                        0.00000   -2.99992
                                       -0.76856    0.00000    1.69898
multipole   180  179  181              -0.21776
                                        0.01805    0.00000    0.06009
                                       -0.28800
                                        0.00000    0.12124
                                       -0.18981    0.00000    0.16676
multipole   181  180  179               0.05672
                                        0.01469    0.00000   -0.10715
                                        0.09749
                                        0.00000   -0.02112
                                        0.03039    0.00000   -0.07637
multipole   182    8  184              -0.08975
                                        0.20388    0.00000    0.25551
                                        0.08393
                                        0.00000   -0.69546
                                       -0.19022    0.00000    0.61153
multipole   183  182  184               0.08202
                                        0.03852    0.00000   -0.09266
                                        0.00460
                                        0.00000    0.05897
                                       -0.05937    0.00000   -0.06357
multipole   184  186  182              -0.21565
                                        0.30101    0.00000    0.17456
                                        0.15633
                                        0.00000   -0.62215
                                       -0.22099    0.00000    0.46582
multipole   185  184  185               0.08265
                                        0.03106    0.00000   -0.09591
                                        0.04227
                                        0.00000    0.02301
                                        0.08405    0.00000   -0.06528
multipole   186  188  184              -0.02245
                                        0.22603    0.00000    0.33147
                                        0.15933
                                        0.00000   -0.60475
                                       -0.25158    0.00000    0.44542
multipole   187  186  187               0.07982
                                       -0.02778    0.00000   -0.06593
                                        0.07962
                                        0.00000    0.06169
                                        0.09443    0.00000   -0.14131
multipole   188  190  186              -0.01637
                                        0.17686    0.00000    0.14050
                                       -0.14798
                                        0.00000   -0.81782
                                       -0.11856    0.00000    0.96580
multipole   189  188  189               0.09822
                                       -0.04651    0.00000   -0.07618
                                       -0.01229
                                        0.00000   -0.00776
                                       -0.13903    0.00000    0.02005
multipole   190  188  191               0.10000
                                        0.03808    0.00000    0.22325
                                       -0.22293
                                        0.00000   -0.22476
                                        0.00964    0.00000    0.44769
multipole   191  190  188               0.20727
                                        0.03469    0.00000   -0.12995
                                       -0.04317
                                        0.00000   -0.08341
                                        0.18046    0.00000    0.12658
multipole   192    8  194              -0.14610
                                        0.04867    0.00000    0.34432
                                        0.09079
                                        0.00000   -0.55511
                                       -0.18572    0.00000    0.46432
multipole   193  192  194               0.06260
                                        0.02182    0.00000   -0.14723
                                        0.09446
                                        0.00000    0.10407
                                       -0.12330    0.00000   -0.19853
multipole   194  196  192              -0.13934
                                        0.20906    0.00000    0.01228
                                        0.22665
                                        0.00000   -0.42949
                                       -0.27135    0.00000    0.20284
multipole   195  194  195               0.07212
                                        0.00697    0.00000   -0.07338
                                        0.14505
                                        0.00000   -0.06944
                                        0.16099    0.00000   -0.07561
multipole   196  198  194              -0.12642
                                        0.20234    0.00000    0.13874
                                        0.10501
                                        0.00000   -0.54234
                                       -0.24399    0.00000    0.43733
multipole   197  196  197               0.06261
                                       -0.00129    0.00000   -0.07060
                                        0.02999
                                        0.00000    0.03807
                                        0.03535    0.00000   -0.06806
multipole   198  200  196               0.04761
                                        0.09262    0.00000   -0.03821
                                        0.29739
                                        0.00000   -0.34473
                                       -0.15365    0.00000    0.04734
multipole   199  198  199               0.02822
                                       -0.09076    0.00000   -0.01870
                                        0.06435
                                        0.00000   -0.07401
                                       -0.07089    0.00000    0.00966
multipole   200  198  201              -0.20360
                                        0.12136    0.00000    0.24914
                                       -0.14128
                                        0.00000   -0.79272
                                        0.08516    0.00000    0.93400
multipole   201  200  198               0.08988
                                        0.10885    0.00000   -0.33633
                                        0.25980
                                        0.00000   -0.00123
                                        0.38791    0.00000   -0.25857
multipole   202    8  204              -0.09097
                                        0.31208    0.00000    0.34711
                                        0.26802
                                        0.00000   -0.71303
                                        0.08629    0.00000    0.44501
multipole   203  202  204               0.07770
                                       -0.02034    0.00000   -0.11842
                                       -0.07729
                                        0.00000    0.11869
                                       -0.15759    0.00000   -0.04140
multipole   204  202  206              -0.21722
                                        0.33699    0.00000    0.20660
                                        0.18535
                                        0.00000   -0.29250
                                       -0.24900    0.00000    0.10715
multipole   205  204  205               0.08258
                                        0.04941    0.00000   -0.09810
                                        0.14471
                                        0.00000   -0.09703
                                        0.04718    0.00000   -0.04768
multipole   206  204  208              -0.01121
                                        0.30142    0.00000    0.15967
                                        0.51002
                                        0.00000   -0.64074
                                       -0.39142    0.00000    0.13072
multipole   207  206  207               0.07564
                                       -0.04299    0.00000   -0.06409
                                        0.03411
                                        0.00000    0.16042
                                        0.14615    0.00000   -0.19453
multipole   208  206  210              -0.27936
                                       -0.19230    0.00000    0.47650
                                        0.08293
                                        0.00000   -0.98827
                                        0.10755    0.00000    0.90534
multipole   209  208  210               0.15532
                                       -0.11954    0.00000   -0.12082
                                        0.01544
                                        0.00000   -0.19834
                                       -0.04507    0.00000    0.18290
multipole   210  211 -211               0.98761
                                       -0.04689    0.00000   -0.18100
                                       -0.19176
                                        0.00000   -0.18380
                                       -0.09197    0.00000    0.37556
multipole   211  210  211              -0.29224
                                       -0.02579    0.00000   -0.18729
                                        0.46099
                                        0.00000   -0.80088
                                       -0.06689    0.00000    0.33989
multipole   212  211  210               0.15787
                                        0.00500    0.00000   -0.16579
                                        0.13325
                                        0.00000   -0.16718
                                       -0.06599    0.00000    0.03393
multipole   213    8  215              -0.08587
                                        0.20388    0.00000    0.25551
                                        0.08393
                                        0.00000   -0.69546
                                       -0.19022    0.00000    0.61153
multipole   214  213  215               0.08588
                                        0.03852    0.00000   -0.09266
                                        0.00460
                                        0.00000    0.05897
                                       -0.05937    0.00000   -0.06357
multipole   215  217  213              -0.01858
                                        0.22603    0.00000    0.33147
                                        0.15933
                                        0.00000   -0.60475
                                       -0.25158    0.00000    0.44542
multipole   216  215  216               0.08368
                                       -0.02778    0.00000   -0.06593
                                        0.07962
                                        0.00000    0.06169
                                        0.09443    0.00000   -0.14131
multipole   217  219  215              -0.01249
                                        0.17686    0.00000    0.14050
                                       -0.14798
                                        0.00000   -0.81782
                                       -0.11856    0.00000    0.96580
multipole   218  217  218               0.10208
                                       -0.04651    0.00000   -0.07618
                                       -0.01229
                                        0.00000   -0.00776
                                       -0.13903    0.00000    0.02005
multipole   219  217  220               0.10000
                                        0.03808    0.00000    0.22325
                                       -0.22293
                                        0.00000   -0.22476
                                        0.00964    0.00000    0.44769
multipole   220  219  217               0.20727
                                        0.03469    0.00000   -0.12995
                                       -0.04317
                                        0.00000   -0.08341
                                        0.18046    0.00000    0.12658
multipole   221  223  224              -0.19058
                                       -0.06827    0.00000    0.32517
                                       -0.13529
                                        0.00000   -0.15457
                                       -0.05672    0.00000    0.28986
multipole   222  221  223               0.07778
                                       -0.01715    0.00000   -0.10934
                                        0.05355
                                        0.00000    0.05042
                                       -0.02825    0.00000   -0.10397
multipole   223  224  221               0.68484
                                       -0.31279    0.00000    0.38571
                                        0.14709
                                        0.00000   -0.12737
                                        0.07287    0.00000   -0.01972
multipole   224  223  221              -0.72760
                                        0.09427    0.00000   -0.12507
                                       -0.42919
                                        0.00000    0.22290
                                        0.14497    0.00000    0.20629
multipole   225  226 -226              -0.30948
                                       -0.07657    0.00000    0.16462
                                        0.51606
                                        0.00000   -0.99579
                                       -0.25470    0.00000    0.47973
multipole   226  225  226               0.15474
                                       -0.06307    0.00000   -0.12148
                                        0.17239
                                        0.00000   -0.07008
                                       -0.03484    0.00000   -0.10231
multipole   227  229  228              -0.26804
                                        0.11821    0.00000    0.25875
                                        0.09964
                                        0.00000   -0.79062
                                       -0.22323    0.00000    0.69098
multipole   228  227  229               0.11919
                                        0.04956    0.00000   -0.09552
                                        0.10300
                                        0.00000   -0.00278
                                       -0.03565    0.00000   -0.10022
multipole   229  227  228              -0.00970
                                        0.06190    0.00000    0.26240
                                       -0.41933
                                        0.00000   -0.31279
                                        0.05123    0.00000    0.73212
multipole   230  229  227               0.05285
                                       -0.00218    0.00000   -0.11166
                                        0.02154
                                        0.00000    0.04194
                                       -0.06597    0.00000   -0.06348
multipole   231    2  232               0.07891
                                        0.00245    0.00000    0.24784
                                       -0.48795
                                        0.00000   -0.41413
                                        0.01190    0.00000    0.90208
multipole   231    8  232               0.11164
                                        0.27906    0.00000    0.31455
                                       -0.62632
                                        0.00000   -0.62629
                                        1.11853    0.00000    1.25261
multipole   231   48  232               0.11164
                                        0.27906    0.00000    0.31455
                                       -0.62632
                                        0.00000   -0.62629
                                        1.11853    0.00000    1.25261
multipole   232  231    2               0.26670
                                        0.06933    0.00000    0.04172
                                       -0.17154
                                        0.00000   -0.07853
                                        0.07867    0.00000    0.25007
multipole   232  231    8               0.21240
                                        0.02872    0.00000   -0.04889
                                       -0.22692
                                        0.00000    0.04052
                                        0.01035    0.00000    0.18640
multipole   232  231   48               0.21240
                                        0.00000    0.00000   -0.12490
                                        0.03622
                                        0.00000   -0.01437
                                        0.00000    0.00000   -0.02185
multipole   233  234 -234               1.02670
                                       -0.11581    0.00000   -0.04116
                                        0.39082
                                        0.00000   -0.32381
                                        0.38691    0.00000   -0.06701
multipole   234  233  234              -0.88956
                                        0.17189    0.00000   -0.22995
                                       -0.43857
                                        0.00000    0.34721
                                        0.38935    0.00000    0.09136
multipole   235  236    2               0.92336
                                       -0.20027    0.00000   -0.03856
                                        0.31193
                                        0.00000   -0.40576
                                        0.16648    0.00000    0.09383
multipole   235  236    8               0.83309
                                       -0.40323    0.00000    0.21520
                                        0.22432
                                        0.00000   -0.07777
                                        0.17848    0.00000   -0.14655
multipole   235  236   44               0.83309
                                       -0.40323    0.00000    0.21520
                                        0.22432
                                        0.00000   -0.07777
                                        0.17848    0.00000   -0.14655
multipole   236  235    2              -0.58827
                                        0.04503    0.00000   -0.13615
                                       -0.44066
                                        0.00000    0.28331
                                        0.05541    0.00000    0.15735
multipole   236  235    8              -0.61662
                                        0.11446    0.00000    0.01734
                                       -0.55506
                                        0.00000    0.08419
                                        0.16576    0.00000    0.47087
multipole   236  235   44              -0.61662
                                        0.11446    0.00000    0.01734
                                       -0.55506
                                        0.00000    0.08419
                                        0.16576    0.00000    0.47087
multipole   237  235  238              -0.46459
                                        0.15200    0.00000    0.02356
                                        0.42562
                                        0.00000   -0.67227
                                       -0.48912    0.00000    0.24665
multipole   238  237  235               0.24812
                                       -0.03635    0.00000   -0.05034
                                       -0.06195
                                        0.00000   -0.06380
                                       -0.05902    0.00000    0.12575
multipole   239  241  245               0.07219
                                        0.28913    0.00000    0.29472
                                        0.20132
                                        0.00000   -0.54286
                                        0.05155    0.00000    0.34154
multipole   240  239  240               0.23262
                                       -0.02789    0.00000   -0.04288
                                        0.00231
                                        0.00000    0.03311
                                       -0.02313    0.00000   -0.03542
multipole   241  239  242              -0.16428
                                        0.16276    0.20457    0.47894
                                       -0.28250
                                       -0.26022   -0.34943
                                        0.25272   -0.20303    0.63193
multipole   242  243  241               0.94625
                                        0.00000    0.00000    0.32359
                                        0.14171
                                        0.00000   -0.31370
                                        0.03825    0.00000    0.17199
multipole   243  242  241              -0.74794
                                        0.02444    0.00000   -0.13267
                                       -0.59982
                                        0.00000    0.24285
                                       -0.01926    0.00000    0.35697
multipole   244  241   55               0.12611
                                        0.01772    0.00000    0.00501
                                        0.02059
                                        0.00000   -0.01260
                                        0.01759    0.00000   -0.00799
multipole   245   57  239              -0.04567
                                        0.56406    0.00000    0.09733
                                        0.69818
                                        0.00000   -0.78665
                                       -0.18468    0.00000    0.08847
multipole   246  245  239               0.10559
                                       -0.01503    0.00000   -0.07122
                                        0.00434
                                        0.00000   -0.00248
                                       -0.01523    0.00000   -0.00186
multipole   247 -248 -248              -0.51966
                                        0.00000    0.00000    0.14279
                                        0.37928
                                        0.00000   -0.41809
                                        0.00000    0.00000    0.03881
multipole   248  247  248               0.25983
                                       -0.03859    0.00000   -0.05818
                                       -0.03673
                                        0.00000   -0.10739
                                       -0.00203    0.00000    0.14412
multipole   249    0    0               1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   250    0    0               1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   251    0    0               1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   252    0    0               1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   253    0    0               1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   254    0    0               2.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   255    0    0               2.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   256    0    0               2.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   257    0    0               2.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   258    0    0              -1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   259    0    0              -1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   260    0    0              -1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000
multipole   261    0    0              -1.00000
                                        0.00000    0.00000    0.00000
                                        0.00000
                                        0.00000    0.00000
                                        0.00000    0.00000    0.00000


      ########################################
      ##                                    ##
      ##  Dipole Polarizability Parameters  ##
      ##                                    ##
      ########################################


polarize      1          1.0730     0.3900      3    4    9   52  223  242
polarize      2          1.3340     0.3900      6
polarize      3          1.3340     0.3900      1    5    7   50  225  227
polarize      4          0.4960     0.3900      1
polarize      5          0.8370     0.3900      3
polarize      6          0.4960     0.3900      2
polarize      7          1.0730     0.3900      3    9   10   52  223  242
polarize      8          1.3340     0.3900     12
polarize      9          1.3340     0.3900      1    7   11   50  225  227
polarize     10          0.4960     0.3900      7
polarize     11          0.8370     0.3900      9
polarize     12          0.4960     0.3900      8   48
polarize     13          1.3340     0.3900     14
polarize     14          0.4960     0.3900     13
polarize     15          1.3340     0.3900     16   17
polarize     16          0.4960     0.3900     15
polarize     17          1.3340     0.3900     15   18
polarize     18          0.4960     0.3900     17
polarize     19          1.3340     0.3900     20   21
polarize     20          0.4960     0.3900     19
polarize     21          1.3340     0.3900     19   22   23
polarize     22          0.4960     0.3900     21
polarize     23          1.3340     0.3900     21   24
polarize     24          0.4960     0.3900     23
polarize     25          1.3340     0.3900     26   29
polarize     26          0.4960     0.3900     25
polarize     27          1.3340     0.3900     28   31
polarize     28          0.4960     0.3900     27
polarize     29          1.3340     0.3900     25   30
polarize     30          0.4960     0.3900     29
polarize     31          1.3340     0.3900     27   32
polarize     32          0.4960     0.3900     31
polarize     33          1.3340     0.3900     34   35
polarize     34          0.4960     0.3900     33
polarize     35          0.8370     0.3900     33   36
polarize     36          0.4960     0.3900     35
polarize     37          1.3340     0.3900     38   39
polarize     38          0.4960     0.3900     37
polarize     39          0.8370     0.3900     37   40
polarize     40          0.4960     0.3900     39
polarize     41          1.3340     0.3900     42
polarize     42          0.4960     0.3900     41
polarize     43          1.3340     0.3900     44   48   49
polarize     44          0.4960     0.3900     43
polarize     45          2.8000     0.3900     46
polarize     46          0.4960     0.3900     45
polarize     47          2.8000     0.3900     47
polarize     48          1.3340     0.3900     12   43
polarize     49          4.0000     0.3900     43
polarize     50          1.0730     0.3900      3    9   52   59  223  242
polarize     51          1.3340     0.3900     54
polarize     52          1.3340     0.3900      1    7   50   53  225  227
polarize     53          0.8370     0.3900     52
polarize     54          0.4960     0.3900     51
polarize     55          1.3340     0.3900     56
polarize     56          0.4960     0.3900     55
polarize     57          1.3340     0.3900     58
polarize     58          0.4960     0.3900     57
polarize     59          1.3340     0.3900     50   60
polarize     60          0.4960     0.3900     59
polarize     61          1.3340     0.3900     62   63
polarize     62          0.4960     0.3900     61
polarize     63          1.7500     0.3900     61   64
polarize     64          1.7500     0.3900     63   65   66
polarize     65          0.6960     0.3900     64
polarize     66          1.7500     0.3900     64   67   68
polarize     67          0.6960     0.3900     66
polarize     68          1.7500     0.3900     66   69
polarize     69          0.6960     0.3900     68
polarize     70          1.3340     0.3900     71   72
polarize     71          0.4960     0.3900     70
polarize     72          1.7500     0.3900     70   73
polarize     73          1.7500     0.3900     72   74   75
polarize     74          0.6960     0.3900     73
polarize     75          1.7500     0.3900     73   76   77
polarize     76          0.6960     0.3900     75
polarize     77          1.7500     0.3900     75   78
polarize     78          0.8370     0.3900     77   79
polarize     79          0.4960     0.3900     78
polarize     80          1.3340     0.3900     81
polarize     81          0.4960     0.3900     80
polarize     82          1.7500     0.3900     83
polarize     83          1.7500     0.3900     82   84   85
polarize     84          0.6960     0.3900     83
polarize     85          1.7500     0.3900     83   86   87
polarize     86          0.6960     0.3900     85
polarize     87          1.7500     0.3900     85   88
polarize     88          0.8370     0.3900     87
polarize     89          1.3340     0.3900     90
polarize     90          0.4960     0.3900     89
polarize     91          1.7500     0.3900     92   94
polarize     92          1.7500     0.3900     91   93   95
polarize     93          0.6960     0.3900     92
polarize     94          1.7500     0.3900     91   97   98
polarize     95          1.0730     0.3900     92   96   97
polarize     96          0.4960     0.3900     95
polarize     97          1.7500     0.3900     94   95  100
polarize     98          1.7500     0.3900     94   99  102
polarize     99          0.6960     0.3900     98
polarize    100          1.7500     0.3900     97  101  104
polarize    101          0.6960     0.3900    100
polarize    102          1.7500     0.3900     98  103  104
polarize    103          0.6960     0.3900    102
polarize    104          1.7500     0.3900    100  102  105
polarize    105          0.6960     0.3900    104
polarize    106          1.3340     0.3900    107  108
polarize    107          0.4960     0.3900    106
polarize    108          1.7500     0.3900    106  109  111
polarize    109          1.5000     0.3900    108  110  113
polarize    110          0.6960     0.3900    109
polarize    111          1.7500     0.3900    108  112  115
polarize    112          0.6960     0.3900    111
polarize    113          1.7500     0.3900    109  114  115
polarize    114          0.6960     0.3900    113
polarize    115          1.5000     0.3900    111  113  116
polarize    116          0.6960     0.3900    115
polarize    117          1.3340     0.3900    118  119
polarize    118          0.4960     0.3900    117
polarize    119          1.7500     0.3900    117  120  122
polarize    120          1.5000     0.3900    119  121  122  124
polarize    121          0.6960     0.3900    120
polarize    122          1.7500     0.3900    119  120  123  126
polarize    123          0.6960     0.3900    122
polarize    124          1.7500     0.3900    120  125  126
polarize    125          0.6960     0.3900    124
polarize    126          1.5000     0.3900    122  124
polarize    127          1.3340     0.3900    128  129
polarize    128          0.4960     0.3900    127
polarize    129          1.7500     0.3900    127  130  131
polarize    130          1.5000     0.3900    129  133
polarize    131          1.7500     0.3900    129  132  135
polarize    132          0.6960     0.3900    131
polarize    133          1.7500     0.3900    130  134  135
polarize    134          0.6960     0.3900    133
polarize    135          1.5000     0.3900    131  133  136
polarize    136          0.6960     0.3900    135
polarize    137          1.3340     0.3900    138  139
polarize    138          0.4960     0.3900    137
polarize    139          1.3340     0.3900    137  140
polarize    140          1.2000     0.3900    139
polarize    141          1.3340     0.3900    142  143
polarize    142          0.4960     0.3900    141
polarize    143          1.3340     0.3900    141  144  145
polarize    144          0.8370     0.3900    143
polarize    145          0.8370     0.3900    143  146
polarize    146          0.4960     0.3900    145
polarize    147          1.3340     0.3900    148  149
polarize    148          0.4960     0.3900    147
polarize    149          1.3340     0.3900    147  150  151
polarize    150          0.8340     0.3900    149
polarize    151          1.0730     0.3900    149  152
polarize    152          0.4960     0.3900    151
polarize    153          1.3340     0.3900    154
polarize    154          0.4960     0.3900    153
polarize    155          1.3340     0.3900    156  157
polarize    156          0.4960     0.3900    155
polarize    157          1.3340     0.3900    155  158
polarize    158          1.2000     0.3900    157
polarize    159          1.3340     0.3900    160
polarize    160          0.4960     0.3900    159
polarize    161          1.3340     0.3900    162  163
polarize    162          0.4960     0.3900    161
polarize    163          1.3340     0.3900    161  164  165
polarize    164          0.8370     0.3900    163
polarize    165          0.8370     0.3900    163  166
polarize    166          0.4960     0.3900    165
polarize    167          1.3340     0.3900    168
polarize    168          0.4960     0.3900    167
polarize    169          1.3340     0.3900    170  171
polarize    170          0.4960     0.3900    169
polarize    171          1.3340     0.3900    169  172  173
polarize    172          0.8340     0.3900    171
polarize    173          1.0730     0.3900    171  174
polarize    174          0.4960     0.3900    173
polarize    175          1.3340     0.3900    176
polarize    176          0.4960     0.3900    175
polarize    177          1.3340     0.3900    178  179
polarize    178          0.4960     0.3900    177
polarize    179          3.3000     0.3900    177  180
polarize    180          1.3340     0.3900    179  181
polarize    181          0.4960     0.3900    180
polarize    182          1.3340     0.3900    183  184
polarize    183          0.4960     0.3900    182
polarize    184          1.3340     0.3900    182  185  186
polarize    185          0.4960     0.3900    184
polarize    186          1.3340     0.3900    184  187
polarize    187          0.4960     0.3900    186
polarize    188          1.3340     0.3900    189  190
polarize    189          0.4960     0.3900    188
polarize    190          1.0730     0.3900    188  191
polarize    191          0.4960     0.3900    190
polarize    192          1.3340     0.3900    193  194
polarize    193          0.4960     0.3900    192
polarize    194          1.3340     0.3900    192  195  196
polarize    195          0.4960     0.3900    194
polarize    196          1.3340     0.3900    194  197
polarize    197          0.4960     0.3900    196
polarize    198          1.3340     0.3900    199  200
polarize    199          0.4960     0.3900    198
polarize    200          1.0730     0.3900    198  201
polarize    201          0.4960     0.3900    200
polarize    202          1.3340     0.3900    203  204
polarize    203          0.4960     0.3900    202
polarize    204          1.3340     0.3900    202  205  206
polarize    205          0.4960     0.3900    204
polarize    206          1.3340     0.3900    204  207
polarize    207          0.4960     0.3900    206
polarize    208          1.0730     0.3900    209  210
polarize    209          0.4960     0.3900    208
polarize    210          1.3340     0.3900    208  211
polarize    211          1.0730     0.3900    210  212
polarize    212          0.4960     0.3900    211
polarize    213          1.3340     0.3900    214  215
polarize    214          0.4960     0.3900    213
polarize    215          1.3340     0.3900    213  216
polarize    216          0.4960     0.3900    215
polarize    217          1.3340     0.3900    218  219
polarize    218          0.4960     0.3900    217
polarize    219          1.0730     0.3900    217  220
polarize    220          0.4960     0.3900    219
polarize    221          1.3340     0.3900    222
polarize    222          0.4960     0.3900    221
polarize    223          1.3340     0.3900      1    7   50  224
polarize    224          0.8370     0.3900    223
polarize    225          1.0730     0.3900      3    9   52  226  242
polarize    226          0.4960     0.3900    225
polarize    227          1.0730     0.3900      3    9   52  228  242
polarize    228          0.4960     0.3900    227
polarize    229          1.3340     0.3900    230
polarize    230          0.4960     0.3900    229
polarize    231          1.0730     0.3900    232
polarize    232          0.4960     0.3900    231
polarize    233          1.3340     0.3900    234
polarize    234          0.8370     0.3900    233
polarize    235          1.3340     0.3900    236  237
polarize    236          0.8370     0.3900    235
polarize    237          0.8370     0.3900    235  238
polarize    238          0.4960     0.3900    237
polarize    239          1.0730     0.3900    240
polarize    240          0.4960     0.3900    239
polarize    241          1.3340     0.3900    244
polarize    242          1.3340     0.3900      1    7   50  225  227  243
polarize    243          0.8370     0.3900    242
polarize    244          0.4960     0.3900    241
polarize    245          1.3340     0.3900    246
polarize    246          0.4960     0.3900    245
polarize    247          0.8370     0.3900    248
polarize    248          0.4960     0.3900    247
polarize    249          0.0280     0.3900
polarize    250          0.1200     0.3900
polarize    251          0.7800     0.3900
polarize    252          1.3500     0.3900
polarize    253          2.2600     0.3900
polarize    254          0.0100     0.0650
polarize    255          0.0800     0.1150
polarize    256          0.5500     0.1800
polarize    257          0.2600     0.2096
polarize    258          1.3500     0.3900
polarize    259          4.0000     0.3900
polarize    260          5.6500     0.3900
polarize    261          7.2500     0.3900


      ########################################
      ##                                    ##
      ##  Biopolymer Atom Type Conversions  ##
      ##                                    ##
      ########################################


biotype       1    N       "Glycine"                           1
biotype       2    CA      "Glycine"                           2
biotype       3    C       "Glycine"                           3
biotype       4    HN      "Glycine"                           4
biotype       5    O       "Glycine"                           5
biotype       6    HA      "Glycine"                           6
biotype       7    N       "Alanine"                           7
biotype       8    CA      "Alanine"                           8
biotype       9    C       "Alanine"                           9
biotype      10    HN      "Alanine"                          10
biotype      11    O       "Alanine"                          11
biotype      12    HA      "Alanine"                          12
biotype      13    CB      "Alanine"                          13
biotype      14    HB      "Alanine"                          14
biotype      15    N       "Valine"                            7
biotype      16    CA      "Valine"                            8
biotype      17    C       "Valine"                            9
biotype      18    HN      "Valine"                           10
biotype      19    O       "Valine"                           11
biotype      20    HA      "Valine"                           12
biotype      21    CB      "Valine"                           15
biotype      22    HB      "Valine"                           16
biotype      23    CG1     "Valine"                           17
biotype      24    HG1     "Valine"                           18
biotype      25    CG2     "Valine"                           17
biotype      26    HG2     "Valine"                           18
biotype      27    N       "Leucine"                           7
biotype      28    CA      "Leucine"                           8
biotype      29    C       "Leucine"                           9
biotype      30    HN      "Leucine"                          10
biotype      31    O       "Leucine"                          11
biotype      32    HA      "Leucine"                          12
biotype      33    CB      "Leucine"                          19
biotype      34    HB      "Leucine"                          20
biotype      35    CG      "Leucine"                          21
biotype      36    HG      "Leucine"                          22
biotype      37    CD1     "Leucine"                          23
biotype      38    HD1     "Leucine"                          24
biotype      39    CD2     "Leucine"                          23
biotype      40    HD2     "Leucine"                          24
biotype      41    N       "Isoleucine"                        7
biotype      42    CA      "Isoleucine"                        8
biotype      43    C       "Isoleucine"                        9
biotype      44    HN      "Isoleucine"                       10
biotype      45    O       "Isoleucine"                       11
biotype      46    HA      "Isoleucine"                       12
biotype      47    CB      "Isoleucine"                       25
biotype      48    HB      "Isoleucine"                       26
biotype      49    CG1     "Isoleucine"                       27
biotype      50    HG1     "Isoleucine"                       28
biotype      51    CG2     "Isoleucine"                       29
biotype      52    HG2     "Isoleucine"                       30
biotype      53    CD      "Isoleucine"                       31
biotype      54    HD      "Isoleucine"                       32
biotype      55    N       "Serine"                            7
biotype      56    CA      "Serine"                            8
biotype      57    C       "Serine"                            9
biotype      58    HN      "Serine"                           10
biotype      59    O       "Serine"                           11
biotype      60    HA      "Serine"                           12
biotype      61    CB      "Serine"                           33
biotype      62    HB      "Serine"                           34
biotype      63    OG      "Serine"                           35
biotype      64    HG      "Serine"                           36
biotype      65    N       "Threonine"                         7
biotype      66    CA      "Threonine"                         8
biotype      67    C       "Threonine"                         9
biotype      68    HN      "Threonine"                        10
biotype      69    O       "Threonine"                        11
biotype      70    HA      "Threonine"                        12
biotype      71    CB      "Threonine"                        37
biotype      72    HB      "Threonine"                        38
biotype      73    OG1     "Threonine"                        39
biotype      74    HG1     "Threonine"                        40
biotype      75    CG2     "Threonine"                        41
biotype      76    HG2     "Threonine"                        42
biotype      77    N       "Cysteine (SH)"                     7
biotype      78    CA      "Cysteine (SH)"                     8
biotype      79    C       "Cysteine (SH)"                     9
biotype      80    HN      "Cysteine (SH)"                    10
biotype      81    O       "Cysteine (SH)"                    11
biotype      82    HA      "Cysteine (SH)"                    12
biotype      83    CB      "Cysteine (SH)"                    43
biotype      84    HB      "Cysteine (SH)"                    44
biotype      85    SG      "Cysteine (SH)"                    45
biotype      86    HG      "Cysteine (SH)"                    46
biotype      87    N       "Cystine (SS)"                      7
biotype      88    CA      "Cystine (SS)"                      8
biotype      89    C       "Cystine (SS)"                      9
biotype      90    HN      "Cystine (SS)"                     10
biotype      91    O       "Cystine (SS)"                     11
biotype      92    HA      "Cystine (SS)"                     12
biotype      93    CB      "Cystine (SS)"                     43
biotype      94    HB      "Cystine (SS)"                     44
biotype      95    SG      "Cystine (SS)"                     47
biotype      96    N       "Cysteine (S-)"                     7
biotype      97    CA      "Cysteine (S-)"                    48
biotype      98    C       "Cysteine (S-)"                     9
biotype      99    HN      "Cysteine (S-)"                    10
biotype     100    O       "Cysteine (S-)"                    11
biotype     101    HA      "Cysteine (S-)"                    12
biotype     102    CB      "Cysteine (S-)"                    43
biotype     103    HB      "Cysteine (S-)"                    44
biotype     104    SG      "Cysteine (S-)"                    49
biotype     105    N       "Proline"                          50
biotype     106    CA      "Proline"                          51
biotype     107    C       "Proline"                          52
biotype     108    O       "Proline"                          53
biotype     109    HA      "Proline"                          54
biotype     110    CB      "Proline"                          55
biotype     111    HB      "Proline"                          56
biotype     112    CG      "Proline"                          57
biotype     113    HG      "Proline"                          58
biotype     114    CD      "Proline"                          59
biotype     115    HD      "Proline"                          60
biotype     116    N       "Phenylalanine"                     7
biotype     117    CA      "Phenylalanine"                     8
biotype     118    C       "Phenylalanine"                     9
biotype     119    HN      "Phenylalanine"                    10
biotype     120    O       "Phenylalanine"                    11
biotype     121    HA      "Phenylalanine"                    12
biotype     122    CB      "Phenylalanine"                    61
biotype     123    HB      "Phenylalanine"                    62
biotype     124    CG      "Phenylalanine"                    63
biotype     125    CD      "Phenylalanine"                    64
biotype     126    HD      "Phenylalanine"                    65
biotype     127    CE      "Phenylalanine"                    66
biotype     128    HE      "Phenylalanine"                    67
biotype     129    CZ      "Phenylalanine"                    68
biotype     130    HZ      "Phenylalanine"                    69
biotype     131    N       "Tyrosine"                          7
biotype     132    CA      "Tyrosine"                          8
biotype     133    C       "Tyrosine"                          9
biotype     134    HN      "Tyrosine"                         10
biotype     135    O       "Tyrosine"                         11
biotype     136    HA      "Tyrosine"                         12
biotype     137    CB      "Tyrosine"                         70
biotype     138    HB      "Tyrosine"                         71
biotype     139    CG      "Tyrosine"                         72
biotype     140    CD      "Tyrosine"                         73
biotype     141    HD      "Tyrosine"                         74
biotype     142    CE      "Tyrosine"                         75
biotype     143    HE      "Tyrosine"                         76
biotype     144    CZ      "Tyrosine"                         77
biotype     145    OH      "Tyrosine"                         78
biotype     146    HH      "Tyrosine"                         79
biotype     147    N       "Tyrosine (O-)"                     7
biotype     148    CA      "Tyrosine (O-)"                     8
biotype     149    C       "Tyrosine (O-)"                     9
biotype     150    HN      "Tyrosine (O-)"                    10
biotype     151    O       "Tyrosine (O-)"                    11
biotype     152    HA      "Tyrosine (O-)"                    12
biotype     153    CB      "Tyrosine (O-)"                    80
biotype     154    HB      "Tyrosine (O-)"                    81
biotype     155    CG      "Tyrosine (O-)"                    82
biotype     156    CD      "Tyrosine (O-)"                    83
biotype     157    HD      "Tyrosine (O-)"                    84
biotype     158    CE      "Tyrosine (O-)"                    85
biotype     159    HE      "Tyrosine (O-)"                    86
biotype     160    CZ      "Tyrosine (O-)"                    87
biotype     161    OH      "Tyrosine (O-)"                    88
biotype     162    N       "Tryptophan"                        7
biotype     163    CA      "Tryptophan"                        8
biotype     164    C       "Tryptophan"                        9
biotype     165    HN      "Tryptophan"                       10
biotype     166    O       "Tryptophan"                       11
biotype     167    HA      "Tryptophan"                       12
biotype     168    CB      "Tryptophan"                       89
biotype     169    HB      "Tryptophan"                       90
biotype     170    CG      "Tryptophan"                       91
biotype     171    CD1     "Tryptophan"                       92
biotype     172    HD1     "Tryptophan"                       93
biotype     173    CD2     "Tryptophan"                       94
biotype     174    NE1     "Tryptophan"                       95
biotype     175    HE1     "Tryptophan"                       96
biotype     176    CE2     "Tryptophan"                       97
biotype     177    CE3     "Tryptophan"                       98
biotype     178    HE3     "Tryptophan"                       99
biotype     179    CZ2     "Tryptophan"                      100
biotype     180    HZ2     "Tryptophan"                      101
biotype     181    CZ3     "Tryptophan"                      102
biotype     182    HZ3     "Tryptophan"                      103
biotype     183    CH2     "Tryptophan"                      104
biotype     184    HH2     "Tryptophan"                      105
biotype     185    N       "Histidine (+)"                     7
biotype     186    CA      "Histidine (+)"                     8
biotype     187    C       "Histidine (+)"                     9
biotype     188    HN      "Histidine (+)"                    10
biotype     189    O       "Histidine (+)"                    11
biotype     190    HA      "Histidine (+)"                    12
biotype     191    CB      "Histidine (+)"                   106
biotype     192    HB      "Histidine (+)"                   107
biotype     193    CG      "Histidine (+)"                   108
biotype     194    ND1     "Histidine (+)"                   109
biotype     195    HD1     "Histidine (+)"                   110
biotype     196    CD2     "Histidine (+)"                   111
biotype     197    HD2     "Histidine (+)"                   112
biotype     198    CE1     "Histidine (+)"                   113
biotype     199    HE1     "Histidine (+)"                   114
biotype     200    NE2     "Histidine (+)"                   115
biotype     201    HE2     "Histidine (+)"                   116
biotype     202    N       "Histidine (HD)"                    7
biotype     203    CA      "Histidine (HD)"                    8
biotype     204    C       "Histidine (HD)"                    9
biotype     205    HN      "Histidine (HD)"                   10
biotype     206    O       "Histidine (HD)"                   11
biotype     207    HA      "Histidine (HD)"                   12
biotype     208    CB      "Histidine (HD)"                  117
biotype     209    HB      "Histidine (HD)"                  118
biotype     210    CG      "Histidine (HD)"                  119
biotype     211    ND1     "Histidine (HD)"                  120
biotype     212    HD1     "Histidine (HD)"                  121
biotype     213    CD2     "Histidine (HD)"                  122
biotype     214    HD2     "Histidine (HD)"                  123
biotype     215    CE1     "Histidine (HD)"                  124
biotype     216    HE1     "Histidine (HD)"                  125
biotype     217    NE2     "Histidine (HD)"                  126
biotype     218    N       "Histidine (HE)"                    7
biotype     219    CA      "Histidine (HE)"                    8
biotype     220    C       "Histidine (HE)"                    9
biotype     221    HN      "Histidine (HE)"                   10
biotype     222    O       "Histidine (HE)"                   11
biotype     223    HA      "Histidine (HE)"                   12
biotype     224    CB      "Histidine (HE)"                  127
biotype     225    HB      "Histidine (HE)"                  128
biotype     226    CG      "Histidine (HE)"                  129
biotype     227    ND1     "Histidine (HE)"                  130
biotype     228    CD2     "Histidine (HE)"                  131
biotype     229    HD2     "Histidine (HE)"                  132
biotype     230    CE1     "Histidine (HE)"                  133
biotype     231    HE1     "Histidine (HE)"                  134
biotype     232    NE2     "Histidine (HE)"                  135
biotype     233    HE2     "Histidine (HE)"                  136
biotype     234    N       "Aspartic Acid"                     7
biotype     235    CA      "Aspartic Acid"                     8
biotype     236    C       "Aspartic Acid"                     9
biotype     237    HN      "Aspartic Acid"                    10
biotype     238    O       "Aspartic Acid"                    11
biotype     239    HA      "Aspartic Acid"                    12
biotype     240    CB      "Aspartic Acid"                   137
biotype     241    HB      "Aspartic Acid"                   138
biotype     242    CG      "Aspartic Acid"                   139
biotype     243    OD      "Aspartic Acid"                   140
biotype     244    N       "Aspartic Acid (COOH)"              7
biotype     245    CA      "Aspartic Acid (COOH)"              8
biotype     246    C       "Aspartic Acid (COOH)"              9
biotype     247    HN      "Aspartic Acid (COOH)"             10
biotype     248    O       "Aspartic Acid (COOH)"             11
biotype     249    HA      "Aspartic Acid (COOH)"             12
biotype     250    CB      "Aspartic Acid (COOH)"            141
biotype     251    HB      "Aspartic Acid (COOH)"            142
biotype     252    CG      "Aspartic Acid (COOH)"            143
biotype     253    OD1     "Aspartic Acid (COOH)"            144
biotype     254    OD2     "Aspartic Acid (COOH)"            145
biotype     255    HD2     "Aspartic Acid (COOH)"            146
biotype     256    N       "Asparagine"                        7
biotype     257    CA      "Asparagine"                        8
biotype     258    C       "Asparagine"                        9
biotype     259    HN      "Asparagine"                       10
biotype     260    O       "Asparagine"                       11
biotype     261    HA      "Asparagine"                       12
biotype     262    CB      "Asparagine"                      147
biotype     263    HB      "Asparagine"                      148
biotype     264    CG      "Asparagine"                      149
biotype     265    OD1     "Asparagine"                      150
biotype     266    ND2     "Asparagine"                      151
biotype     267    HD2     "Asparagine"                      152
biotype     268    N       "Glutamic Acid"                     7
biotype     269    CA      "Glutamic Acid"                     8
biotype     270    C       "Glutamic Acid"                     9
biotype     271    HN      "Glutamic Acid"                    10
biotype     272    O       "Glutamic Acid"                    11
biotype     273    HA      "Glutamic Acid"                    12
biotype     274    CB      "Glutamic Acid"                   153
biotype     275    HB      "Glutamic Acid"                   154
biotype     276    CG      "Glutamic Acid"                   155
biotype     277    HG      "Glutamic Acid"                   156
biotype     278    CD      "Glutamic Acid"                   157
biotype     279    OE      "Glutamic Acid"                   158
biotype     280    N       "Glutamic Acid (COOH)"              7
biotype     281    CA      "Glutamic Acid (COOH)"              8
biotype     282    C       "Glutamic Acid (COOH)"              9
biotype     283    HN      "Glutamic Acid (COOH)"             10
biotype     284    O       "Glutamic Acid (COOH)"             11
biotype     285    HA      "Glutamic Acid (COOH)"             12
biotype     286    CB      "Glutamic Acid (COOH)"            159
biotype     287    HB      "Glutamic Acid (COOH)"            160
biotype     288    CG      "Glutamic Acid (COOH)"            161
biotype     289    HG      "Glutamic Acid (COOH)"            162
biotype     290    CD      "Glutamic Acid (COOH)"            163
biotype     291    OE1     "Glutamic Acid (COOH)"            164
biotype     292    OE2     "Glutamic Acid (COOH)"            165
biotype     293    HE2     "Glutamic Acid (COOH)"            166
biotype     294    N       "Glutamine"                         7
biotype     295    CA      "Glutamine"                         8
biotype     296    C       "Glutamine"                         9
biotype     297    HN      "Glutamine"                        10
biotype     298    O       "Glutamine"                        11
biotype     299    HA      "Glutamine"                        12
biotype     300    CB      "Glutamine"                       167
biotype     301    HB      "Glutamine"                       168
biotype     302    CG      "Glutamine"                       169
biotype     303    HG      "Glutamine"                       170
biotype     304    CD      "Glutamine"                       171
biotype     305    OE1     "Glutamine"                       172
biotype     306    NE2     "Glutamine"                       173
biotype     307    HE2     "Glutamine"                       174
biotype     308    N       "Methionine"                        7
biotype     309    CA      "Methionine"                        8
biotype     310    C       "Methionine"                        9
biotype     311    HN      "Methionine"                       10
biotype     312    O       "Methionine"                       11
biotype     313    HA      "Methionine"                       12
biotype     314    CB      "Methionine"                      175
biotype     315    HB      "Methionine"                      176
biotype     316    CG      "Methionine"                      177
biotype     317    HG      "Methionine"                      178
biotype     318    SD      "Methionine"                      179
biotype     319    CE      "Methionine"                      180
biotype     320    HE      "Methionine"                      181
biotype     321    N       "Lysine"                            7
biotype     322    CA      "Lysine"                            8
biotype     323    C       "Lysine"                            9
biotype     324    HN      "Lysine"                           10
biotype     325    O       "Lysine"                           11
biotype     326    HA      "Lysine"                           12
biotype     327    CB      "Lysine"                          182
biotype     328    HB      "Lysine"                          183
biotype     329    CG      "Lysine"                          184
biotype     330    HG      "Lysine"                          185
biotype     331    CD      "Lysine"                          186
biotype     332    HD      "Lysine"                          187
biotype     333    CE      "Lysine"                          188
biotype     334    HE      "Lysine"                          189
biotype     335    NZ      "Lysine"                          190
biotype     336    HZ      "Lysine"                          191
biotype     337    N       "Lysine (NH2)"                      7
biotype     338    CA      "Lysine (NH2)"                      8
biotype     339    C       "Lysine (NH2)"                      9
biotype     340    HN      "Lysine (NH2)"                     10
biotype     341    O       "Lysine (NH2)"                     11
biotype     342    HA      "Lysine (NH2)"                     12
biotype     343    CB      "Lysine (NH2)"                    192
biotype     344    HB      "Lysine (NH2)"                    193
biotype     345    CG      "Lysine (NH2)"                    194
biotype     346    HG      "Lysine (NH2)"                    195
biotype     347    CD      "Lysine (NH2)"                    196
biotype     348    HD      "Lysine (NH2)"                    197
biotype     349    CE      "Lysine (NH2)"                    198
biotype     350    HE      "Lysine (NH2)"                    199
biotype     351    NZ      "Lysine (NH2)"                    200
biotype     352    HZ      "Lysine (NH2)"                    201
biotype     353    N       "Arginine"                          7
biotype     354    CA      "Arginine"                          8
biotype     355    C       "Arginine"                          9
biotype     356    HN      "Arginine"                         10
biotype     357    O       "Arginine"                         11
biotype     358    HA      "Arginine"                         12
biotype     359    CB      "Arginine"                        202
biotype     360    HB      "Arginine"                        203
biotype     361    CG      "Arginine"                        204
biotype     362    HG      "Arginine"                        205
biotype     363    CD      "Arginine"                        206
biotype     364    HD      "Arginine"                        207
biotype     365    NE      "Arginine"                        208
biotype     366    HE      "Arginine"                        209
biotype     367    CZ      "Arginine"                        210
biotype     368    NH      "Arginine"                        211
biotype     369    HH      "Arginine"                        212
biotype     370    N       "Ornithine"                         7
biotype     371    CA      "Ornithine"                         8
biotype     372    C       "Ornithine"                         9
biotype     373    HN      "Ornithine"                        10
biotype     374    O       "Ornithine"                        11
biotype     375    HA      "Ornithine"                        12
biotype     376    CB      "Ornithine"                       213
biotype     377    HB      "Ornithine"                       214
biotype     378    CG      "Ornithine"                       215
biotype     379    HG      "Ornithine"                       216
biotype     380    CD      "Ornithine"                       217
biotype     381    HD      "Ornithine"                       218
biotype     382    NE      "Ornithine"                       219
biotype     383    HE      "Ornithine"                       220
biotype     384    N       "MethylAlanine (AIB)"              -1
biotype     385    CA      "MethylAlanine (AIB)"              -1
biotype     386    C       "MethylAlanine (AIB)"              -1
biotype     387    HN      "MethylAlanine (AIB)"              -1
biotype     388    O       "MethylAlanine (AIB)"              -1
biotype     389    CB      "MethylAlanine (AIB)"              -1
biotype     390    HB      "MethylAlanine (AIB)"              -1
biotype     391    N       "Pyroglutamic Acid"                -1
biotype     392    CA      "Pyroglutamic Acid"                -1
biotype     393    C       "Pyroglutamic Acid"                -1
biotype     394    HN      "Pyroglutamic Acid"                -1
biotype     395    O       "Pyroglutamic Acid"                -1
biotype     396    HA      "Pyroglutamic Acid"                -1
biotype     397    CB      "Pyroglutamic Acid"                -1
biotype     398    HB      "Pyroglutamic Acid"                -1
biotype     399    CG      "Pyroglutamic Acid"                -1
biotype     400    HG      "Pyroglutamic Acid"                -1
biotype     401    CD      "Pyroglutamic Acid"                -1
biotype     402    OE      "Pyroglutamic Acid"                -1
biotype     403    N       "N-Terminal GLY"                  231
biotype     404    CA      "N-Terminal GLY"                    2
biotype     405    C       "N-Terminal GLY"                    3
biotype     406    HN      "N-Terminal GLY"                  232
biotype     407    O       "N-Terminal GLY"                    5
biotype     408    HA      "N-Terminal GLY"                    6
biotype     409    N       "N-Terminal ALA"                  231
biotype     410    CA      "N-Terminal ALA"                    8
biotype     411    C       "N-Terminal ALA"                    9
biotype     412    HN      "N-Terminal ALA"                  232
biotype     413    O       "N-Terminal ALA"                   11
biotype     414    HA      "N-Terminal ALA"                   12
biotype     415    N       "N-Terminal VAL"                  231
biotype     416    CA      "N-Terminal VAL"                    8
biotype     417    C       "N-Terminal VAL"                    9
biotype     418    HN      "N-Terminal VAL"                  232
biotype     419    O       "N-Terminal VAL"                   11
biotype     420    HA      "N-Terminal VAL"                   12
biotype     421    N       "N-Terminal LEU"                  231
biotype     422    CA      "N-Terminal LEU"                    8
biotype     423    C       "N-Terminal LEU"                    9
biotype     424    HN      "N-Terminal LEU"                  232
biotype     425    O       "N-Terminal LEU"                   11
biotype     426    HA      "N-Terminal LEU"                   12
biotype     427    N       "N-Terminal ILE"                  231
biotype     428    CA      "N-Terminal ILE"                    8
biotype     429    C       "N-Terminal ILE"                    9
biotype     430    HN      "N-Terminal ILE"                  232
biotype     431    O       "N-Terminal ILE"                   11
biotype     432    HA      "N-Terminal ILE"                   12
biotype     433    N       "N-Terminal SER"                  231
biotype     434    CA      "N-Terminal SER"                    8
biotype     435    C       "N-Terminal SER"                    9
biotype     436    HN      "N-Terminal SER"                  232
biotype     437    O       "N-Terminal SER"                   11
biotype     438    HA      "N-Terminal SER"                   12
biotype     439    N       "N-Terminal THR"                  231
biotype     440    CA      "N-Terminal THR"                    8
biotype     441    C       "N-Terminal THR"                    9
biotype     442    HN      "N-Terminal THR"                  232
biotype     443    O       "N-Terminal THR"                   11
biotype     444    HA      "N-Terminal THR"                   12
biotype     445    N       "N-Terminal CYS (SH)"             231
biotype     446    CA      "N-Terminal CYS (SH)"               8
biotype     447    C       "N-Terminal CYS (SH)"               9
biotype     448    HN      "N-Terminal CYS (SH)"             232
biotype     449    O       "N-Terminal CYS (SH)"              11
biotype     450    HA      "N-Terminal CYS (SH)"              12
biotype     451    N       "N-Terminal CYX (SS)"             231
biotype     452    CA      "N-Terminal CYX (SS)"               8
biotype     453    C       "N-Terminal CYX (SS)"               9
biotype     454    HN      "N-Terminal CYX (SS)"             232
biotype     455    O       "N-Terminal CYX (SS)"              11
biotype     456    HA      "N-Terminal CYX (SS)"              12
biotype     457    N       "N-Terminal CYD (S-)"             231
biotype     458    CA      "N-Terminal CYD (S-)"              48
biotype     459    C       "N-Terminal CYD (S-)"               9
biotype     460    HN      "N-Terminal CYD (S-)"             232
biotype     461    O       "N-Terminal CYD (S-)"              11
biotype     462    HA      "N-Terminal CYD (S-)"              12
biotype     463    N       "N-Terminal PRO"                  239
biotype     464    CA      "N-Terminal PRO"                  241
biotype     465    C       "N-Terminal PRO"                  242
biotype     466    HN      "N-Terminal PRO"                  240
biotype     467    O       "N-Terminal PRO"                  243
biotype     468    HA      "N-Terminal PRO"                  244
biotype     469    CD      "N-Terminal PRO"                  245
biotype     470    HD      "N-Terminal PRO"                  246
biotype     471    N       "N-Terminal PHE"                  231
biotype     472    CA      "N-Terminal PHE"                    8
biotype     473    C       "N-Terminal PHE"                    9
biotype     474    HN      "N-Terminal PHE"                  232
biotype     475    O       "N-Terminal PHE"                   11
biotype     476    HA      "N-Terminal PHE"                   12
biotype     477    N       "N-Terminal TYR"                  231
biotype     478    CA      "N-Terminal TYR"                    8
biotype     479    C       "N-Terminal TYR"                    9
biotype     480    HN      "N-Terminal TYR"                  232
biotype     481    O       "N-Terminal TYR"                   11
biotype     482    HA      "N-Terminal TYR"                   12
biotype     483    N       "N-Terminal TYD (O-)"             231
biotype     484    CA      "N-Terminal TYD (O-)"               8
biotype     485    C       "N-Terminal TYD (O-)"               9
biotype     486    HN      "N-Terminal TYD (O-)"             232
biotype     487    O       "N-Terminal TYD (O-)"              11
biotype     488    HA      "N-Terminal TYD (O-)"              12
biotype     489    N       "N-Terminal TRP"                  231
biotype     490    CA      "N-Terminal TRP"                    8
biotype     491    C       "N-Terminal TRP"                    9
biotype     492    HN      "N-Terminal TRP"                  232
biotype     493    O       "N-Terminal TRP"                   11
biotype     494    HA      "N-Terminal TRP"                   12
biotype     495    N       "N-Terminal HIS (+)"              231
biotype     496    CA      "N-Terminal HIS (+)"                8
biotype     497    C       "N-Terminal HIS (+)"                9
biotype     498    HN      "N-Terminal HIS (+)"              232
biotype     499    O       "N-Terminal HIS (+)"               11
biotype     500    HA      "N-Terminal HIS (+)"               12
biotype     501    N       "N-Terminal HIS (HD)"             231
biotype     502    CA      "N-Terminal HIS (HD)"               8
biotype     503    C       "N-Terminal HIS (HD)"               9
biotype     504    HN      "N-Terminal HIS (HD)"             232
biotype     505    O       "N-Terminal HIS (HD)"              11
biotype     506    HA      "N-Terminal HIS (HD)"              12
biotype     507    N       "N-Terminal HIS (HE)"             231
biotype     508    CA      "N-Terminal HIS (HE)"               8
biotype     509    C       "N-Terminal HIS (HE)"               9
biotype     510    HN      "N-Terminal HIS (HE)"             232
biotype     511    O       "N-Terminal HIS (HE)"              11
biotype     512    HA      "N-Terminal HIS (HE)"              12
biotype     513    N       "N-Terminal ASP"                  231
biotype     514    CA      "N-Terminal ASP"                    8
biotype     515    C       "N-Terminal ASP"                    9
biotype     516    HN      "N-Terminal ASP"                  232
biotype     517    O       "N-Terminal ASP"                   11
biotype     518    HA      "N-Terminal ASP"                   12
biotype     519    N       "N-Terminal ASH (COOH)"           231
biotype     520    CA      "N-Terminal ASH (COOH)"             8
biotype     521    C       "N-Terminal ASH (COOH)"             9
biotype     522    HN      "N-Terminal ASH (COOH)"           232
biotype     523    O       "N-Terminal ASH (COOH)"            11
biotype     524    HA      "N-Terminal ASH (COOH)"            12
biotype     525    N       "N-Terminal ASN"                  231
biotype     526    CA      "N-Terminal ASN"                    8
biotype     527    C       "N-Terminal ASN"                    9
biotype     528    HN      "N-Terminal ASN"                  232
biotype     529    O       "N-Terminal ASN"                   11
biotype     530    HA      "N-Terminal ASN"                   12
biotype     531    N       "N-Terminal GLU"                  231
biotype     532    CA      "N-Terminal GLU"                    8
biotype     533    C       "N-Terminal GLU"                    9
biotype     534    HN      "N-Terminal GLU"                  232
biotype     535    O       "N-Terminal GLU"                   11
biotype     536    HA      "N-Terminal GLU"                   12
biotype     537    N       "N-Terminal GLH (COOH)"           231
biotype     538    CA      "N-Terminal GLH (COOH)"             8
biotype     539    C       "N-Terminal GLH (COOH)"             9
biotype     540    HN      "N-Terminal GLH (COOH)"           232
biotype     541    O       "N-Terminal GLH (COOH)"            11
biotype     542    HA      "N-Terminal GLH (COOH)"            12
biotype     543    N       "N-Terminal GLN"                  231
biotype     544    CA      "N-Terminal GLN"                    8
biotype     545    C       "N-Terminal GLN"                    9
biotype     546    HN      "N-Terminal GLN"                  232
biotype     547    O       "N-Terminal GLN"                   11
biotype     548    HA      "N-Terminal GLN"                   12
biotype     549    N       "N-Terminal MET"                  231
biotype     550    CA      "N-Terminal MET"                    8
biotype     551    C       "N-Terminal MET"                    9
biotype     552    HN      "N-Terminal MET"                  232
biotype     553    O       "N-Terminal MET"                   11
biotype     554    HA      "N-Terminal MET"                   12
biotype     555    N       "N-Terminal LYS"                  231
biotype     556    CA      "N-Terminal LYS"                    8
biotype     557    C       "N-Terminal LYS"                    9
biotype     558    HN      "N-Terminal LYS"                  232
biotype     559    O       "N-Terminal LYS"                   11
biotype     560    HA      "N-Terminal LYS"                   12
biotype     561    N       "N-Terminal LYD (NH2)"            231
biotype     562    CA      "N-Terminal LYD (NH2)"              8
biotype     563    C       "N-Terminal LYD (NH2)"              9
biotype     564    HN      "N-Terminal LYD (NH2)"            232
biotype     565    O       "N-Terminal LYD (NH2)"             11
biotype     566    HA      "N-Terminal LYD (NH2)"             12
biotype     567    N       "N-Terminal ARG"                  231
biotype     568    CA      "N-Terminal ARG"                    8
biotype     569    C       "N-Terminal ARG"                    9
biotype     570    HN      "N-Terminal ARG"                  232
biotype     571    O       "N-Terminal ARG"                   11
biotype     572    HA      "N-Terminal ARG"                   12
biotype     573    N       "N-Terminal ORN"                  231
biotype     574    CA      "N-Terminal ORN"                    8
biotype     575    C       "N-Terminal ORN"                    9
biotype     576    HN      "N-Terminal ORN"                  232
biotype     577    O       "N-Terminal ORN"                   11
biotype     578    HA      "N-Terminal ORN"                   12
biotype     579    N       "N-Terminal AIB"                   -1
biotype     580    CA      "N-Terminal AIB"                   -1
biotype     581    C       "N-Terminal AIB"                   -1
biotype     582    HN      "N-Terminal AIB"                   -1
biotype     583    O       "N-Terminal AIB"                   -1
biotype     584    N       "C-Terminal GLY"                    1
biotype     585    CA      "C-Terminal GLY"                    2
biotype     586    C       "C-Terminal GLY"                  233
biotype     587    HN      "C-Terminal GLY"                    4
biotype     588    OXT     "C-Terminal GLY"                  234
biotype     589    HA      "C-Terminal GLY"                    6
biotype     590    N       "C-Terminal ALA"                    7
biotype     591    CA      "C-Terminal ALA"                    8
biotype     592    C       "C-Terminal ALA"                  233
biotype     593    HN      "C-Terminal ALA"                   10
biotype     594    OXT     "C-Terminal ALA"                  234
biotype     595    HA      "C-Terminal ALA"                   12
biotype     596    N       "C-Terminal VAL"                    7
biotype     597    CA      "C-Terminal VAL"                    8
biotype     598    C       "C-Terminal VAL"                  233
biotype     599    HN      "C-Terminal VAL"                   10
biotype     600    OXT     "C-Terminal VAL"                  234
biotype     601    HA      "C-Terminal VAL"                   12
biotype     602    N       "C-Terminal LEU"                    7
biotype     603    CA      "C-Terminal LEU"                    8
biotype     604    C       "C-Terminal LEU"                  233
biotype     605    HN      "C-Terminal LEU"                   10
biotype     606    OXT     "C-Terminal LEU"                  234
biotype     607    HA      "C-Terminal LEU"                   12
biotype     608    N       "C-Terminal ILE"                    7
biotype     609    CA      "C-Terminal ILE"                    8
biotype     610    C       "C-Terminal ILE"                  233
biotype     611    HN      "C-Terminal ILE"                   10
biotype     612    OXT     "C-Terminal ILE"                  234
biotype     613    HA      "C-Terminal ILE"                   12
biotype     614    N       "C-Terminal SER"                    7
biotype     615    CA      "C-Terminal SER"                    8
biotype     616    C       "C-Terminal SER"                  233
biotype     617    HN      "C-Terminal SER"                   10
biotype     618    OXT     "C-Terminal SER"                  234
biotype     619    HA      "C-Terminal SER"                   12
biotype     620    N       "C-Terminal THR"                    7
biotype     621    CA      "C-Terminal THR"                    8
biotype     622    C       "C-Terminal THR"                  233
biotype     623    HN      "C-Terminal THR"                   10
biotype     624    OXT     "C-Terminal THR"                  234
biotype     625    HA      "C-Terminal THR"                   12
biotype     626    N       "C-Terminal CYS (SH)"               7
biotype     627    CA      "C-Terminal CYS (SH)"               8
biotype     628    C       "C-Terminal CYS (SH)"             233
biotype     629    HN      "C-Terminal CYS (SH)"              10
biotype     630    OXT     "C-Terminal CYS (SH)"             234
biotype     631    HA      "C-Terminal CYS (SH)"              12
biotype     632    N       "C-Terminal CYX (SS)"               7
biotype     633    CA      "C-Terminal CYX (SS)"               8
biotype     634    C       "C-Terminal CYX (SS)"             233
biotype     635    HN      "C-Terminal CYX (SS)"              10
biotype     636    OXT     "C-Terminal CYX (SS)"             234
biotype     637    HA      "C-Terminal CYX (SS)"              12
biotype     638    N       "C-Terminal CYD (S-)"               7
biotype     639    CA      "C-Terminal CYD (S-)"              48
biotype     640    C       "C-Terminal CYD (S-)"             233
biotype     641    HN      "C-Terminal CYD (S-)"              10
biotype     642    OXT     "C-Terminal CYD (S-)"             234
biotype     643    HA      "C-Terminal CYD (S-)"              12
biotype     644    N       "C-Terminal PRO"                   50
biotype     645    CA      "C-Terminal PRO"                   51
biotype     646    C       "C-Terminal PRO"                  233
biotype     647    OXT     "C-Terminal PRO"                  234
biotype     648    HA      "C-Terminal PRO"                   54
biotype     649    N       "C-Terminal PHE"                    7
biotype     650    CA      "C-Terminal PHE"                    8
biotype     651    C       "C-Terminal PHE"                  233
biotype     652    HN      "C-Terminal PHE"                   10
biotype     653    OXT     "C-Terminal PHE"                  234
biotype     654    HA      "C-Terminal PHE"                   12
biotype     655    N       "C-Terminal TYR"                    7
biotype     656    CA      "C-Terminal TYR"                    8
biotype     657    C       "C-Terminal TYR"                  233
biotype     658    HN      "C-Terminal TYR"                   10
biotype     659    OXT     "C-Terminal TYR"                  234
biotype     660    HA      "C-Terminal TYR"                   12
biotype     661    N       "C-Terminal TYD (O-)"               7
biotype     662    CA      "C-Terminal TYD (O-)"               8
biotype     663    C       "C-Terminal TYD (O-)"             233
biotype     664    HN      "C-Terminal TYD (O-)"              10
biotype     665    OXT     "C-Terminal TYD (O-)"             234
biotype     666    HA      "C-Terminal TYD (O-)"              12
biotype     667    N       "C-Terminal TRP"                    7
biotype     668    CA      "C-Terminal TRP"                    8
biotype     669    C       "C-Terminal TRP"                  233
biotype     670    HN      "C-Terminal TRP"                   10
biotype     671    OXT     "C-Terminal TRP"                  234
biotype     672    HA      "C-Terminal TRP"                   12
biotype     673    N       "C-Terminal HIS (+)"                7
biotype     674    CA      "C-Terminal HIS (+)"                8
biotype     675    C       "C-Terminal HIS (+)"              233
biotype     676    HN      "C-Terminal HIS (+)"               10
biotype     677    OXT     "C-Terminal HIS (+)"              234
biotype     678    HA      "C-Terminal HIS (+)"               12
biotype     679    N       "C-Terminal HIS (HD)"               7
biotype     680    CA      "C-Terminal HIS (HD)"               8
biotype     681    C       "C-Terminal HIS (HD)"             233
biotype     682    HN      "C-Terminal HIS (HD)"              10
biotype     683    OXT     "C-Terminal HIS (HD)"             234
biotype     684    HA      "C-Terminal HIS (HD)"              12
biotype     685    N       "C-Terminal HIS (HE)"               7
biotype     686    CA      "C-Terminal HIS (HE)"               8
biotype     687    C       "C-Terminal HIS (HE)"             233
biotype     688    HN      "C-Terminal HIS (HE)"              10
biotype     689    OXT     "C-Terminal HIS (HE)"             234
biotype     690    HA      "C-Terminal HIS (HE)"              12
biotype     691    N       "C-Terminal ASP"                    7
biotype     692    CA      "C-Terminal ASP"                    8
biotype     693    C       "C-Terminal ASP"                  233
biotype     694    HN      "C-Terminal ASP"                   10
biotype     695    OXT     "C-Terminal ASP"                  234
biotype     696    HA      "C-Terminal ASP"                   12
biotype     697    N       "C-Terminal ASH (COOH)"             7
biotype     698    CA      "C-Terminal ASH (COOH)"             8
biotype     699    C       "C-Terminal ASH (COOH)"           233
biotype     700    HN      "C-Terminal ASH (COOH)"            10
biotype     701    OXT     "C-Terminal ASH (COOH)"           234
biotype     702    HA      "C-Terminal ASH (COOH)"            12
biotype     703    N       "C-Terminal ASN"                    7
biotype     704    CA      "C-Terminal ASN"                    8
biotype     705    C       "C-Terminal ASN"                  233
biotype     706    HN      "C-Terminal ASN"                   10
biotype     707    OXT     "C-Terminal ASN"                  234
biotype     708    HA      "C-Terminal ASN"                   12
biotype     709    N       "C-Terminal GLU"                    7
biotype     710    CA      "C-Terminal GLU"                    8
biotype     711    C       "C-Terminal GLU"                  233
biotype     712    HN      "C-Terminal GLU"                   10
biotype     713    OXT     "C-Terminal GLU"                  234
biotype     714    HA      "C-Terminal GLU"                   12
biotype     715    N       "C-Terminal GLH (COOH)"             7
biotype     716    CA      "C-Terminal GLH (COOH)"             8
biotype     717    C       "C-Terminal GLH (COOH)"           233
biotype     718    HN      "C-Terminal GLH (COOH)"            10
biotype     719    OXT     "C-Terminal GLH (COOH)"           234
biotype     720    HA      "C-Terminal GLH (COOH)"            12
biotype     721    N       "C-Terminal GLN"                    7
biotype     722    CA      "C-Terminal GLN"                    8
biotype     723    C       "C-Terminal GLN"                  233
biotype     724    HN      "C-Terminal GLN"                   10
biotype     725    OXT     "C-Terminal GLN"                  234
biotype     726    HA      "C-Terminal GLN"                   12
biotype     727    N       "C-Terminal MET"                    7
biotype     728    CA      "C-Terminal MET"                    8
biotype     729    C       "C-Terminal MET"                  233
biotype     730    HN      "C-Terminal MET"                   10
biotype     731    OXT     "C-Terminal MET"                  234
biotype     732    HA      "C-Terminal MET"                   12
biotype     733    N       "C-Terminal LYS"                    7
biotype     734    CA      "C-Terminal LYS"                    8
biotype     735    C       "C-Terminal LYS"                  233
biotype     736    HN      "C-Terminal LYS"                   10
biotype     737    OXT     "C-Terminal LYS"                  234
biotype     738    HA      "C-Terminal LYS"                   12
biotype     739    N       "C-Terminal LYD (NH2)"              7
biotype     740    CA      "C-Terminal LYD (NH2)"              8
biotype     741    C       "C-Terminal LYD (NH2)"            233
biotype     742    HN      "C-Terminal LYD (NH2)"             10
biotype     743    OXT     "C-Terminal LYD (NH2)"            234
biotype     744    HA      "C-Terminal LYD (NH2)"             12
biotype     745    N       "C-Terminal ARG"                    7
biotype     746    CA      "C-Terminal ARG"                    8
biotype     747    C       "C-Terminal ARG"                  233
biotype     748    HN      "C-Terminal ARG"                   10
biotype     749    OXT     "C-Terminal ARG"                  234
biotype     750    HA      "C-Terminal ARG"                   12
biotype     751    N       "C-Terminal ORN"                    7
biotype     752    CA      "C-Terminal ORN"                    8
biotype     753    C       "C-Terminal ORN"                  233
biotype     754    HN      "C-Terminal ORN"                   10
biotype     755    OXT     "C-Terminal ORN"                  234
biotype     756    HA      "C-Terminal ORN"                   12
biotype     757    N       "C-Terminal AIB"                   -1
biotype     758    CA      "C-Terminal AIB"                   -1
biotype     759    C       "C-Terminal AIB"                   -1
biotype     760    HN      "C-Terminal AIB"                   -1
biotype     761    OXT     "C-Terminal AIB"                   -1
biotype     762    N       "Deprotonated N-Terminus"          -1
biotype     763    H       "Deprotonated N-Terminus"          -1
biotype     764    C       "Formyl N-Terminus"                -1
biotype     765    H       "Formyl N-Terminus"                -1
biotype     766    O       "Formyl N-Terminus"                -1
biotype     767    CH3     "Acetyl N-Terminus"               221
biotype     768    H       "Acetyl N-Terminus"               222
biotype     769    C       "Acetyl N-Terminus"               223
biotype     770    O       "Acetyl N-Terminus"               224
biotype     771    C       "Protonated C-Terminus"           235
biotype     772    O       "Protonated C-Terminus"           236
biotype     773    OH      "Protonated C-Terminus"           237
biotype     774    HO      "Protonated C-Terminus"           238
biotype     775    N       "Amide C-Terminus"                225
biotype     776    HN      "Amide C-Terminus"                226
biotype     777    N       "N-MeAmide C-Terminus"            227
biotype     778    HN      "N-MeAmide C-Terminus"            228
biotype     779    CH3     "N-MeAmide C-Terminus"            229
biotype     780    H       "N-MeAmide C-Terminus"            230
biotype    2001    O       "Water"                           247
biotype    2002    H       "Water"                           248
biotype    2003    LI      "Lithium Ion"                     249
biotype    2004    NA      "Sodium Ion"                      250
biotype    2005    K       "Potassium Ion"                   251
biotype    2006    RB      "Rubidium Ion"                    252
biotype    2007    CS      "Cesium Ion"                      253
biotype    2008    MG      "Magnesium Ion"                   255
biotype    2009    CA      "Calcium Ion"                     256
biotype    2012    F       "Fluoride Ion"                    258
biotype    2013    CL      "Chloride Ion"                    259
biotype    2014    BR      "Bromide Ion"                     260
biotype    2015    I       "Iodide Ion"                      261
)**";
}
