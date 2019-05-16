const char* nacl_xyz1 =
R"**(     2  NaCl
     1  Na+    0.000000    0.000000    0.000000     7
     2  Cl-    2.200000    0.000000    0.000000    15
)**";

const char* nacl_xyz2 =
R"**(     2  NaCl
     1  Na+    0.000000    0.000000    0.000000     7
     2  Cl-    2.380000    0.000000    0.000000    15
)**";

const char* nacl_xyz3 =
R"**(     2  NaCl
     1  Na+    0.000000    0.000000    0.000000     7
     2  Cl-    2.500000    0.000000    0.000000    15
)**";

const char* nacl_key =
R"**(
parameters                 amoeba09
verbose

a-axis                       18.643
integrator                    RESPA
thermostat                    BUSSI
barostat                 MONTECARLO

neighbor-list
list-buffer                     0.1
vdw-cutoff                      2.6
mpole-cutoff                    5.0
)**";
