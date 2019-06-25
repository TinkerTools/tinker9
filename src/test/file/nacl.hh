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

a-axis                       18.643
integrator                    RESPA
thermostat                    BUSSI
barostat                 MONTECARLO

neighbor-list
list-buffer                     0.1
vdw-cutoff                      2.6
mpole-cutoff                    5.0
)**";

const char* nacl_xyz4 =
R"**(     2  NaCl
     1  Na+    0.123000    0.134000    0.145000     7     2
     2  Cl-    1.123000    1.234000    1.347000    15     1
)**";

const char* nacl_key4 =
R"**(
parameters                 amoeba09

integrator                    RESPA
thermostat                    BUSSI
barostat                 MONTECARLO

neighbor-list
list-buffer                     0.1
vdw-cutoff                      2.6
mpole-cutoff                    5.0

a-axis                       20.000
pme-grid                   40 40 40
ewald
ewald-cutoff                    7.0

# originally N2 parameters
multipole     7   15    0               1.00000
                                        0.00000    0.00000    0.12578
                                        0.16329
                                        0.00000    0.16329
                                        0.00000    0.00000   -0.32658
multipole    15    7    0              -1.00000
                                        0.00000    0.00000    0.12578
                                        0.16329
                                        0.00000    0.16329
                                        0.00000    0.00000   -0.32658
)**";
