Command Line GUI
================

Tinker programs can take interactive input. The prompt messages try to be
self-explanatory.
Following is an example to use the interactive interface for a short DHFR
simulation with an *AMOEBA* force field.

.. warning::

   The simulation will not run unless more keywords are added to the key file.
   This is only a demonstration to the interactive interface.

=======================  ===================  ===================
Item                     Value                Input
=======================  ===================  ===================
Coordinate File          dhfr2.xyz            dhfr2.xyz
Simulations Steps        1000                 1000
Time Step                2 fs                 2
Time Between Saves       1 ps                 1
Ensemble                 NVT                  2
Thermostat Temperature   298 K                298
=======================  ===================  ===================

.. code-block:: text

   zw@Blade:~/tutorial/example$ tinker9 dynamic

        ###############################################################
      ###################################################################
     ###                                                               ###
    ###        Tinker9  ---  Software Tools for Molecular Design        ###
    ##                                                                   ##
    ##                    Version 1.0.0-rc   Jan 2021                    ##
    ###                       All Rights Reserved                       ###
     ###                                                               ###
      ###################################################################
        ###############################################################
    Compiled at:  20:54:14  Feb  3 2021
    Commit Date:  Wed Feb 3 20:52:06 2021 -0600
    Commit:       3b0791c0


    Enter Cartesian Coordinate File Name :  dhfr2.xyz

    ################################################################
    Joint Amber-CHARMM Benchmark on Dihydrofolate Reductase in Water
    23558 Atoms, 62.23 Ang Cube, 9 Ang Nonbond Cutoffs, 64x64x64 PME
    ################################################################

    Enter the Number of Dynamics Steps to be Taken :  1000

    Enter the Time Step Length in Femtoseconds [1.0] :  2

    Enter Time between Saves in Picoseconds [0.1] :  1

    Available Statistical Mechanical Ensembles :
       (1) Microcanonical (NVE)
       (2) Canonical (NVT)
       (3) Isoenthalpic-Isobaric (NPH)
       (4) Isothermal-Isobaric (NPT)
    Enter the Number of the Desired Choice  [1] :  2

    Enter the Desired Temperature in Degrees K [298] :  298

Once you are familar with the interactive interface, you can simply append the
interactive input to the program name.

.. code-block:: text

   zw@Blade:~/tutorial/example$ tinker9 dynamic dhfr2.xyz 1000 2 1 2 298

        ###############################################################
      ###################################################################
     ###                                                               ###
    ###        Tinker9  ---  Software Tools for Molecular Design        ###
    ##                                                                   ##
    ##                    Version 1.0.0-rc   Jan 2021                    ##
    ###                       All Rights Reserved                       ###
     ###                                                               ###
      ###################################################################
        ###############################################################
    Compiled at:  20:54:14  Feb  3 2021
    Commit Date:  Wed Feb 3 20:52:06 2021 -0600
    Commit:       3b0791c0


    ################################################################
    Joint Amber-CHARMM Benchmark on Dihydrofolate Reductase in Water
    23558 Atoms, 62.23 Ang Cube, 9 Ang Nonbond Cutoffs, 64x64x64 PME
    ################################################################

You can also use *stdin redirection* to run Tinker programs. Save the
command line arguments in a file and use it as follows.

.. code-block:: text

   zw@Blade:~/tutorial/example$ cat args.txt 
   dhfr2.xyz
   1000
   2
   1
   2
   298
   zw@Blade:~/tutorial/example$ tinker9 dynamic < args.txt
