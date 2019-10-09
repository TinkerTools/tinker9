# Procedures

## Main
   - Initialize the Fortran runtime library.
   - Initialize the canonical Tinker library.
      - `call initial`, `call getxyz`, `call mechanic`, etc.
   - Initialize `tinker.gpu`.
      - Set the global flag (`rc_flag`) for the resource management.
      - `initialize();` (see section **Initialize**)
   - Stop `tinker.gpu`.
      - `finish();` (see section **Finish**)
   - Stop the canonical Tinker library.
      - `call final`
   - Stop the Fortran runtime library.


## Initialize


## Finish

