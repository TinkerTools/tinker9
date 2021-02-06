Tutorials
=========

Let's create a directory *tutorial* under the home directory and copy the
directories *tinker9/example* and *tinker9/params* here.

.. code-block:: text

   zw@Blade:~$ mkdir ~/tutorial; cd ~/tutorial
   zw@Blade:~/tutorial$ cp -r ~/tinker9/{example,params} .
   zw@Blade:~/tutorial$ ls
   example/  params/
   zw@Blade:~/tutorial$ cd example; ls
   ar94.key  dhfr2.key  dhfr.key  dhfr.seq
   ar94.xyz  dhfr2.xyz  dhfr.pdb  dhfr.xyz

This tutorial also assumes the executable *tinker9* is in your *PATH*. If not,
prefix the directory to the *tinker9* executable.

.. toctree::

   filetypes
   cmdgui
   analyze
   minimize
   dynamic
