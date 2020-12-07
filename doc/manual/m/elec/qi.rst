Quasi-Internal Frame
====================

Rotation Matrix
---------------

Consider two vectors *u*, *v* and two reference frames *A*, *B*.
*R* is the rotation matrix of the *axes* such that

.. math::

   R u_A = u_B,

   R v_A = v_B.

Since :math:`u_A^T v_A=u_B^T v_B`,

.. math::

   R^T R=I.

A 2-D tensor, e.g., quadrupole moment *Q*, in two reference frames are
associated by

.. math::

   u_A^T Q_A v_A = u_B^T Q_B v_B.

It is easy to prove that

.. math::

   R Q_A R^T = Q_B.

Two common transformations used in Tinker are:

- From (A) *Local Frame* (in which parameters are provided)
  to (B) *Global Frame* (in which the calculation is done);
- From (A) *Global Frame* (for direct pairwise electrostatics)
  to (B) *Quasi-Internal (QI) Frame* (for optimized algebra),
  as shown in :numref:`fig-qiframe`.

.. _fig-qiframe:
.. figure:: ../fig/qiframe.*
   :width: 300 px
   :align: center

   Global frame *g* and QI frame *i* of atoms 1 and 2.
   The z direction of this QI frame is chosen along the distance vector.

Optimized Algebra in QI Frame
-----------------------------
