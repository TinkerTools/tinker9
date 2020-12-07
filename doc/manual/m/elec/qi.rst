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

Multipole Interaction in QI Frame
---------------------------------

==================================  ==================  =============
Gradients                           Global Frame        QI Frame
==================================  ==================  =============
:math:`\partial f(r)/\partial x_2`  :math:`f'(r)r_x/r`  0
:math:`\partial f(r)/\partial y_2`  :math:`f'(r)r_y/r`  0
:math:`\partial f(r)/\partial z_2`  :math:`f'(r)r_z/r`  :math:`f'(r)`
==================================  ==================  =============

where :math:`f'(r)=\partial f(r)/\partial r`.

================================  =================================================
Potential Terms                   Notes
================================  =================================================
:math:`\phi_1`                    :math:`\phi_1`
:math:`\phi'_{1x}`                :math:`\partial\phi_1/\partial x_1`
:math:`\phi'_{1y}`                :math:`\partial\phi_1/\partial y_1`
:math:`\phi'_{1z}`                :math:`\partial\phi_1/\partial z_1`
:math:`\phi''_{1xx}`              :math:`\partial^2\phi_1/\partial x_1^2`
:math:`\phi''_{1yy}`              :math:`\partial^2\phi_1/\partial y_1^2`
:math:`\phi''_{1zz}`              :math:`\partial^2\phi_1/\partial z_1^2`
:math:`\phi''_{1xy}`              :math:`\partial^2\phi_1/\partial x_1\partial y_1`
:math:`\phi''_{1xz}`              :math:`\partial^2\phi_1/\partial x_1\partial z_1`
:math:`\phi''_{1yz}`              :math:`\partial^2\phi_1/\partial y_1\partial z_1`
:math:`\phi_2`                    :math:`\phi_2`
:math:`\phi'_{2x}`                :math:`\partial\phi_2/\partial x_2`
:math:`\phi'_{2y}`                :math:`\partial\phi_2/\partial y_2`
:math:`\phi'_{2z}`                :math:`\partial\phi_2/\partial z_2`
:math:`\phi''_{2xx}`              :math:`\partial^2\phi_2/\partial x_2^2`
:math:`\phi''_{2yy}`              :math:`\partial^2\phi_2/\partial y_2^2`
:math:`\phi''_{2zz}`              :math:`\partial^2\phi_2/\partial z_2^2`
:math:`\phi''_{2xy}`              :math:`\partial^2\phi_2/\partial x_2\partial y_2`
:math:`\phi''_{2xz}`              :math:`\partial^2\phi_2/\partial x_2\partial z_2`
:math:`\phi''_{2yz}`              :math:`\partial^2\phi_2/\partial y_2\partial z_2`
================================  =================================================

Charge Terms
~~~~~~~~~~~~

:math:`A \leftarrow B` stands for :math:`A = A + B`.

.. math::

   \phi_1 &\leftarrow T_{12}^{(1,1)} C_2 = C_2 B_0, \\
   \phi'_1 &\leftarrow T_{12}^{(2:4,1)} C_2 = C_2 \begin{pmatrix}
   0 \\
   0 \\
   r B_1 \end{pmatrix}, \\
   \phi''_1 &\leftarrow T_{12}^{(5:13,1)} C_2 = \begin{pmatrix}
   0 \\
   0 \\
   0 \\
   0 \\
   0 \\
   0 \end{pmatrix}.

.. math::

   \phi_2 &\leftarrow T_{21}^{(1,1)} C_1 = C_1 B_0, \\
   \phi'_2 &\leftarrow T_{21}^{(2:4,1)} C_1 = C_1 \begin{pmatrix}
   0 \\
   0 \\
   -r B_1 \end{pmatrix}, \\
   \phi''_2 &\leftarrow T_{21}^{(5:13,1)} C_1 = \begin{pmatrix}
   0 \\
   0 \\
   0 \\
   0 \\
   0 \\
   0 \end{pmatrix}.

Dipole Terms
~~~~~~~~~~~~

Quadrupole Terms
~~~~~~~~~~~~~~~~
