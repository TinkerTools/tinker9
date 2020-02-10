# Permanent Multipole (Up to Quadrupole) and Polarization Interactions
This file documented the pairwise quardupole interactions between atoms `i` and
`k`. Sometimes atoms `i` and `k` were referred to as atoms 1 and 2,
respectively.


## Units and Multipole Moments
The electrostatic potential at $\bm{r}$ due to the charge distribution nearby is

$$
\varphi(\bm{r}) = \frac{1}{4\pi\epsilon_0}
                  \int d\bm{s} \frac{\rho(\bm{s})}{|\bm{r}-\bm{s}|}.
$$

Tinker uses a variable `electric` (of the `chgpot` module) to represent the
factor $(4\pi\epsilon_0)^{-1}$.
Its unit is (Å/(e·e))(kcal/mol), and its current default
magnitude is 332.063713, which is a constant defined by `coulomb` of the `units`
module. The default value is editable by the `ELECTRIC` keyword.

Expanding $|\bm{r}-\bm{s}|^{-1}$ by Taylor series, the potential can be
rewritten as

$$
\frac{1}{4\pi\epsilon_0}
\left(
             \left[             \int d\bm{s} \rho(\bm{s})        \right]
             \frac{1}{r}
- \sum_i     \left[             \int d\bm{s} \rho(\bm{s})s_i     \right]
             \nabla_i \frac{1}{r}
+ \sum_{i,j} \left[ \frac{1}{2} \int d\bm{s} \rho(\bm{s})s_i s_j \right]
             \nabla_i \nabla_j \frac{1}{r}
- \cdots
\right),
$$
where three pairs of square brackets define monopole ($C$, charge), dipole
($D_i$), and quadrupole ($Q_{ij}^*$) moments, respectively.


Tinker parameters use slightly different definitions though.
The parameters directly came from quantum packages, where for historical
reasons, *traceless quadrupoles* and atomic units were widely used,
so Bohr must be converted to Å. Moreover, Tinker will scale the
*traceless quadrupoles* by 1/3 before energy evaluation.
The reason is as follows.

Potential due to the quadrupole is

$$
\frac{1}{4\pi\epsilon_0} \sum_{i,j} \left[
\frac{1}{2} \int d\bm{s} \rho(\bm{s})s_i s_j\right]
\frac{3r_i r_j - r^2\delta_{ij}}{r^5},
$$
which can be rewritten as

$$
\frac{1}{4\pi\epsilon_0} \sum_{i,j} \left[
\frac{1}{2} \int d\bm{s} \rho(\bm{s})(3s_i s_j - s^2\delta_{ij})\right]
\frac{r_i r_j}{r^5}.
$$

So the traceless quadrupole tensor can be defined based on the equation above as

$$
\Theta_{ij} = \frac{1}{2} \int d\bm{s} \rho(\bm{s})(3s_i s_j - s^2\delta_{ij}). 
$$
We can easily confirm that $\sum_i^{x,y,z}(3s_i s_i-s^2) = 0$, therefore

$$
\Theta_{ij} = 3Q_{ij}^* - \delta_{ij}\sum_k^{x,y,z}Q_{kk}^*.
$$
Tinker defines and uses $Q_{ij} = \Theta_{ij}/3$ in the energy evaluation, thus
the energy expression is the same as the non-traceless quadrupoles, but still
has the advantage to optimize the expressions because the quadrupoles are traceless. 


## Derivatives of 1/r
As we can see in the previous section, high-order gradients of $1/r$ are
necessary for the electrostatic potential.
The displacement $\bm{r}_{12} = \bm{r}_2 - \bm{r}_1 = (r_x,r_y,r_z)$ starts from
atom 1 and ends in atom 2.
For any scalar $f$, it is easy to confirm that

$$
\frac{\partial f}{\partial \bm{r}} =  \frac{\partial f}{\partial \bm{r}_2}
                                   = -\frac{\partial f}{\partial \bm{r}_1},
$$
and for $f(r,r_x,r_y,r_z)$ where $r$ is a function of $r_x,r_y,r_z$,

$$
\nabla f(r,r_x,r_y,r_z) = \frac{\partial f(r)}{\partial r}\nabla r +
\sum_j\frac{\partial f(r_j)}{\partial r_j}\nabla r_j =
\frac{\partial f(r)}{\partial r}\frac{\bm{r}}{r} +
\nabla f(r_x,r_y,r_z).
$$

So the derivatives of $1/r$ with respect to $a$, $b$, $c$, etc., which are
one of the x, y, z directions, are

| Terms | Expressions |
|-------|-------------|
| $T_0$ | $\lambda_1/r$ |
| $T_1$ | $-\lambda_3 r_a/r^3$ |
| $T_2$ | $\lambda_5 3r_a r_b/r^5 -\lambda_3\delta_{ab}/r^3$ |
| $T_3$ | $-\lambda_7 15 r_a r_b r_c/r^7 +\lambda_5 3\Sigma r_a\delta_{bc}/r^5$ |
| $T_4$ | $\lambda_9 105 r_a r_b r_c r_d/r^9 -\lambda_7 15\Sigma r_a r_b\delta_{cd}/r^7 +\lambda_5 3\Sigma\delta_{ab}\delta_{cd}/r^5$ |
| $T_5$ | $-\lambda_{11}945 r_a r_b r_c r_d r_e/r^{11} +\lambda_9 105\Sigma r_a r_b r_c\delta_{de}/r^9 -\lambda_7 15\Sigma r_a\delta_{bc}\delta_{de}/r^7$ |

where

| Dampings   | Values of $\lambda$s   |
|------------|------------------------|
| non-EWALD  | 1                      |
| EWALD      | `erfc` damping factors |
| Thole      | [J.Phys.Chem.B, 107, 5933 (2003)](https://dx.doi.org/10.1021/jp027815+) |

In the following sections, $\lambda_k (k-2)!! / r^k$ will be denoted by $R_k$.


## Rotation Matrix
Rotation matrix $R$ maps a vector (1D tensor) $\bm{x}$ from frame a to b,
$R\bm{x}_a = \bm{x}_b$. To preserve the dot product of any two vectors $\bm{x}$
and $\bm{y}$ in either frame, $\bm{x}_b^t \bm{y}_b = \bm{x}_a^t \bm{y}_a$, the
rotation matrix must be orthogonal, i.e., $R^t R = I$.

Rotating a 2D tensor (matrix) $Q_a$ from frame a to b, we have $Q_b$, and the
quadratic form must be preserved in both frames.
From $\bm{x}_b^t Q_b \bm{y}_b = \bm{x}_a^t Q_a \bm{y}_a$, we get

$$
Q_b = R Q_a R^t.
$$

Rotation matrix also preserves the trace of matrix $Q$.


## Pairwise Energy, Torque, and Force
For each atom, $M$ denotes its charge ($C$), dipole ($D$), traceless quadrupole
($Q$), and induced dipole ($\mu$) in the *global* or *quasi-internal (QI)*
frame, depending on the context. The electrostatic energy of atom 1 due to atom
2 is given by

$$
U_1 = \varphi_1(2) C_1 + \varphi_1'(2) D_1 + \varphi_1''(2) Q_1 + \cdots,
$$
and must be equal to $U_2$. Therefore, the pairwise energy of atoms 1 and 2 can
alternatively be defined as

$$
U = M_2 T(\bm{r}) M_1,
$$
where the $T$ matrix is only a function of $\bm{r}$.


## About This File
The pdf version was generated by
```bash
pandoc mpole.md -o mpole.pdf \
-Vfontfamily=fouriernc \
-Vheader-includes='\hypersetup{colorlinks=true}' \
-Vheader-includes='\usepackage{bm}'
```


Zhi Wang


Feb 9, 2020.
