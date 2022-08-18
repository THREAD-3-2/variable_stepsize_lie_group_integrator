.. _vi:

========================
 Variational integrator
========================

.. _continuous_vi:


Continuous setting
-------------------

In this work, the electromechanically coupled beam dynamics is approximated within the constrained discrete variational scheme 
with the null space projection. The Lagrange-d'Alembert principle can be extended to constrained systems by enforcing the constraints via Lagrange multipliers as

.. math::
    
    \begin{align}
        \delta \int_{0}^{T} \left[  L( \mathbf{q}, \dot{\mathbf{q}}) - \mathbf{g}^T(\mathbf{q})\cdot \boldsymbol\lambda \right] dt + \int_{0}^{T}\mathbf{f}^{\rm ext}( t) \cdot \delta \mathbf{q}dt=0,
    \end{align}

where :math:`\mathbf{q}` is the configuration, :math:`L( \mathbf{q}, \dot{\mathbf{q}})` is the Lagrangian, :math:`\mathbf{g}` represents holonomic constraints, 
:math:`\boldsymbol\lambda` is the Lagrangian multiplier and :math:`\mathbf{f}^{\rm ext}(t)` is the external force. By considering 
the electrical effect in geometrically exact beam, the electric potential :math:`\phi_o` and the incremental variables :math:`(\alpha, \beta )` are 
treated as the electrical degrees of freedom :math:`\boldsymbol\phi=\begin{bmatrix}  \phi_o& \alpha& \beta \end{bmatrix}` such that 
the configuration of the beam model is extended to

.. math::
    
    \begin{align}
        \mathbf{q}=\begin{bmatrix} 
                    \boldsymbol{\varphi}& \mathbf{d}_1
                    &\mathbf{d}_2&\mathbf{d}_3& \boldsymbol\phi
                \end{bmatrix}^T.
    \end{align}

According to the kinematic assumptions in geometrically exact beams, the directors have to fulfill the orthogonal constraints

.. math::
    
    \begin{align}
        \mathbf{g}( \mathbf{q})=\begin{bmatrix}                           
                            \frac{1}{2}(\mathbf{d}_1^T \mathbf{d}_1-1) \\                                               
                            \frac{1}{2}(\mathbf{d}_2^T \mathbf{d}_2-1)\\   
                            \frac{1}{2}(\mathbf{d}_3^T \mathbf{d}_3-1) \\  
                            \mathbf{d}_1^T\mathbf{d}_2\\
                            \mathbf{d}_1^T\mathbf{d}_3\\ 
                            \mathbf{d}_2^T\mathbf{d}_3
                            \end{bmatrix}=\mathbf{0}.
    \end{align}

The continuous Lagrangian contains the difference between the kinetic energy :math:`T(\dot{\mathbf{q}})` and the internal 
potential energy :math:`V(\mathbf{q})`

.. math::
    
    \begin{align}
        L(\mathbf{q}, \dot{\mathbf{q}})= T(\dot{\mathbf{q}}) -V(\mathbf{q}).
    \end{align}

Since the electrical variables do not contribute to the kinetic energy, the kinetic energy for geometrically exact beams is computed as

.. math::

    \begin{align}
        T=  \int_c \left( \frac{1}{2} A_{\rho} \left| \dot{\boldsymbol  \varphi} \right| ^2 + \frac{1}{2} \sum_{i=1}^{2} M^i_{\rho}\left| \dot{\mathbf{d}}_i \right| ^2 \right) ds, \label{T}
    \end{align}

where :math:`A_{\rho}` is the mass density per reference arc-length and :math:`M^i_{\rho}` are the principle mass moments of 
inertia of cross section.   In accordance with the configuration for beam, the component of the consistent mass matrix corresponding 
to the electrical degree of freedom :math:`\boldsymbol\phi` will be zero.

For the coupled hyperelastic material in DEA, the internal potential energy is computed  by an integration of the beam strain 
energy density :math:`\Omega_b` over the beam center line

.. math::
    
    \begin{align}
        V(\mathbf{q}) = \int_c \Omega_b (s) ds.
    \end{align}

The external force :math:`\mathbf{f}^{\rm ext}` contains all non-conservative forces, such as the viscoelastic effect in this work. 
Based on the Kelvin-Voigt model, the non-conservative work contributed by the viscoelastic effect is given by

.. math::
    
    \begin{align}
        W^{\rm vis}=\int_{B_0}  \mathbf{P}^{\rm vis}: \mathbf{F}dV,
    \end{align}

where the work is computed from the two conjugate quantities being the first Piola-Kirchhoff stress :math:`\mathbf{P}^{\rm vis}` from 
the Kelvin-Voigt model and the deformation gradient :math:`\mathbf{F}`. In this case, the external force corresponding to the viscoelastic 
effect can be formulated as

.. math::
    
    \begin{align}
        \mathbf{f}^{\rm vis}(\mathbf{q},\dot{\mathbf{q}})=\frac{\partial W^{\rm vis}}{\partial \mathbf{q}}= \int_{B_0} \frac{\partial W^{\rm vis}}{\partial \mathbf{F} } : \frac{\partial \mathbf{F}}{\partial \mathbf{q} } dV =  \int_c \int_{\Sigma} \mathbf{P}^{\rm vis}: \frac{\partial \mathbf{F}}{\partial \mathbf{q} }dAds.
    \end{align}

.. _discrete_vi:


Discrete Euler-Lagrange equations
---------------------------------

The beam is first spatially discretized with the 1D finite elements, where one-dimensional Lagrange-type linear shape functions are 
applied in the discretization of beam configuration :math:`\mathbf{q}`. In this case, the beam directors are directly discretized in 
space together with the beam centroids. Then the variational integration scheme is applied to temporally discretize the action of 
the dynamic system, by which the good long term energy behavior can be obtained. In the variational integration scheme, the action 
integral within the time interval :math:`(t_n,t_{n+1})` is approximated with the discrete Lagrangian :math:`L_d` as

.. math::
    
    \begin{align}
        \int_{t_n}^{t_{n+1}} L(\mathbf{q},  \dot{\mathbf{q}})dt \approx L_d(\mathbf{q}_n,\mathbf{q}_{n+1}) = \Delta t L(\frac{\mathbf{q}_{n+1}+\mathbf{q}_n}{2},\frac{\mathbf{q}_{n+1}-\mathbf{q}_n}{\Delta t}),
    \end{align}

where the discrete Lagrangian :math:`L_d` is computed by applying the finite difference approximation to the velocity :math:`\dot{\mathbf{q}}` 
and the midpoint rule to the configuration :math:`\mathbf{q}`, i.e.

.. math::
    
    \begin{align}
        \dot{\mathbf{q}}\approx \frac{\mathbf{q}_{n+1}-\mathbf{q}_n}{\Delta t}, \;\;\;\;\; \mathbf{q}\approx\frac{\mathbf{q}_{n+1}+\mathbf{q}_n}{2} \label{mid}.
    \end{align}

After the temporal discretization, the discrete Euler-Lagrange equations can be obtained by taking the variation of the discrete 
action and requiring stationarity. To eliminate the constraint forces :math:`\boldsymbol\lambda` from the system, the nodal 
reparametrization :math:`\mathbf{q}_{n+1} = \mathbf{F}_d (\mathbf{u}_{n+1}, \mathbf{q}_{n})` and the discrete null space 
matrix :math:`\mathbf{P}_d` are applied to the discrete Euler-Lagrange equations leading to

.. math::
    
    \begin{align}
        \mathbf{P}_d^T(\mathbf{q}_n) \left[ \frac{\partial L_d(\mathbf{q}_{n-1}, \mathbf{q}_{n})}{\partial \mathbf{q}_{n}} + \frac{\partial L_d\left( \mathbf{q}_{n}, \mathbf{F}_d (\mathbf{u}_{n+1}, \mathbf{q}_{n})\right) }{\partial \mathbf{q}_{n}} + \mathbf{f}_n^{\rm ext-} + \mathbf{f}_{n-1}^{\rm ext+} \right] = \mathbf{0},
    \end{align}

where :math:`\mathbf{u}_{n+1}` is the generalized configuration acting as the unknown variable, :math:`\mathbf{f}_n^{\rm ext-}` and 
:math:`\mathbf{f}_{n-1}^{\rm ext+}` are the discrete generalized external forces evaluated as

.. math::
    
    \begin{align}
        \mathbf{f}_n^{\rm ext-}=\frac{\Delta t}{2} \mathbf{f}^{\rm vis} (\frac{\mathbf{q}_{n+1}+\mathbf{q}_n}{2},\frac{\mathbf{q}_{n+1}-\mathbf{q}_n}{\Delta t}),
        \mathbf{f}_{n-1}^{\rm ext+}=\frac{\Delta t}{2} \mathbf{f}^{\rm vis} (\frac{\mathbf{q}_{n-1}+\mathbf{q}_n}{2},\frac{\mathbf{q}_{n-1}-\mathbf{q}_n}{\Delta t}).
    \end{align}

The internal null space matrix at time :math:`t_n` is written as

.. math::
    
    \begin{align}
        \mathbf{P}_{\rm int}(\mathbf{q}_n)=
        \begin{bmatrix}                           
            \mathbf{I} & \mathbf{0} & \mathbf{0}\\                                               
            \mathbf{0} & -\hat{\mathbf{d}}_{1,n} & \mathbf{0}\\   
            \mathbf{0} & -\hat{\mathbf{d}}_{2,n} & \mathbf{0}\\  
            \mathbf{0} & -\hat{\mathbf{d}}_{3,n} & \mathbf{0}\\
            \mathbf{0}& \mathbf{0}& \mathbf{I}
        \end{bmatrix}, \label{null}
    \end{align}

where :math:`\hat{\mathbf{d}}_{i,n}` denotes the skew-symmetric matrix corresponding to the director vector 
:math:`\mathbf{d}_{i,n}` at :math:`t_n` and :math:`\mathbf{I}` is the 3 by 3 identity matrix. For a multibody dynamic system 
composed of flexible beam actuators, rigid bodies, joints and constraints, the null space matrix can be designed by considering 
the electric potential as extra degree of freedom as well.

To solve the discrete Euler-Lagrange equations efficiently, the system can be reduced further into the minimal possible dimensions 
by use of the nodal reparametrization. The generalized configuration of the electromechanically coupled beam is specified by

.. math::
    
    \begin{align}
        \mathbf{u}=\begin{bmatrix} \mathbf{u}_\varphi& \boldsymbol\theta&\mathbf{v}\end{bmatrix}^T
    \end{align}

with :math:`\mathbf{u}_\varphi,  \boldsymbol\theta$ and $\mathbf{v}` characterizing the incremental displacement, the incremental 
rotation and the incremental electric potential, respectively. In this case, the nodal configuration for the next time step can 
be updated as

.. math::
    
    \begin{align}
        \mathbf{q}_{n+1} = \mathbf{F}_d ({\mathbf{u}_{n+1}}, \mathbf{q}_{n})=
        \begin{bmatrix}
            \boldsymbol{\varphi}_n+\mathbf{u}_\varphi\\
            {\rm exp}(\hat{\boldsymbol\theta})\cdot \mathbf{d}_{1,n}\\
            {\rm exp}(\hat{\boldsymbol\theta})\cdot \mathbf{d}_{2,n}\\
            {\rm exp}(\hat{\boldsymbol\theta})\cdot \mathbf{d}_{3,n}\\
            \boldsymbol{\phi} + \mathbf{v}
        \end{bmatrix}.
    \end{align}

By means of the nodal reparametrization, the unknowns of the discrete Euler-Lagrange equation are changed from :math:`\mathbf{q}_{n+1}` to the generalized 
variables :math:`\mathbf{u}_{n+1}`. The nonlinear equation system is solved by use of the  Newton-Rapson algorithm with the tangent 
matrix at iteration :math:`i`

.. math::
    
    \begin{align}
        \mathbf{K}_T^i
        =\mathbf{P}^T(\mathbf{q}_n)  \frac{\partial \mathbf{R}^L(\mathbf{q}_{n+1}^i)}{\partial \mathbf{q}_{n+1}^i} \frac{\partial \mathbf{q}_{n+1}(\mathbf{u}_{n+1}^i)}{\partial \mathbf{u}_{n+1}^i},
    \end{align}

in which :math:`\mathbf{R}^L(\mathbf{q}_{n+1})` is the residual of the discrete Euler-Lagrange equation.