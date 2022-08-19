.. _RKMK_var_step:

====================================
RKMK methods with variable step size
====================================

The underlying idea of RKMK methods is to express a vector field :math:`F\in\mathfrak{X}(\mathcal{M})` as :math:`F\vert_m = \psi_*(f(m))\vert_m` , where :math:`\psi_*` is the infinitesimal generator of :math:`\psi`, a transitive action on :math:`\mathcal{M}`, and :math:`f:\M\rightarrow\mathfrak{g}`. This allows us to transform the problem from the manifold :math:`\mathcal{M}` to the Lie algebra :math:`\g`, on which we can perform a time step integration. We then map the result back to :math:`\mathcal{M}`, and repeat this up to the final integration time.  More explicitly, let $h_n$ be the size of the $n-$th time step, we then update :math:`y_n\in\mathcal{M}` to :math:`y_{n+1}` by

.. math::
    :name: eq:1
    
    \begin{align}
        \begin{cases}
        \sigma(0) = 0\in\mathfrak{g},\\
        \dot{\sigma}(t) = \dexp_{\sigma(t)}^{-1}\circ f\circ \psi (\exp(\sigma(t)),y_n)\in T_{\sigma(t)}\mathfrak{g}, \\
        y_{n+1} = \psi(\exp(\sigma_1),y_n)\in \mathcal{M},
        \end{cases}
    \end{align}

where :math:`\sigma_1\approx \sigma(h_n)\in\mathfrak{g}` is computed with a Runge-Kutta method.


One approach for varying the step size is based on embedded Runge--Kutta pairs for vector spaces. This approach consists of a principal method of order :math:`p`, used to propagate the numerical solution, together with some auxiliary method, of order :math:`\tilde{p}<p`, that is only used to obtain an estimate of the local error. This local error estimate is in turn used to derive a step size adjustment formula that attempts to keep the local error estimate approximately equal to some user-defined tolerance :math:`\tol` in every step.
Both methods are applied to solve the ODE for :math:`\sigma(t)` in \eqref{286eq:update}, yielding two approximations :math:`\sigma_1` and :math:`\tilde{\sigma}_1` respectively, using the same step size :math:`h_n`. Now, some distance measure between :math:`\sigma_1` and  :math:`\tilde{\sigma}_1` provides an estimate :math:`e_{n+1}$ for the size of the local truncation error. Thus,
:math:`e_{n+1}=C
h_{n+1}^{\aux{p}+1}+\mathcal{O}(h^{\aux{p}+2})`. Aiming at :math:`e_{n+1}\approx\tol` in every step, one may use a formula of the type

.. math::
    :name: eq:2
    \begin{align}
        h_{n+1} = \theta\left(\frac{\tol}{e_{n+1}}\right)^{\tfrac{1}{\aux{p}+1}}\, h_n
    \end{align}

where :math:`\theta` is typically chosen between :math:`0.8` and :math:`0.9`.
If :math:`e_n>\tol`, the step is rejected. Hence, we can redo the step with the step size obtained by the same formula.