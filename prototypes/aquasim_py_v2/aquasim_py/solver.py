"""Replacement for dassl.cpp/das.cpp: implicit DAE  G(t, y, y') = 0 with banded Jacobian.
Preferred: scikits.odes 'ida' (SUNDIALS IDA = DASSL successor, same BDF/variable order).
Fallback: index-1 reduction -> scipy.integrate.solve_ivp(method='BDF', jac_sparsity=...).
"""
def integrate(residual, t_span, y0, yp0, algebraic_mask, t_out, rtol=1e-6, atol=1e-8):
    try:
        from scikits.odes import dae
    except ImportError as e:
        raise ImportError("pip install scikits.odes (needs SUNDIALS)") from e
    def res(t, y, yp, out):
        out[:] = residual(t, y, yp)
    solver = dae("ida", res, rtol=rtol, atol=atol,
                 algebraic_vars_idx=[i for i, a in enumerate(algebraic_mask) if a],
                 compute_initcond="yp0", old_api=False)
    return solver.solve(t_out, y0, yp0)
