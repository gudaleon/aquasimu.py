"""DREAM(ZS) sampler -- own numpy implementation.

Algorithm after Vrugt (2016), "Markov chain Monte Carlo simulation using the DREAM software
package: Theory, concepts, and MATLAB implementation", Environ. Model. Softw. 75, 273-316, and
ter Braak & Vrugt (2008) for the DE-MC_ZS update:

  * N parallel chains; proposals are built from differences of states drawn from an ARCHIVE Z of
    past states (thinned every K generations), so few chains suffice.
  * jump:   x_p = x + (1 + lambda) * gamma(delta, d') * [ sum_{j=1}^{delta} Z_r1(j) - sum Z_r2(j) ] + eps
            lambda ~ U(-c, c),  eps ~ N(0, c*^2),  delta ~ U{1..delta_max}
            gamma = 2.38 / sqrt(2 delta d')  (d' = number of dimensions updated), gamma = 1 with
            probability 0.2 (mode jumps)
  * randomised subspace sampling: each dimension is updated with crossover probability CR drawn
    from nCR values; the CR distribution is adapted during burn-in from the normalised squared
    jumping distance.
  * Metropolis acceptance; reflection at the prior bounds.
  * Gelman-Rubin R-hat on the second half of the chains, monitored during the run.

NOT implemented (documented simplification): the snooker update and the outlier-chain
replacement of the full DREAM(ZS) package.  Priors are uniform boxes in the sampling space.
"""
import numpy as np
import multiprocessing as mp


def gelman_rubin(chains):
    """chains: (n_samples, n_chains, d) -> R-hat per parameter (uses all supplied samples)."""
    n, m, d = chains.shape
    means = chains.mean(axis=0)                       # (m, d)
    W = chains.var(axis=0, ddof=1).mean(axis=0)       # within
    B = n * means.var(axis=0, ddof=1)                 # between
    var_hat = (n - 1) / n * W + B / n
    return np.sqrt(var_hat / np.maximum(W, 1e-300))


def _eval(args):
    f, x = args
    return f(x)


def dream_zs(logpost, lo, hi, n_chains=8, n_gen=6000, m0_factor=10, thin_K=10, n_cr=3, delta_max=3,
             p_gamma1=0.2, c=0.1, cstar=1e-6, adapt_frac=0.1, seed=1, n_proc=8, rhat_every=100,
             verbose=True, x0=None):
    """Run DREAM(ZS).  logpost(x) -> log posterior (use -inf outside support / for failures).

    Returns dict(chains (n_gen+1, n_chains, d), logp (n_gen+1, n_chains), accept_rate, rhat_trace,
    rhat_gens, archive, cr_probs)."""
    rng = np.random.default_rng(seed)
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    d = len(lo)
    span = hi - lo
    # ---- archive Z initialised from the prior ------------------------------------------------------
    m0 = m0_factor * d
    Z = np.empty((m0 + (n_gen // thin_K + 2) * n_chains, d))
    Z[:m0] = lo + span * rng.random((m0, d))
    nz = m0
    # ---- chains initialised from the prior ----------------------------------------------------------
    X = lo + span * rng.random((n_chains, d)) if x0 is None else np.array(x0, float)
    pool = mp.get_context("fork").Pool(n_proc) if n_proc > 1 else None
    pmap = (lambda xs: pool.map(_eval, [(logpost, x) for x in xs])) if pool else (lambda xs: [logpost(x) for x in xs])
    try:
        lp = np.array(pmap(list(X)))
        chains = np.empty((n_gen + 1, n_chains, d)); logp = np.empty((n_gen + 1, n_chains))
        chains[0], logp[0] = X, lp
        cr_vals = np.arange(1, n_cr + 1) / n_cr
        pcr = np.full(n_cr, 1.0 / n_cr)
        sq_jump = np.zeros(n_cr); cr_count = np.zeros(n_cr)
        n_acc = 0; rhat_trace, rhat_gens = [], []
        n_adapt = int(adapt_frac * n_gen)
        for t in range(1, n_gen + 1):
            props = np.empty_like(X); crs = np.empty(n_chains, int)
            zstd = Z[:nz].std(axis=0) + 1e-300
            for i in range(n_chains):
                delta = rng.integers(1, delta_max + 1)
                idx = rng.choice(nz, size=2 * delta, replace=False)
                diff = Z[idx[:delta]].sum(axis=0) - Z[idx[delta:]].sum(axis=0)
                k = rng.choice(n_cr, p=pcr); crs[i] = k
                mask = rng.random(d) < cr_vals[k]
                if not mask.any():
                    mask[rng.integers(d)] = True
                dprime = mask.sum()
                gamma = 1.0 if rng.random() < p_gamma1 else 2.38 / np.sqrt(2 * delta * dprime)
                lam = rng.uniform(-c, c, d)
                eps = rng.normal(0, cstar, d) * span
                step = np.where(mask, (1 + lam) * gamma * diff + eps, 0.0)
                p = X[i] + step
                # reflect at bounds (repeat once; clip as last resort)
                for _ in range(2):
                    p = np.where(p < lo, 2 * lo - p, p); p = np.where(p > hi, 2 * hi - p, p)
                props[i] = np.clip(p, lo, hi)
            lpp = np.array(pmap(list(props)))
            u = np.log(rng.random(n_chains))
            acc = (u < lpp - lp) & np.isfinite(lpp)
            if t <= n_adapt:                              # CR adaptation (burn-in only)
                dist = (((props - X) / zstd) ** 2).sum(axis=1)
                for i in range(n_chains):
                    cr_count[crs[i]] += 1
                    if acc[i]:
                        sq_jump[crs[i]] += dist[i]
                if t % 50 == 0 and cr_count.min() > 0 and sq_jump.sum() > 0:
                    w = sq_jump / cr_count; pcr = w / w.sum(); pcr = np.maximum(pcr, 0.02); pcr /= pcr.sum()
            X = np.where(acc[:, None], props, X); lp = np.where(acc, lpp, lp)
            n_acc += acc.sum()
            chains[t], logp[t] = X, lp
            if t % thin_K == 0:
                Z[nz:nz + n_chains] = X; nz += n_chains
            if t % rhat_every == 0 and t >= 2 * rhat_every:
                half = chains[t // 2:t + 1]
                r = gelman_rubin(half); rhat_trace.append(r); rhat_gens.append(t)
                if verbose and t % (rhat_every * 10) == 0:
                    print(f"  gen {t:6d}  max R-hat {r.max():.3f}  acc {n_acc / (t * n_chains):.2f}  best logp {logp[:t+1].max():.2f}", flush=True)
    finally:
        if pool:
            pool.close(); pool.join()
    return dict(chains=chains, logp=logp, accept_rate=n_acc / (n_gen * n_chains), rhat_trace=np.array(rhat_trace),
                rhat_gens=np.array(rhat_gens), archive=Z[:nz], cr_probs=pcr)
