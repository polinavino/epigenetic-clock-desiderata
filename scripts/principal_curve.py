"""
principal_curve.py

Implements the Hastie-Stuetzle (1989) principal curve algorithm.

Reference:
    Hastie, T. and Stuetzle, W. (1989). Principal Curves.
    Journal of the American Statistical Association, 84(406), 502-516.
"""

import numpy as np
from sklearn.decomposition import PCA
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import uniform_filter1d


def fit_principal_curve(X, n_iter=10, smoothing_factor=None, random_state=42):
    """
    Fit a principal curve to X using the Hastie-Stuetzle algorithm.

    Parameters
    ----------
    X : (n_samples, n_features) array, centered
    n_iter : int
    smoothing_factor : float or None
        Passed to UnivariateSpline as s. None = auto.
    random_state : int

    Returns
    -------
    curve_points : (n_samples, n_features)
        Projection of each point onto the curve.
    tau : (n_samples,)
        Arc-length parameter per point.
    converged : bool
    """
    n, p = X.shape
    converged = False

    # Initialization: arc-length = projection onto first PC
    pca = PCA(n_components=1, random_state=random_state)
    t = pca.fit_transform(X)[:, 0]   # (n,)

    def make_splines(t, X, sf):
        """
        Fit one smoothing spline x_j(t) per dimension j.
        Returns list of fitted UnivariateSpline objects and sort index.
        """
        sort_idx = np.argsort(t)
        t_sorted = t[sort_idx]
        splines  = []
        for j in range(p):
            x_sorted = X[sort_idx, j]
            s = sf if sf is not None else max(n * x_sorted.var(), 1e-6)
            try:
                spl = UnivariateSpline(t_sorted, x_sorted, s=s, ext=3)
            except Exception:
                # Fallback: fit linear spline (k=1)
                try:
                    spl = UnivariateSpline(t_sorted, x_sorted, k=1, s=0, ext=3)
                except Exception:
                    spl = None
            splines.append(spl)
        return splines, sort_idx, t_sorted

    def eval_splines(splines, t_sorted, t_query):
        """Evaluate all splines at query points t_query -> (len(t_query), p)."""
        out = np.zeros((len(t_query), p))
        for j, spl in enumerate(splines):
            if spl is not None:
                out[:, j] = spl(t_query)
            else:
                # constant fallback
                out[:, j] = 0.0
        return out

    def project_points(X, splines, t_sorted, t_current):
        """
        For each point X[i], find the t* that minimizes ||X[i] - f(t*)||^2
        by evaluating the curve on a fine grid and taking the argmin.
        Then refine with a local search.
        """
        # Evaluate curve on a fine grid spanning current t range
        t_min, t_max = t_sorted.min(), t_sorted.max()
        t_grid = np.linspace(t_min, t_max, min(5 * n, 2000))
        curve_grid = eval_splines(splines, t_sorted, t_grid)  # (grid, p)

        t_new = np.zeros(n)
        for i in range(n):
            dists = ((curve_grid - X[i]) ** 2).sum(axis=1)   # (grid,)
            t_new[i] = t_grid[np.argmin(dists)]
        return t_new

    # ── Main loop ─────────────────────────────────────────────────────────────
    for iteration in range(n_iter):
        t_old = t.copy()

        # Smoothing step: fit splines f_j(t) for each dimension
        splines, sort_idx, t_sorted = make_splines(t, X, smoothing_factor)

        # Projection step: find t* for each point
        t = project_points(X, splines, t_sorted, t)

        delta = np.mean(np.abs(t - t_old))
        print(f"    Iteration {iteration+1}/{n_iter}  delta={delta:.8f}")
        if delta < 1e-4:
            print(f"    Converged.")
            converged = True
            break

    # Compute final curve points (projections)
    splines, sort_idx, t_sorted = make_splines(t, X, smoothing_factor)
    curve_points = eval_splines(splines, t_sorted, t)

    return curve_points, t, converged


def test_principal_curve():
    """Smoke test on synthetic helical data."""
    from scipy.stats import pearsonr
    import numpy as np

    np.random.seed(42)
    t_true = np.linspace(0, 4, 300)
    X_test = np.column_stack([
        t_true + np.random.randn(300) * 0.1,
        np.sin(t_true * 2) + np.random.randn(300) * 0.1,
        np.cos(t_true) + np.random.randn(300) * 0.1,
    ])
    X_test -= X_test.mean(axis=0)

    print("Testing principal curve on synthetic data (300 x 3)...")
    curve, tau, converged = fit_principal_curve(X_test, n_iter=10)

    r, _ = pearsonr(tau, t_true)
    print(f"  Tau correlation with true parameter: r={abs(r):.3f}")
    print(f"  Converged: {converged}")
    print(f"  {'PASS' if abs(r) > 0.9 else 'FAIL: r < 0.9'}")
    return abs(r) > 0.9


if __name__ == "__main__":
    test_principal_curve()
