import numpy as np

def crossing_times(P):
    k = len(P) // 2
    Q = P[:k, :k]
    t = np.linalg.solve(np.eye(k) - Q, np.ones(k)) # solves (I - Q)t = 1
    # check that expected values of the interior nodes give k^2
    assert np.allclose(t[:k-1], k**2)
    # check that expected value from the bridge node gives k^2 - k + 1
    assert np.isclose(t[k - 1], k**2 - k + 1)
    return t

def lambda2(P):
    """
    Returns the second largest by the absolute value eigenvalue
    Requires k >= 3
    """
    assert len(P) // 2 >= 3
    evs = np.linalg.eigvals(P).real
    return np.sort(abs(evs))[-2]

def exact_distance(P, pi, tol = 1e-3, start = 0, n_max = 100000):
    """Return delta(n) = d_TV(P^n(s_node, ), pi) until Delta <= tol"""
    mu = np.zeros(len(P))
    mu[start] = 1.0
    Delta = []
    while True:
        d = 0.5 * np.abs(mu - pi).sum()
        Delta += [d]
        if d <= tol:
            break
        mu = mu @ P
        assert len(Delta) < n_max, "no mixing or chain is periodic"
    return np.array(Delta)

def mixing_time(Delta, eps = 0.25):
    """First n with Delta(n) <= eps."""
    assert Delta[-1] <= eps, "distances history too short for this epsilon"
    return int(np.argmax(Delta <= eps))