import numpy as np

def one_replicate_d(pi, nbrs, rng, N, n_max, start = 0):
    """Return array of D_hat(n)'s for 0 <= n <= n_max, simulating the walk of 
    N independent walkers with starting node = start
    """
    n_nodes = len(pi)
    D_hat = np.zeros(n_max + 1)
    pos = np.full(N, start)

    for n in range(n_max + 1):
        # keep track of how many walkers are there in each of the nodes
        counts = np.bincount(pos, minlength=n_nodes) 
        # use now that value to get estimate for mu
        mu_hat = counts / N
        D_hat[n] = 0.5 * np.abs(mu_hat - pi).sum()
        if n == n_max:
            break
        for i in range(N):
            next_pos = rng.integers(len(nbrs[pos[i]]))
            pos[i] = nbrs[pos[i]][next_pos]
    return D_hat

def one_replicate_cross(nbrs, rng, start = 0):
    """Returns T_cross for the start at node = start, simulating the walk"""
    k = len(nbrs) // 2
    x = start 
    t = 0
    while x < k: # while walker still in community A
        next_pos = rng.integers(len(nbrs[x]))
        x = nbrs[x][next_pos]
        t += 1
    return t

def t_cross_sample(nbrs, N, rng, start=0):
    """Returns array of N independent T_cross starting at start"""
    return np.array([one_replicate_cross(nbrs, rng, start) for _ in range(N)])

def dist_results(pi, nbrs, rng, N, R, n_max, start=0):
    """
    Returns mean of Delta_hat(n) and its SE over R independent 
    replicates at each n
    """
    walks = np.array([one_replicate_d(pi, nbrs, rng, N, n_max, start) for _ in range(R)])
    mean = walks.mean(axis = 0)
    se = walks.std(axis = 0, ddof=1)/np.sqrt(R) # set ddof = 1 for sample sd
    return mean, se
