import numpy as np

def build_barbell(k):
    n = 2 * k
    A = np.zeros((n, n))
    #first community
    for i in range(k):
        A[i, :k] = 1
        A[i, i] = 0

    #second community
    for i in range(k, n):
        A[i, k:] = 1
        A[i, i] = 0

    #bridge
    A[k - 1, k] = 1
    A[k, k - 1] = 1

    return A

def transition_matrix(A):
    """Returns transition matrix P and the stationary vector pi"""
    deg = A.sum(axis = 1)
    P = A / deg[:, None] # dividing row v by deg(v)
    pi = deg / deg.sum() # deg.sum() = 2 |E|
    return P, pi

def build_neighbor_list(A):
    return [np.flatnonzero(A[v]) for v in range(len(A))]

def check_graph(A, P, pi, k):
    """Sanity checks"""
    assert np.array_equal(A, A.T), "A not symmetric"
    assert np.all(np.diag(A) == 0), "self-loops"
    assert A.sum() / 2 == k * (k - 1) + 1, "wrong number of edges"
    assert np.allclose(P.sum(axis=1), 1), "rows of P do not sum to 1"
    assert np.allclose(pi @ P, pi), "pi not stationary"
    assert np.allclose(pi[:, None] * P, (pi[:, None] * P).T), "detailed balance fails"
