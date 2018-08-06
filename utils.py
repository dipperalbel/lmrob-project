import numpy as np


def qr_qy(
    q=None, 
    r=None,
    X=None, 
    y):
    if X:
        #R: qr(X)$qr
        q,r = np.linalg.qr(X)
    #R: r * (row(r) == col(r))
    z = np.diag((np.diagonal(r)))  
    # R: Z = qr.qy(QR, z)
    Zq, Zr = np.linalg.qr(q)
    Z = np.matmul(Zq, y)
    return Z



def ghq(n=1, modify=True):
    n=int(n)
    if n<0:
        print('Need non-negative number of nodes')
        return
    elif n==0:
        return {'nodes'[]:, 'weights':[]}
    il = np.arange(1, n)
    muzero = np.sqrt(np.pi)

    b = np.sqrt(il/2)

    A = np.zeros(n**2)
    A[(n+1)*(il-1)+1]=b
    A[(n+1)*il-1]=b
    A = A.reshape(n,n)
    vd = np.linalg.eig(A,)
    n_rev = np.arange(n-1,-1,-1)
    x = vd[0]
    w = vd[1][0, n_rev,]
    w = np.sqrt(np.pi)*w*w
    if modify:
        w = w*np.exp(x**2)
    return {'nodes': x, 'weights': w}


#Not completed
def _lmrob_hat(
    x, 
    w=None, 
    wqr = None,
    names=True):
    if not w:
        w = [1]*x.shape[0]
    if not wqr:
        wqr = qr_func(np.sqrt(w)*x)
        #need from wqr: "compact" Q, row num of qr, rank
    # a = qr.qy(wqr, diag(1, NROW(wqr$qr), wqr$rank)) - R code
    h = np.minimum(1, np.sum(a**2, axis=1))
    return h