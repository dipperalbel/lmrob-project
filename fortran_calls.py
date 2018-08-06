from numpy import f2py



def rllarsbi_call(
    X,#matrix 
    Y):#matrix
    with open("rllarsbi.f", 'rb') as sourcefile:
        sourcecode = sourcefile.read()
    f2py.compile(sourcecode, modulename='rllarsbi', verbose=False)
    import rllarsbi
    #in case it is patsy.dmatrix
    x = X*1
    y = Y*1
    args={'x':x,
        'y':y,
        'tol':1e-7,
        'nit':0,
        'k':0,
        'kode':0,
        'sigma':0.0,
        'theta':np.zeros(x.shape[0]),
        'rs':np.zeros(x.shape[0]),
        'sc1':np.zeros(x.shape[0]),
        'sc2':np.zeros(x.shape[1]),
        'sc3':np.zeros(x.shape[1]),
        'sc4':np.zeros(x.shape[1]),
        'bet0':0.773372647623,
        'n':x.shape[0],
        'np':x.shape[1],
        'mdx':x.shape[0],
        'mdt':x.shape[0]}

    rllarsbi.rllarsbi(**args)#ARGS ARE CHANGED INPLACE
    return args