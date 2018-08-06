"""
C entry points:
    lmrob.R:
        C_Cdqrls
    lmrob.MM.R:
        R_psifun
        R_chifun
        R_wgtfun
        R_rho_inf

INPUT
    formula, data, subset : pointers towards the data to work with
    weights               : optional vector speciying vectors to be used in addition to robustness ones
    na.action             : action towards NAs
    method = 'MM'         : specification of estimator-chain. MM==SM, recommended - 'KS2014'
    model, x, y           : bool, what to return with results (model frame, model matrix, response respectively)
    singular.ok = TRUE,   : indicator if a singular fit is accepted
    contrasts = NULL,     :
    offset = NULL         : vector of offsets, added to the response during fitting (len==len(y))
    control = NULL        :  set of parameters for tuning, dict? 
    init = NULL           : initial estimator. string, function or list.
        str      - spec of built-in estimator ('S' or 'M-S')
        function - takes x, y, control, model.frame, returns initial coeffs
        list     - contains coeffs and scale
    ...)

OUTPUT
lmrob object, usable by other methods, having several values accessible
(R) - returns value used in the function call, missing if none provided
    coefficients  : estimate of the coeff vector
    scale         : scale as used in the M-estimator
    residuals     : residuals, associated with the estimator
    converged     : bool, if IRWLS (iteratively reweighted least squares) converged
    iter          : number of IRWLS iterations
    rweights      : 'robustness weights'
    fitted.values : fitted values, associated w/ the estimator
    init.S        : list, returned by lmrob.S or lmrob.M.S (for MM only)
    init          : list with results of intermediate estimates (not MM)
    rank          : rank of the fitted linear model
    cov           : estimated covariance matrix for the coeffs
    df.residual   : residual degrees of freedom
    weights       : (R)
    na.action     : (R)
    offset        : (R)
    contrasts     : (R), contrasts used
    xlevels       : record of the level of the factors used
    call          : matched call
    terms         : 'terms' object used
    model, x, y   : data as requested on function call

"""



import pandas as pd
import numpy as np
import numpy.linalg
import patsy
from itertools import chain
import scipy as sp
import scipy.linalg
import logging
import time

logging.basicConfig(format="[%(asctime)s] [%(levelname)s]: %(message)s",
                   level=logging.DEBUG,
                   filename='lmrob {}.log'.format(time.asctime()))


def lmrob(
    formula, # string in appropriate format
    data=None, # dict of named list-likes or pd.DataFrame
    subset=None, 
    weights=None, #list-like
    na_action=None, 
    method = 'MM',
    model = True, 
    x = False, 
    y = False, 
    singular_ok = True, 
    contrasts = None,
    offset = None, #list-like, must match y in length
    control = None, 
    init = None, #either a string {'M-S', 'S'} or func (TODO)
    *args):

    if not control:
        if method:
            control = create_control(method=method)
        else:
            control = create_control()
    return_x, return_y = x, y
    logging.info('Getting column names')
    if type(data)==pd.DataFrame:
        var_names=list(data.columns)
    elif type(data)==dict:
        var_names=list(data.keys())
    # basic functionality - construct data, response, response offset.
    logging.info('Preparing data')
    offset = np.array(offset)
    predictors_data, response_data = patsy.dmatrices(formula, data=data)
    assert len(offset)==response_data.shape[0], "Offset's length doesnt match response length"
    if weights:
        logging.info('Weights detected, preparing')
        weights_dict={name:weight for name, weight in zip(var_names, weights)}
        weights_transformed = np.array(*chain(patsy.dmatrix(formula.split('~')[1], data=weights_dict).to_list()))
        assert len(weights_transformed)==predictors_data.shape[1], "Weights length doesn't match number of predictors" 
        assert sum(weights_transformed<0)==0, 'Negative weights are not allowed'
        assert sum(weights_transformed==0)==0, 'Zero weights are not implemented'
    else:
        weights_transformed = np.ones(predictors_data.shape[1])
    #TODO allow zero weights
    #if the model is empty, error would be thrown by now, so we assume it is not.

    X = predictors_data*np.sqrt(weights_transformed)
    Y = response_data - offset.reshape(-1,1)

    #Data is prepared, beginning computation

    #singular fit
    #we need to find rank; 
    #TODO Tau from Household reflections are not currently implemented.
    K = min(X.shape)
    X_rank = np.linalg.matrix_rank(X)
    Q,R,pivot = sp.linalg.qr(X, pivoting=True, )
    is_singular = bool(X_rank-K)
    if is_singular&(not singular_ok):
        print('Singular fit encountered!\nInterrupting.')
        return
    if not X_rank:
        print('Zero rank matrix encountered!\nInterrupting.')


    #do initial estimation
    if init:
        if type(init)==str:
            assert init in ['M-S', 'S'], 'init must be either S or M-S'
            if init=='M-S':
                initialized = lmrob_M_S(X,Y,control,split)
                pass #do lmrob.m.s. here TODO
            elif init=='S':
                initialized = lmrob_S(X,Y,control)
                pass #do lmrob.s here TODO
        elif callable(init):
            print('Function initiation not implemented, passing')
            #TODO func initiation
        #check if coefficients and scale after init are numeric
        if (control['method']=='MM')|(control['method'][0]=='S'):
            control['method'] = control['method'][1:]
        """
        change control arguments:
        ## modify (default) control$method, possibly dropping first letter:
        if (control$method == "MM" || substr(control$method, 1, 1) == "S")
            control$method <- substring(control$method, 2)
        ## check for control$cov argument
        if (class(init)[1] != "lmrob.S" && control$cov == '.vcov.avar1')
            control$cov <- ".vcov.w"
        """






    results = {} 

    return results

def create_control(method=None):
    pass
def lmrob_M_S(X, Y, control, split):
    pass
def lmrob_S(X, Y, control):
    pass
def lmrob_fit(X, Y, control, initialized):
    pass

