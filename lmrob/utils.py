import numpy as np
import re
import patsy
import ctypes
import copy
import os
import sys
import inspect
import math
import copy
import scipy.integrate as integrate
from scipy.optimize import brentq
from scipy.stats import norm

_Mpsi_R_names  = ['bisquare', 'lqq', 'welsh', 'optimal', 'hampel', 'ggw']
_Mpsi_M_names  = ['huber']
_Mpsi_names  = _Mpsi_R_names + _Mpsi_M_names
_path_of_script = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe() ))[0]))
_LMROBLIBPATH = os.path.join(_path_of_script,"./liblmrob.so")
lmroblib = ctypes.cdll.LoadLibrary(_LMROBLIBPATH)


def wgtHighMedian(x, weights=None):
    """

    Compute the weighted Hi-Median of x.


    :Example:

    
    >>> module.wgtHighMedian([1,2,4,5,6,8,10,12])        

    :param x: numeric vector
    :type x: list
    :param weights: numeric vector of weights; of the same length as x.
    :type weights: list
    :return: The weighted Hi-Median of x.
    :rtype: float
    :raises: ValueError
    """
    # make np array
    # np.array(x) create an array from x and it assaign to  the variable x
    x = np.array(x)

    # sort x
    # after the array is created, it is sorted
    x_sort = np.argsort(x)
    x = x[x_sort]
    # parse weights and sort
    if (weights != None):
        if(len(x) != len(weights)):
            raise ValueError('x and weights cannot have different length')
    if(weights is None):
        weights = np.array([1]*len(x))
    weights = np.array(weights)
    weights = weights[x_sort]

    # cumSum and interpolate
    # cumsum return an array with the cumulative sum of the elements. The array obtain it is substracted with  0.5 *  the array weights
    cumW = np.cumsum(weights) - 0.5 * weights
    cumW /= np.sum(weights)

    # function interp
    res = np.interp(0.5, cumW, x)

    # find closest value
    abs_diff = np.absolute(np.array(x-res))
    abs_idx = np.where(abs_diff == np.min(abs_diff))
    fin = np.max(x[abs_idx])
    return(fin)

# function for mad (median absolute deviation)
def mad(x, center=None, constant=1.4826, na_rm=False, low=False, high=False):
    """
    Compute the median absolute deviation, i.e., the (lo-/hi-) median of the absolute deviations from the median, and (by default) adjust by a factor for asymptotically normal consistency.
    
    :Example:

    >>> module.mad([1,2,4,5,6,8,10,12])

    :param x: numeric vector.
    :type x: list
    :param center: optionally, the centre: defaults to the median.
    :type center: float
    :param constant: scale factor.
    :type constant: float
    :param na_rm: if True then NaN values are stripped from x before computation takes place.
    :type na_rm: bool
    :param low: if True, compute the ‘lo-median’, i.e., for even sample size, do not average the two middle values, but take the smaller one.
    :type low: bool
    :param high: if True, compute the ‘hi-median’, i.e., take the larger of the two middle values for even sample size.
    :type high: bool
    :return: the median absolute deviation
    :rtype: float
    :raises: ValueError
    """
    axis = None
    if (na_rm == True):
        i = np.isnan(x)
        x = np.array(x)[np.where(i == False)]
    if(center == None):
        center = np.median(x, axis=None)
    # Fin qua tutto apposto
    n = len(x)
    if ((low or high) and n % 2 == 0):
        if (low and high):
            raise ValueError('low and high cannot be both true')
        n2 = ((n//2) + high) - 1
        return constant * np.partition((np.absolute(np.array(x) - center)), n2)[n2]
    return constant*np.median(np.absolute(np.array(x) - center), axis)

def tauHuber(x, mu, k=1.5, s=None, resid=None):
    """
    Calculate correction factor Tau for the variance of Huber-M-Estimators.
    
    :Example:

    >>> module.tauHuber([1,2,4,5,6,8,10,12])

    :param x: numeric vector.
    :type x: list
    :param mu: location estimator.
    :type mu: float
    :param k: tuning parameter of Huber Psi-function.
    :type k: float
    :param s: scale estimator held constant through the iterations.
    :type s: float
    :param resid: the value of (x - mu)/s .
    :type resid: float
    :return: The correction factor Tau for the variance of Huber-M-Estimators.
    :rtype: float
    """
    # get needed parameters
    if(s == None):
        #s = mad(x)
        # mimic R default parameters: set constant to 1.4826
        s = mad(x, constant=1.4826)
    if(resid == None):
        resid = (np.array(x) - mu)/s
    # perform calculation
    inta = np.absolute(resid) <= k
    psi = np.copy(resid)
    psi[inta == False] = np.sign(resid[inta == False])*k
    psiP = np.bool_(inta)
    res = len(x) * np.sum(np.square(psi)) / np.square(np.sum(psiP))
    # return
    return(res)

def sumU(x, weights):
    """
    Calculate weighted sum
    
    :Example:

    >>> module.sumU([1,2,4,5,6,8,10,12], [52,44,82,24,68,42,82,20])
    
    :param x: numeric vector.
    :type x: list
    :param weights: numeric vector of weights; of the same length as x.
    :type weights: list
    :return: The weighted sum.
    :rtype: float
    :raises: ValueError
    """
    if(len(x) != len(weights)):
        raise ValueError('x and weights not of same length')
    sums = np.sum(np.array(x) * np.array(weights))
    return(sums)

def huberM(x,
           k=1.5,
           weights=None,
           tol=0.000001,
           mu=None,
           s=None,
           se=False,
           warn0scale=True):
    """
    (Generalized) Huber M-estimator of location with MAD scale.

    :Example:

    
    >>> module.huberM([1,2,4,5,6,8,10,12])



    :param x: numeric vector.
    :type x: list
    :param k: positive factor; the algorithm winsorizes at k standard deviations.
    :type k: float
    :param weights: numeric vector of non-negative weights of same length as x, or None.
    :type weights: list
    :param tol: convergence tolerance.
    :type tol: float
    :param mu: initial location estimator.
    :type mu: float
    :param s: scale estimator held constant through the iterations.
    :type s: float
    :param se: logical indicating if the standard error should be computed and returned (as SE component). Currently only available when weights is None.
    :type se: bool
    :param warn0scale: logical; if True, and s is 0 and len(x) > 1, this will be warned about.
    :type warn0scale: bool
    :return: The tuple (mu, s , it , se ) containing the location, scale parameters, number of iterations used and the se component.
    :rtype: tuple
    """
    # parse x
    x = np.array(x)
    # parse mu
    if(mu == None):
        if(weights == None):
            mu = np.median(x)
        else:
            mu = wgtHighMedian(x, weights)
    # parse s
    if(s == None):
        if(weights == None):
            #s = mad(x, center = mu)
            # mimic R default parameters: set constant to 1.4826
            s = mad(x, center=mu, constant=1.4826)
        else:
            s = wgtHighMedian(np.abs(x - mu), weights)

    # get only valid data
    if(np.sum(np.isnan(x)) > 0):
        i = np.isnan(x)
        x = np.array(x)[np.where(i == False)]
        if(weights != None):
            weights = np.array(weights)[np.where(i == False)]
    # get length
    n = len(x)
    # check weights
    if(weights == None):
        sumW = np.copy(n)
    else:
        if(np.sum(np.array(weights) <= 0) or len(weights) != n):
            raise ValueError('Something wrong with weights')
        sumW = np.sum(weights)
    it = 0
    # if sum of weights 0
    if(sumW == 0):
        # return(list(mu = NA., s = NA., it = it, se = NA.)) #R
        return((np.NaN, np.NaN, it, np.NaN))
    if(se and weights != None):
        raise ValueError(
            'Std.error computation not yet available for the case of weights')
    if(s <= 0):
        if(s < 0):
            raise ValueError('negative scale s')
        if(warn0scale == True):
            warnings.warn("scale 's' is zero -- returning initial 'mu'")
    else:
        while(1):
            it = it + 1

            y = np.minimum(np.maximum(mu - k * s, x), mu + k * s)
            if(weights == None):
                mu1 = np.sum(y)/sumW
            else:
                mu1 = sumU(y, weights)/sumW
            if(np.absolute(mu - mu1) < (tol * s)):
                break
            mu = np.copy(mu1)

    if(se):
        se2 = s * np.sqrt(tauHuber(x, mu, k, s)/n)
    else:
        se2 = np.NaN
    # return
    mu2 = mu.item(0)
    return((mu2, s, it, se2))



class psi_func:


    """
    The class "psi_func" is used to store ψ (psi) functions for M-estimation. In particular, an object of the class contains ρ(x) (rho), its derivative ψ(x) (psi), the weight function ψ(x)/x, and ﬁrst derivative of ψ, Dpsi = ψ0(x)
    """

    def __init__(self, rho, psi,wgt, Dpsi, Dwgt, tDefs, Erho, Epsi2, EDpsi, name, xtras):
        """Example of docstring on the __init__ method.

        The __init__ method may be documented in either the class level
        docstring, or as a docstring on the __init__ method itself.

        Either form is acceptable, but the two should not be mixed. Choose one
        convention to document the __init__ method and be consistent with it.

        Note:
            Do not include the `self` parameter in the ``Args`` section.

        Args:
            rho (function): the ρ() function. This is used to formulate the objective function; ρ() can be regarded as generalized negative log-likelihood.
            psi (function): ψ() is the derivative of ρ.
            wgt (function): The weight function ψ(x)/x.
            Dpsi (function): the derivative of ψ, Dpsi(x) = psi0(x).
            Dwgt (function): the derivative of the weight function.
            tDefs (list): named numeric vector oftuning parameterDefault values.
            Epsi2 (function): The function for computing E[ψ2(X)] when X is standard normal.
            EDpsi (function): The function for computing E[ψ0(X)] when X is standard normal.
            name (string): Name of ψ-function used for printing.
            xtras (string): Potentially further information.
        """
        self.rho = rho
        self.psi = psi
        self.wgt = wgt
        self.Dpsi = Dpsi
        self.Dwgt = Dwgt
        self.tDefs = tDefs
        self.Erho = Erho
        self.Epsi2 = Epsi2
        self.EDpsi = EDpsi
        self.name = name
        self.xtras = xtras
    def __str__(self):
        """ 
        The function returns the name of the ψ-function class when you want to print the class.

        Returns:
            The name of the ψ-function class.        
        """
        return str(self.name) + ' psi function'

def getStrVars(strTerms):
    strVars = []
    strCombs = []

    for f in strTerms:
        terms = f.split("*")
        if len(terms) > 1:
            strCombs.append(f)
        for term in terms:
            if not term in strVars:
                strVars.append(term)
    return strVars, strCombs

def model_frame_(formula, data, subset, weights, na_action, offset):
    #TODOS: 1. subset
    #       2. weights
    #       3. na_action: drop equivalent but dont validated
    #       4. offset
    if not formula is None:
        formula = "".join(formula.split())
        splitFormula = formula.split("~")
        if len(splitFormula) < 2:
            raise Exception("Maybe a response into the formula")

        vResponse = splitFormula[0]
        vIntercept = splitFormula[1].split("+")
        vIndependent = getIndependent(vIntercept)
    else:
        vResponse, vIntercept, vIndependent = [], [], []

    ind = 0
    if not data is None and len(vResponse) > 0:
        for v in vIndependent:
            vData = data.get(v, None)
            if not vData is None:
                if ind == 0:
                    nanColums = np.empty(vData.shape, dtype='bool')
                    nanColums.fill(False)
                nanColums += np.isnan(vData)
                ind += 1
            else:
                raise Exception("the variable %s dont exist" %(v))

        vData = data.get(vResponse, None)
        if not vData is None:
            nanColums += np.isnan(vData)
        else:
            raise Exception("the variable response %s dont exist" %(v))

        nRow = len([False for i in nanColums if not i])
        nCol = ind + 1
        dataIndependent = np.ones((nRow, nCol))
        dataResponse = np.ones((nRow, 1))

        for i, v in enumerate(vIndependent):
            dataIndependent[:, i + 1] = data.get(v)[nanColums == False]

        dataResponse[:, 0] = vData[nanColums == False]
        dataIndependent[:, 0] = dataResponse[:, 0]

        dictReturn = {
            "response" : np.asfortranarray(dataResponse),
            "independent" : np.asfortranarray(dataIndependent),
            "terms" : {
                "variableNames" : [vResponse] + vIndependent,
                "interceptVariables" : vIntercept
            },
            "weights" : None,
            "offset" : None
        }
        return dictReturn

def model_frame(formula, data, subset, weights, na_action, offset):
    # Make model -----
    Y, X = patsy.dmatrices(formula, data, NA_action=na_action, return_type='matrix')
    W = weights
    O = offset
    n = Y.shape[0]
    condition_weights = weights.size > 0 if isinstance(weights, np.ndarray) else len(weights) > 0
    condition_offset = offset.size > 0 if isinstance(offset, np.ndarray) else len(offset) > 0
    if condition_weights:
        y, W = patsy.dmatrices(formula + '+' + W, data, NA_action=na_action, return_type='matrix')
        i = W.design_info.column_name_indexes[weights]
        W = W[:, i].reshape(n, 1)
        if not subset is None:
            W = W[subset, :]

    if condition_offset:
        y, O = patsy.dmatrices(formula + '+' + O, data, NA_action=na_action, return_type='matrix')
        i = O.design_info.column_name_indexes[offset]
        O = O[:, i].reshape(n, 1)
        if not subset is None:
            O = O[subset, :]

    if not subset is None:
        Y = Y[subset, :]
        X = X[subset, :]

    # Get assign -----
    Xdi = X.design_info
    termNames = Xdi.column_names
    termNamesFiltred = ["".join(re.split(r"\[T\.[a-zA-Z0-9]*\]", term)) for term in termNames]
    strFactors = [term for i, term in enumerate(termNamesFiltred) if termNamesFiltred.index(term) == i]
    assign = [strFactors.index(term) for term in termNamesFiltred]
    assign = np.array(assign, dtype="int")

    # Get factors names -----
    strFactors = strFactors[1:]
    strVars = [factor for factor in strFactors if not ":" in factor]

    # Get equivalent to 'factors' -----
    factors = np.zeros((len(strVars), len(strFactors)), dtype="int")
    for i, var in enumerate(strVars):
        factors[ i, [var in factor for factor in strFactors] ] = 1

    # Get data Class -----
    Xfi = Xdi.factor_infos
    dataClass = [Xfi[patsy.EvalFactor(var)].type for var in strVars]
    dataClass = np.array(dataClass)

    dictReturn = {
        "response" : np.asfortranarray(Y),
        "model_matrix" : np.asfortranarray(X),
        "terms" : {
            "info" : Xdi,
            "factors" : factors,
            "dataClass" : dataClass,
            "assign" : assign
        },
        "weights" : W,
        "offset" : O
    }
    return dictReturn

def model_matrix_(_obj, _dict):
    #TODO: adapt to work with categorical variables
    x = _dict.get("independent")
    x[:, 0] = 1.0

    assign = np.arange(1, x.shape[1])
    dictReturn = {
        "x" : x,
        "assign" : assign
    }
    return dictReturn


def _psi_conv_cc(psi, cc):
    if psi.lower() == 'ggw':
        if np.all(cc == np.array([-.5, 1  , 0.95, np.nan])):
            return 1
        elif np.all(cc == np.array([-.5, 1  , 0.85, np.nan])):
            return 2
        elif np.all(cc == np.array([-.5, 1. , np.nan, 0.5])):
            return 3
        elif np.all(cc == np.array([-.5, 1.5, 0.95, np.nan])):
            return 4
        elif np.all(cc == np.array([-.5, 1.5, 0.85, np.nan])):
            return 5
        elif np.all(cc == np.array([-.5, 1.5, np.nan, 0.5])):
            return 6
        else:
            return cc
    elif psi.lower() == 'lqq':
        l = np.isnan(cc) == False
        if np.all(cc[l] == np.array([-.5, 1.5, 0.95, np.nan])[l]):
            return np.array([1.4734061, 0.9822707, 1.5])
        elif np.all(cc[l] == np.array([-.5, 1.5, np.nan, 0.5])[l]):
            return np.array([0.4015457, 0.2676971, 1.5])
        else:
            return cc
    return cc

def _psi2ipsi(psi):
    psi = _regularize_Mpsi(psi, redescending=False)
    dictData = ['huber', 'bisquare','welsh','optimal','hampel','ggw','lqq']
    i = dictData.index(psi)
    if i is None:
        raise Exception("internal logic error in psi() function name: ", psi)
    return i

def _regularize_Mpsi(psi, redescending=True):

    psi = psi.lower()
    psi = {"tukey": "bisquare", "biweight": "bisquare"}.get(psi, psi)

    nms = _Mpsi_R_names if redescending  else _Mpsi_names
    if psi in nms:
        i = nms.index(psi)
    else:
        raise

    return nms[i]

def _convSS(ss):
    if ss == "single":
        return 0
    elif ss == "nonsingular":
        return 1
    else:
        raise Exception("unknown setting for 'subsampling': %s" % ss)

def lmrob_rweights(resid, scale, cc, psi, eps = 16 * np.finfo(float).eps):
    if scale == 0:
        m = np.max(resid)
        if m == 0:
            return np.zeros_like(resid)
        else:
            return eps * m
    else:
        return Mwgt(resid / scale, cc, psi)


def lmrob_kappa(obj={"control":None}, control=None):
    
    if control is None:
        control = obj.get("control", None)
        if control is None:
            Exception('control is missing')
    controlKappa = control.copy()
    if controlKappa.method in ['S', 'SD']:
        controlKappa.tuning_psi = controlKappa.tuning_chi
    else:
        controlKappa.tuning_psi = controlKappa.tuning_psi
    gh = ghq(13 if controlKappa.numpoints is None else controlKappa.numpoints)
    psi = lambda r: Mpsi(r, controlKappa.tuning_psi, controlKappa.psi, deriv=0)
    wgt = lambda r: Mwgt(r, controlKappa.tuning_psi, controlKappa.psi)
    expr = lambda kappa: psi(gh.get("nodes")) * gh.get("nodes") - kappa * wgt(gh.get("nodes"))
    FF = lambda kappa: expr(kappa) * norm.pdf(gh.get("nodes"))
    lmrob_E = lambda kappa: np.sum(gh.get('weights') * FF(kappa))
    fun_min = lambda kappa: lmrob_E(kappa)
    return np.array(brentq(fun_min,0.1, 1))


def lmrob_tau(obj={"control": None, "x":None}, x=np.array([]), control=None, h=np.array([]), fast=True):
    if control is None:
        control = obj.get("control")
        if control is None:
            Exception("control' is missing")

    if h.size == 0:
        if obj.get("qr"):
            h = _lmrob_hat(wqr=obj.get("wqr"))
        else:
            h = _lmrob_hat(x, obj.get("rweights"))
    ## speed up: use approximation if possible
    if fast and not control.method in ['S', 'SD']:
        c_psi = control.tuning_psi
        tfact = None
        tcorr = None
        if control.psi == "optimal":
            if np.all(c_psi == 1.060158):
                tfact = 0.94735878
                tcorr = -0.09444537
        elif control.psi == "bisquare":
            if np.all(c_psi == 4.685061):
               tfact = 0.9473684
               tcorr = -0.0900833
        elif control.psi == "welsh":
            if np.all(c_psi == 2.11):
               tfact = 0.94732953
               tcorr = -0.07569506
        elif control.psi == "ggw":
            if np.array_equal(c_psi, np.array([-.5, 1.0, 0.95, None])):
                tfact = 0.9473787
                tcorr = -0.1143846
            elif np.array_equal(c.psi, np.array([-.5, 1.5, 0.95, None])):
                tfact = 0.94741036
                tcorr = -0.08424648
        elif control.psi == "lqq":
            l = np.isnan(c_psi) == False
            if np.array_equal(c_psi[l], np.array([-.5, 1.5, 0.95, np.nan])[l]):
                tfact = 0.94736359
                tcorr = -0.08594805
        elif control.psi == "hampel":
            if np.array_equal(c_psi, np.array([1.35241275, 3.15562975, 7.212868])):
               tfact = 0.94739770
               tcorr = -0.04103958
    if not tfact is None:
        return np.sqrt(1 - tfact * h) * (tcorr * h  + 1)
    
    ## else "non-fast" -- need to compute the integrals

    # kappa
    kappa =  obj.get("kappa", lmrob_kappa(obj, control))

    psi = control.psi
    if psi is None:
        Exception("parameter psi is not defined")

    cpsi = control.tuning_chi if control.method in ['S', 'SD']  else control.tuning_psi
    cpsi = _psi_conv_cc(psi, cpsi)# has its test
    ipsi = _psi2ipsi(psi)

    # constant for stderr of u_{-i} part and other constants
    inta = lambda r: _Mpsi(r, cpsi, ipsi) ** 2 * norm.pdf(r)
    intb = lambda r: _Mpsi(r, cpsi, ipsi, deriv = 1) * norm.pdf(r)

    ta, aux = integrate.quad(inta, -np.inf, np.inf)
    tb, aux = integrate.quad(intb, -np.inf, np.inf)

    hu = np.unique(h)
    nu = hu.size

    tau = np.zeros_like(hu)

    tc = ta / tb ** 2
    gh = ghq(control.numpoints)
    ghz = gh.get("nodes")
    ghw = gh.get("weights")

    for i, hui in enumerate(hu):
        s = np.sqrt(tc * (hui - hui ** 2))
        tc2 = hui / tb
        def func(w,v, sigma_i): 
            t = (v - tc2 * _Mpsi(v, cpsi, ipsi) + w * s) / sigma_i
            psi_t = _Mpsi(t, cpsi, ipsi)
            return (psi_t * t - kappa * psi_t/t) * norm.pdf(v) * norm.pdf(w)

        def wint(v, sigma_i):
            return np.array(list(map(lambda v_j: np.sum(func(ghz, v_j, sigma_i) * ghw), v)))

        def vint(sigma_i):            
            return np.sum(wint(ghz, sigma_i) * ghw)
        tau[i] = brentq(vint, 3.0/20.0 if (hui < 0.9) else 1.0/16.0, 1.1)
    return tau[[np.where(hu == i)[0][0] for i in h]]



def Mwgt(x, cc, psi):
    byref = ctypes.byref
    ipsi = ctypes.c_int(_psi2ipsi(psi))
    c_chi = _psi_conv_cc(psi, cc)
    c_chi_p = c_chi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    wgt = lmroblib.wgt
    wgt.argtypes = [ctypes.c_double, 
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.c_int]
    wgt.restype = ctypes.c_double
    rweights = np.array([wgt(res, c_chi_p, ipsi) for res in x])
    return rweights


def Mpsi(x, cc, psi, deriv):
    byref = ctypes.byref
    ipsi = ctypes.c_int(_psi2ipsi(psi))
    #c_chi = ctypes.c_double(_psi_conv_cc(psi, cc))
    c_chi = _psi_conv_cc(psi, cc)
    c_chi_p = c_chi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    kernel_switch = {-1: lmroblib.rho,
                      0: lmroblib.psi,
                      1: lmroblib.psip,
                      2: lmroblib.psi2
                    }
    kernel = kernel_switch.get(deriv)
    kernel.argtypes = [ctypes.c_double, 
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.c_int]
    kernel.restype = ctypes.c_double
    values = np.array([kernel(res, c_chi_p, ipsi) for res in x])
    if ipsi.value != 0 and deriv == -1:
        rho_inf = lmroblib.rho_inf
        rho_inf.argtypes = [ctypes.POINTER(ctypes.c_double),
                            ctypes.c_int]
        rho_inf.restype = ctypes.c_double
        return np.array([kernel(res, c_chi_p, ipsi) * rho_inf(c_chi_p, ipsi) for res in x])
    else:
        return np.array([kernel(res, c_chi_p, ipsi) for res in x])


def _Mpsi(x, c_chi_, ipsi_, deriv=0):
    if not isinstance(x, np.ndarray):
        x = np.array([x])
    byref = ctypes.byref
    ipsi = ctypes.c_int(ipsi_)
    c_chi_p = c_chi_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    kernel_switch = {-1: lmroblib.rho,
                      0: lmroblib.psi,
                      1: lmroblib.psip,
                      2: lmroblib.psi2
                    }
    kernel = kernel_switch.get(deriv)
    kernel.argtypes = [ctypes.c_double, 
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.c_int]
    kernel.restype = ctypes.c_double
    
    values = np.array([kernel(res, c_chi_p, ipsi) for res in x])
    if ipsi.value != 0 and deriv == -1:
        rho_inf = lmroblib.rho_inf
        rho_inf.argtypes = [ctypes.POINTER(ctypes.c_double),
                            ctypes.c_int]
        rho_inf.restype = ctypes.c_double
        return np.array([kernel(res, c_chi_p, ipsi) * rho_inf(c_chi_p, ipsi) for res in x])
    else:
        return np.array([kernel(res, c_chi_p, ipsi) for res in x])


def _lmrob_hat(x, w=None, wqr=None, names=True):
    nrow = x.shape[0]
    if w is None:
        w = np.ones((nrow, 1))
    if wqr is None:
        wqr = np.linalg.qr(np.sqrt(w.reshape(nrow, 1)) * x)
    q, r = wqr
    pmin = lambda x: min(x,1)
    h = np.array(list(map(pmin, np.linalg.norm(q, axis=1)**2)))
    return h


def Mchi(x, cc, psi, deriv):
    byref = ctypes.byref
    byref = ctypes.byref
    ipsi = ctypes.c_int(_psi2ipsi(psi))
    c_chi = ctypes.c_double(_psi_conv_cc(psi, cc))
    kernel_switch = { 0: lmroblib.rho,
                      1: lmroblib.psi,
                      2: lmroblib.psip,
                      3: lmroblib.psi2
                    }
    kernel = kernel_switch.get(deriv)
    kernel.argtypes = [ctypes.c_double, 
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.c_int]
    kernel.restype = ctypes.c_double
    values = np.array([kernel(res, byref(c_chi), ipsi) for res in x])
    if deriv > 0:
        rho_inf = lmroblib.rho_inf
        rho_inf.argtypes = [ctypes.POINTER(ctypes.c_double),
                            ctypes.c_int]
        rho_inf.restype = ctypes.c_double
        return np.array([kernel(res, byref(c_chi), ipsi) / rho_inf(byref(c_chi), ipsi) for res in x])
    else:
        return np.array([kernel(res, byref(c_chi), ipsi) for res in x])


def split_frame(mf, x=None, _type=None):
    if x is None:
        x_ = mf.get("model_matrix")
    else:
        x_ = copy.deepcopy(x)
    if _type is None:
        type_ = "f"
    else:
        type_ = copy.deepcopy(_type)

    mt = mf.get("terms")
    p = x_.shape[1]
    ## --- split categorical and interactions of categorical vars.
    ##     from continuous variables
    factors = mt.get("factors")
    factor_idx = mt.get("dataClass") == "categorical"
    #factor_idx = np.array([_class == "categorical" for _class in dataClass])
    if (not np.any(factor_idx)):
        dictReturn = {
            "x1_idx" : np.ndarray(shape=p, dtype=bool) * False,
            "x1" : np.ndarray(shape=(x.shape[0],0), dtype="float") * np.nan,
            "x2" : x_
        }
        return dictReturn

    ## --- include interactions cat * cont in x1:
    if type_ == "fi":
        factor_asgn = np.where( np.dot(factor_idx, factors) > 0 )[0]
    ## --- include also continuous variables that interact with factors in x1:
           ##     make sure to include interactions of continuous variables
           ##     interacting with categorical variables, too
    elif type_ == "fii":
        factor_asgn = np.ndarray(shape=0)
        factors_cat = copy.deepcopy(factors)
        factors_cat[factors_cat > 0] = 1 ## fix triple+ interactions
        factors_cat[:, np.dot(factor_idx, factors) == 0] = 0
        for i in range(factors.shape[1]):
            comp = factors[:, i] > 0
            ## if any of the components is a factor: include in x1 and continue
            if np.any(factor_idx[comp]):
                factor_asgn = np.append(factor_asgn, i)
            ## if there is an interaction of this term with a categorical var.
            else:
                tmp = np.sum(factors_cat[comp, :], 0) >= np.sum(comp)
                if np.any(tmp):
                    ## if no other continuous variables are involved
                    ## include in x1 and continue
                    ## if (identical(factors[!comp, tmp], factors.cat[!comp, tmp]))
                    not_comp = comp == False
                    if (not np.all( np.sum(factors[not_factor.idx * not_comp, tmp]) > 0 ) ):
                        factor_asgn = factor_asgn = np.append(factor_asgn, i)
    ## --- do not include interactions cat * cont in x1:        
    elif type_ == "f":
        not_factor_idx = factor_idx == False
        factor_asgn = np.where( (np.dot(factor_idx, factors) > 0) * (np.dot(not_factor_idx, factors) == False) )[0]
    else:
        raise Exception("unknown split type")

    factor_asgn = np.append(0, factor_asgn + 1)
    x1_idx = np.array([i in factor_asgn for i in mt.get("assign")])

    ## x1: factors and (depending on type) interactions of / with factors
    ## x2: continuous variables
    dictReturn = {
            "x1" : x_[:, x1_idx],
            "x1_idx" : x1_idx,
            "x2" : x_[:, x1_idx == False]
        }
    return dictReturn


def ghq(n=1, modify=True):
    ## Adapted from gauss.quad in statmod package
    ## which itself has been adapted from Netlib routine gaussq.f
    ## Gordon Smyth, Walter and Eliza Hall Institute

    n = int(n)
    if n < 0:
        Exception("need non-negative number of nodes")
    if n == 0:
        return  {"nodes":[], "weights": []}
    ## i <- seq_len(n) # 1 .. n
    i1 = np.arange(1, n)

    muzero = math.sqrt(math.pi)
    ## a = numeric(n)
    b = np.sqrt(i1 / 2)

    A = np.zeros((n, n))
    ## A[(n+1)*(i-1)+1] = a # already 0
    I, J = np.diag_indices(n - 1)
    A[I, J + 1] = b
    A[I + 1, J] = b
    x, w = np.linalg.eigh(A)
    w = muzero * (w[0,::-1] ** 2)

    return {"nodes":x, "weights": w * np.exp(x ** 2) if modify  else w}

def _vcov_avar1(obj, x=None):
    return obj.get("x")

def _vcov_w(obj, x=None):
    return obj.get("x")