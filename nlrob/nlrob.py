################################################################################
# nlrob
#
# This is a translation from nlrob function of robustbase packages in R
#
# Includes the following functions:
#
#     - nlrob
#     - nlrob_MM
#     - nlrob_tau
#     - nlrob_CM
#     - nlrob_mtl
#     - nls
#
#  Version 001 
#       - Date 20180931
#
#
################################################################################

from .utils import *
from .utils import  JDEoptim, _psi_conv_cc, _psi2ipsi, _regularize_Mpsi, _convSS, _Mpsi, _Mwgt
from argparse import Namespace 
from six import string_types



class NlrobControl():
    """
        Class that contains the parameters options of the lmrob functions
    """
    def __init__(self,
                 method,
                 psi ="bisquare",
                 init = "S",
                 optimizer = "JDEoptim",
                 fnscale=None,
                 tuning_chi_tau=None,
                 tuning_chi_scale=None,
                 tuning_chi=None,
                 cutoff=2.5,
                 *args, **kwargs
                ):

        _Mchi_tuning_defaults = {

            ## Here, psi must be redescending! -> 'huber' not possible
            'bisquare': np.array([1.54764]),
            'welsh': np.array([0.5773502]),
            'ggw': np.array([-0.5, 1.5, np.nan, 0.5]),
            'lqq': np.array([-0.5, 1.5, np.nan, 0.5]),
            'optimal': np.array([0.4047]),
            'hampel': np.array([1.5, 3.5, 8]) * 0.2119163
        }

        _Mpsi_tuning_defaults = {
              'huber':np.array([1.345]),
              'bisquare':np.array([4.685061]),
              'welsh':np.array([2.11]),
              'ggw':np.array([-0.5, 1.5, .95, np.nan]),
              'lqq':np.array([-0.5, 1.5, .95, np.nan]),
              'optimal':np.array([1.060158]),
              'hampel':np.array([1.5, 3.5, 8]) * 0.9016085
        }
        self.tuning_psi_M = None
        self.psi = psi
        if method == "M":
            self.method = method
        elif method == "MM":
            self.method = method
            self.init = init
            self.psi = psi
            self.tuning_chi_scale = _psi_conv_cc(psi, _Mchi_tuning_defaults[psi])
            self.tuning_psi_M = _psi_conv_cc(psi, _Mpsi_tuning_defaults[psi])
            self.optimizer = optimizer
            self.fnscale=fnscale
        elif method == "tau":
            self.method = method
            self.psi = psi
            self.tuning_chi_tau = tuning_chi_tau if tuning_chi_tau else None
            self.tuning_chi_scale = tuning_chi_scale if tuning_chi_scale else None
            self.fnscale = fnscale if fnscale else None
        elif method == "CM":
            self.method = method
            self.psi = psi
            self.tuning_chi = tuning_chi if tuning_chi else None
            self.fnscale = fnscale if fnscale else None
        elif method == "mtl":
            self.method = method
            self.fnscale = fnscale if fnscale else None
            self.cutoff = cutoff if cutoff else 2.5
        else:
            raise Exception("Method %s not correctly supported yet" % method)
            
    def copy(self):
        return copy.copy(self)

    def __str__(self):
        if self.method == "MM":
            string = "self.method = {:}\n".format(self.method)
            string += "self.init = {:}\n".format(self.init)
            string += "self.psi = {:}\n".format(self.psi)
            string += "self.tuning_chi_scale  = {:}\n".format(self.tuning_chi_scale)
            string += "self.tuning_psi_M  = {:}\n".format(self.tuning_psi_M)
            string += "self.optimizer = {:}\n".format(self.optimizer)
            string += "self.optArgs = {:}\n".format(self.optArgs)

            return string



def _Mwgt_psi1(psi, cc=None):
    global deriv
    if cc is None:
        cc = _Mpsi_tuning_default[psi]
    ipsi = _psi2ipsi(psi)
    ccc = _psi_conv_cc(psi, cc)
    def return_func(x, deriv=0):
        if deriv:
            return _Mpsi(x, ccc, ipsi, deriv)
        else:
            return _Mwgt(x, ccc, ipsi)
    return return_func



def nlrob(formula, data, start=np.zeros(1),
          lower=np.array([-np.Inf]),
          upper=np.array([np.Inf]),
          weights = None,
          method = "MM",
          psi=None,
          scale = None,
          control=None,
          test_vec = "resid",
          maxit = 20,
          tol = 1e-06,
          algorithm = "lm", doCov=False, trace=False):

    """
        Fits a nonlinear regression model by robust methods. Per default, by an M-estimator, using iterated reweighted least squares (called “IRLS” or also “IWLS”).
        This function returns a dictionary with the results

    Parameters
    ----------

    formula: str
        A nonlinear formula including variables and parameters of the model,
        such as y ~ f(x, theta) (cf. nls). (For some checks: if f(.) is 
        linear, then we need parentheses, e.g., y ~ (a + b * x)
    data: pandas.core.frame.DataFrame
        Data frame containing the variables in the model. If not found in
        data, the variables are taken from environment(formula), typically
        the environment from which nlrob is called.
    start: pandas.core.frame.DataFrame
        A named numeric vector of starting parameters estimates, only for
        method = "M".
    lower: pandas.core.frame.DataFrame
        numeric vectors of lower and upper bounds; if needed, will be
        replicated to be as long as the longest of start, lower or upper. For
        (the default) method = "M", if the bounds are unspecified all
        parameters are assumed to be unconstrained; also, for method "M",
        bounds can only be used with the "port" algorithm. They are ignored,
        with a warning, in cases they have no effect.
        For methods "CM" and "mtl", the bounds must additionally have an entry
        named "sigma" as that is determined simultaneously in the same
        optimization, and hence its lower bound must not be negative.
    upper:  array_like
        numeric vectors of lower and upper bounds; if needed, will be
        replicated to be as long as the longest of start, lower or upper. For
        (the default) method = "M", if the bounds are unspecified all
        parameters are assumed to be unconstrained; also, for method "M",
        bounds can only be used with the "port" algorithm. They are ignored,
        with a warning, in cases they have no effect.
    weights: arrray_like
        An optional vector of weights to be used in the fitting process (for
        intrinsic weights, not the weights w used in the iterative (robust)
        fit). I.e., sum(w * e^2) is minimized with e = residuals,
        e[i] = y[i] - f(xreg[i], theta), where f(x, theta) is the nonlinear
        function, and w are the robust weights from resid * weights.
    method: str
        a character string specifying which method to use. The default is "M", for historical and back-compatibility reasons. For the other methods, primarily see nlrob.algorithms.
        "M"
        Computes an M-estimator, using nls(*, weights=*) iteratively (hence, IRLS) with weights equal to ψ(r_i) / r_i, where r_i is the i-the residual from the previous fit.
        "MM"
        Computes an MM-estimator, starting from init, either "S" or "lts".
        "tau"
        Computes a Tau-estimator.
        "CM"
        Computes a “Constrained M” (=: CM) estimator.
        "mtl"
        Compute as “Maximum Trimmed Likelihood” (=: MTL) estimator.
        Note that all methods but "M" are “random”, hence typically to be preceded by set.seed() in usage.
    psi: func 
        A function of the form g(x, 'tuning constant(s)', deriv) that for deriv=0 returns psi(x)/x and for deriv=1 returns psi'(x). Note that tuning constants can not be passed separately, but directly via the specification of psi, typically via a simple _Mwgt_psi1() call as per default.
    scale: float
        When not None, a positive number specifying a scale kept fixed during
        the iterations (and returned as Scale component).
    test_vec: str
        Character string specifying the convergence criterion. The relative
        change is tested for residuals with a value of "resid" (the default),
        for coefficients with "coef", and for weights with "w".
    maxit: int
        maximum number of iterations in the robust loop.
    tol: float
        non-negative convergence tolerance for the robust fit.        
    algorithm: str
        character string specifying the algorithm to use for nls, see there, 
        only when method = "M". The default algorithm is a Gauss-Newton 
        algorithm.
    doCov: bool
        a logical specifying if nlrob() should compute the asymptotic 
        variance-covariance matrix (see vcov) already. This used to be 
        hard-wired to TRUE; however, the default has been set to FALSE, as vcov
        (obj) and summary(obj) can easily compute it when needed.
    control: obj
        An optional object of control settings.
    trace: bool
        logical value indicating if a “trace” of the nls iteration progress 
        should be printed. Default is False. 
        If True, in each robust iteration, the residual sum-of-squares and the 
        parameter values are printed at the conclusion of each nls iteration.

    Returns
    -------

    coefficients: array_like
        Coefficients of the regressor
    residuals: array_like
        Difference between the real values and the fitted_values
    fitted_values: array_like
        Estimated values by th regressor
    

    Examples
    --------
    >>> import pandas
    >>> from nlrob import *
    >>> 
    >>> formula = "density ~ Asym/(1 + np.exp(( xmid - np.log(conc) )/scal))"
    >>> lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
    ...         index=[0])
    >>> upper = np.array([1])
    >>> data = pandas.read_csv("Function=NLROBInput.csv")
    >>> # M Method
    >>> method = "M"
    >>> Rfit = nlrob(formula, data, lower, lower, upper, method=method)


    See Also
    --------
    nlrob_MM:
    nlrob_tau:
    nlrob_CM:
    nlrob_mtl:
    """
    
    hasWgts = not weights is None
    if method != "M":
        control = NlrobControl(method)
        if (hasWgts):
            raise Exception("specifying 'weights' is not yet supported for method %s " % method)
        if not psi is None:
            print("For method = \"%s\", currently 'psi' must be specified via 'control'" % method )
        def fixAns(mod):
            ctrl  = mod.get("ctrl")
            if isinstance(ctrl.psi, string_types) and isinstance(ctrl.tuning_psi_M, (int,float)):
                psi = _Mwgt_psi1(ctrl.psi, ctrl.tuning_psi_M)
                res_sc = mod.get("residuals") / mod.get("Scale")
                mod.update({"psi":psi})
                mod.update({"w": psi(res_sc)})
                mod.update({"rweights": psi(res_sc)})
            return mod

        if method == "MM":
            return fixAns(nlrob_MM(formula, data, lower=lower, upper=upper,
                                    tol=tol, ctrl=control))

        elif method == "CM":
            return fixAns(nlrob_CM(formula, data, lower=lower, upper=upper,
                                    tol=tol, ctrl=control))

        elif method == "tau":
            return fixAns(nlrob_tau(formula, data, lower=lower, upper=upper,
                                    tol=tol, ctrl=control))

        elif method == "mtl":
            return fixAns(nlrob_CM(formula, data, lower=lower, upper=upper,
                                    tol=tol, ctrl=control))

    else:
        psi = _Mwgt_psi1("huber", cc=1.345)
        updateScale = scale is None
        if not updateScale:
            if isinstance(scale, (float, int)) and scale > 0:
                Scale = scale
            else:
                raise Exception("'scale' must be NULL or a positive number")

        if hasWgts and np.any(weights < 0 or np.isnan(weights).any()):
            raise Exception("'weights' must be nonnegative and not contain NAs")
        data_names = data.keys()
        for var in data_names:
            if re.search("\s", var.strip()):
                continue
            globals().update({var:data[var].values})
        pnames = start.keys()
        p0 = start.values[0]

        y = data[formula.split("~")[0].rstrip()].values
        nobs = y.size
        right_hand_term = formula.split("~")[1]
        for i, pname in enumerate(pnames):
            right_hand_term = right_hand_term.replace(pname, "p0[%d]" % i)
        fit = eval(right_hand_term)
        resid = y - fit
        iris_delta = lambda old, new: np.sqrt(np.sum((old - new)**2)/np.max((1e-20, np.sum(old ** 2))))
        converged = False
        method_exit = False
        status = "converged"
        for iiter in range(maxit):
            if trace:
                print("robust iteration")
            previous = eval(test_vec)
            if updateScale:
                Scale = np.median(np.abs(resid)) / 0.6745
            if Scale == 0:
                convi = 0
                method_exit = True
                status =  "could not compute scale of residuals"
                print(status)
            else:
                w = psi(resid/Scale)
                if hasWgts: 
                    w = w * weights
                data.update({"_nlrob_w":w})
                out = nls(formula=formula, data=data, start=start, algorithm=algorithm,
                        lower=lower, upper=upper)
                coef = out.get("coefficients")
                resid = out.get("residuals")
                convi = iris_delta(previous, eval(test_vec))

            converged = convi <= tol
            if converged:
                break
            elif trace:
                print(" --> irls.delta(previous, %s) = %g -- *not* converged\n" % (test_vec, convi))
                
        if not converged or method_exit:
            st = "failed to converge in %d steps" % maxit
            print(st)
            status = st
        
        if hasWgts:
            tmp = weights != 0
            w[tmp] = w[tmp] / weights[tmp]
        res_sc = resid / Scale
        rw = psi(res_sc)

        if not converged or not doCov:
            asCov = None
        else:
            AtWAinv = np.linalg.inv(out.get("cov"))
            tau = np.mean(rw ** 2) / np.mean(psi(res_sc)) ** 2
            asCov = AtWAinv * Scale ** 2 * tau        
        dictReturn = {"coefficients": coef,
                      "formula": formula,
                      "nobs":nobs,
                      "residuals": resid,
                      "fitted_values": fit,
                      "Scale": Scale,
                      "w": w,
                      "rweights": rw,
                      "cov": asCov,
                      "test_vec": test_vec,
                      "status": status,
                      "iter": iiter,
                      "psi": psi,
                      "data": data
                      }

        return dictReturn


    

def nlrob_MM(formula, data, lower, upper, tol=1e-6, psi="bisquare", init="S",
             ctrl=NlrobControl("MM")):
    """
    Compute an MM-estimator for nonlinear robust (constrained) regression.
    Returns a dictionary with all the variables of interest

    Parameters
    ----------

    formula: str
        A nonlinear formula including variables and parameters of the model,
        such as y ~ f(x, theta)
    data: pandas.core.frame.DataFrame
        Data frame containing the variables in the model. If not found in
        data, the variables are taken from environment(formula), typically
        the environment from which nlrob is called.
    lower: pandas.core.frame.DataFrame
        A vector of starting estimates.
    upper: array_like
        upper bound, the shape could be 1 or the the same of lower.
    psi: str
        A function (possibly by name) of the form g(x, 'tuning constant(s)',
        deriv)
        that for deriv=0 returns (x)/x and for deriv=1 returns  0 (x).
    init: str
    ctrl: object
        NlrobControl Class

    Returns
    -------
    coefficients: array_like
        Numeric vector of coefficient estimates.
    fitted_values : array_like
    residuals: array_like
        numeric vector of the residuals.
    hessian: array_like
        hessian matrix
    ctrl: object
        NlrobControl Class


    Examples
    --------
    >>> import pandas
    >>> from nlrob import *
    >>> 
    >>> formula = "density ~ Asym/(1 + np.exp(( xmid - np.log(conc) )/scal))"
    >>> lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
    ...         index=[0])
    >>> upper = np.array([1])
    >>> data = pandas.read_csv("Submodule=NLROBMM.Data=Input.csv")
    >>> Rfit_MM = nlrob_MM(formula, data, lower, lower, upper, method=method)

    See Also
    --------
    nlrob:
    nlrob_tau:
    nlrob_CM:
    nlrob_mtl:
    """
    if ctrl:
        init = ctrl.init
        psi = ctrl.psi
    c1 = ctrl.tuning_chi_scale
    c2 = ctrl.tuning_psi_M
    
    if psi == "lqq":
        c12 = c1[0] + c2[1]
        lqqMax = (c1[0] * c1[2] - 2 * c12)/( 1 - c1[2]) + c12

    rho1 = lambda t: Mchi(t, c1, psi)
    rho2 = lambda t: Mchi(t, c2, psi)

    if psi == "bisquare":
        rho_inv = lambda y: c1 * np.sqrt(1 - (1 -y) ** (0.3333333))
    elif psi == "lqq":
        rho_inv = lambda y: np.array(brentq(lambda x: rho1(x) - y, 0, lqqMax))
    elif psi == "optimal":
        rho_inv = lambda y: np.sqrt(y / 1.38) * c1 * 3
    elif psi == "hampel":
        def rho_inv(y):
            # TODO:
            C = MrhoInf(c1, psi)
            a = c1[0]
            b = c1[1]
            r = c1[2]
            if a / C > y:
                return np.sqrt(2 * C * y)
            elif (2 * b - a )/ C > y:
                return 0.5 * a + C / a * y
            else:
                return r + np.sqrt(r ** 2 - ((r - b) * (2 *C / a * y + (b - a)) \
                       - b*r ))
    else:
        raise Exception("Psi function '%s' not supported yet" % psi)

    M_scale = lambda sigma, u: np.sum(rho1(u / sigma)) / nobs - 0.5
    globals().update({"M_scale": M_scale})

    # Defining local variables of functions
    data_vars = data.keys()
    for var in data_vars:
        if re.search("\s", var.strip()):
            continue
        globals().update({var:data[var].values})
    

    par = lower.keys()
    
    y = data[formula.split("~")[0].rstrip()].values
    right_hand_term = formula.split("~")[1]
    for i, pname in enumerate(par):
        right_hand_term = right_hand_term.replace(pname, "vector[%d]" % i)
    globals().update({"right_hand_term": right_hand_term})

    if init == "lts":
        def objective_initial(vector):
            global right_hand_term
            y_hat = eval(right_hand_term)
            return np.sum(np.sort(y - y_hat)[:h] ** 2)

    elif init == "S":
        def objective_initial(vector):
            global right_hand_term
            y_hat = eval(right_hand_term)
            res = y - y_hat
            med_abs_res = np.median(np.abs(res))
            return np.array(brentq(M_scale, constant[0] *  med_abs_res, constant[1] * med_abs_res, args=(res)))
    else:
        raise Exception("Initialization 'init = \"%s\"' not supported (yet)" %init)

    def objective_M(vector, sigma):
        global right_hand_term
        y_hat = eval(right_hand_term)
        return np.sum(rho2( (y - y_hat) / sigma))


    def fminfn(p, sigma):
        global fnscale
        global parscale
        return objective_M(p * parscale, sigma) / fnscale

    def fmingr(p, sigma):
        global fnscale
        global parscale

        x = np.zeros_like(p)
        df = np.zeros_like(p)
        for i in range(p.size):
            epsused = eps = 1e-3
            tmp = p[i] + eps
            if tmp > upper[i]:
                tmp = upper[i]
                epsused = tmp - p[i]
            x[i] = tmp * parscale[i]
            s = objective_M(p, sigma)
            val1 = s / fnscale
            tmp = p[i] - eps
            if (tmp < lower[i]):
                tmp = lower[i]
                eps = p[i] - tmp
            x[i] = tmp * parscale[i]
            s = objective_M(p, sigma)
            val2 = s/ fnscale
            df[i] = (val1 - val2)/(epsused + eps)
            if df[i] == np.Inf or df[i] == -np.Inf:
                raise Exception("non-finite finite-difference value [%d]" % i+1)
            x[i] = p[i] * parscale[i]
        return df

    npar = len(par)
    nobs = y.size

    if npar > nobs:
        raise Exception("npar > nobs")
    if ctrl.fnscale:
        fnscale = ctrl.fnscale
    else:
        fnscale = np.sum((y - np.mean(y)) ** 2)
    globals().update({"fnscale": fnscale})
    if psi == "bisquare":
        constant = [1 / c1]
    elif psi == "lqq":
        constant = [1 / lqqMax]
    elif psi == "optimal":
        constant = [1 / c1 * 1 / 3]
    elif psi == "hampel":
        constant = [1 / c1[2]]
    constant.append(2  / rho_inv(2 / (nobs + 2)) if nobs % 2 else 1 / rho_inv( 1 /(nobs + 1)))
    globals().update({"constant": constant})

    if init == "lts":
        h = (nobs + npar + 1) // 2
    if npar > upper.size:
        upper = np.repeat(upper, npar)
    initial = JDEoptim(lower.values[0], upper, objective_initial, tol=tol,
                       fnscale=fnscale)
    parscale = initial.get("par")
    globals().update({"parscale": parscale})
    for var in par:
        exec("%s = lower['%s'].values" % (var, var), globals(), locals())
    res = y - eval(formula.split("~")[1])
    med_abs_res = np.median(np.abs(res))
    sigma = np.array(brentq(M_scale, constant[0] * med_abs_res, constant[1] * med_abs_res, args=(res)))
    lower = lower.values.ravel()
    bounds = Bounds(lb=lower, ub=upper)

    globals().update({"sigma": sigma})
    M = minimize(fminfn, initial.get("par"), jac=fmingr, args=sigma, method='L-BFGS-B',
                 bounds=bounds, tol=tol)
    coef = dict(zip(par, M.x))
    if M.status == 0:
        status = "converged"
    elif M.status == 1:
        status = "maximum number of iterations reached without convergence"
    else:
        status = M.message

    for var in par:
        exec("%s = coef['%s']" % (var, var), globals(), locals())
    try:
        hess = np.linalg.inv(M.hess_inv.todense())
    except:
        hess = None
    vector = M.x
    fit = eval(right_hand_term)
    dictReturn = {"formula": formula,
                  "nobs": nobs,
                  "coefficients": coef,
                  "fitted_values": fit,
                  "residuals": y - fit,
                  "ctrl": ctrl,
                  "crit": M.fun,
                  "initial": initial,
                  "Scale": sigma,
                  "status": status,
                  "hessian": hess}
    return dictReturn


def nlrob_tau(formula, data, lower, upper, tol=1e-6, psi="bisquare",
              ctrl=NlrobControl("tau"), tuning_chi_scale=None,
              tuning_chi_tau=None):
    """
    Computes a Tau-estimator for nonlinear robust (constrained) regression.
    Returns a dictionary with all the variables of interest

    Parameters
    ----------
    formula: str
        A nonlinear formula including variables and parameters of the model,
        such as y ~ f(x, theta)

    data: pandas.core.frame.DataFrame
        Data frame containing the variables in the model. If not found in
        data, the variables are taken from environment(formula), typically
        the environment from which nlrob is called.
    lower: pandas.core.frame.DataFrame
        Dataframe with the initial guesses
    upper: array_like
        upper bound, the shape could be 1 or the the same of lower.
    psi: str
        A function (possibly by name) of the form g(x, 'tuning constant(s)',
        deriv)
        that for deriv=0 returns (x)/x and for deriv=1 returns  0 (x).
    init: str
    ctrl: object
        NlrobControl Class

    Returns
    -------

    coefficients: array_like
        Numeric vector of coefficient estimates.
    fitted_values : array_like
    residuals: array_like
        numeric vector of the residuals.
    hessian: array_like
        hessian matrix
    ctrl: object
        NlrobControl Class

    
    Examples
    --------
    >>> import pandas
    >>> from nlrob import *
    >>> 
    >>> formula = "density ~ Asym/(1 + np.exp(( xmid - np.log(conc) )/scal))"
    >>> lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
    ...         index=[0])
    >>> upper = np.array([1])
    >>> data = pandas.read_csv("Submodule=NLROBTAU.Data=Input.csv")
    >>> Rfit_tau = nlrob_tau(formula, data, lower, lower, upper, method=method)

    See Also
    --------
    nlrob:
    nlrob_MM:
    nlrob_CM:
    nlrob_mtl:
    """
    if ctrl:
        psi = ctrl.psi
    
    # Defining local variables of functions
    data_names = data.keys()
    for var in data_names:
        if re.search("\s", var.strip()):
            continue
        globals().update({var:data[var].values})
    

    if ctrl.tuning_chi_scale is None:
        if psi == "bisquare":
            _chi_s = {"b": 0.2, "cc": 1.55}
        elif psi == "optimal":
            _chi_s = {"b": 0.5, "cc": 0.405}

    if ctrl.tuning_chi_tau is None:
        if psi == "bisquare":
            _chi_t = {"b": 0.46, "cc": 6.04}
        elif psi == "optimal":
            _chi_t = {"b": 0.128, "cc": 1.06}

    b1 = _chi_s.get("b")
    c1 = _chi_s.get("cc")
    b2 = _chi_t.get("b")
    c2 = _chi_t.get("cc")

    if psi == "bisquare":
        b1 = b1 / MrhoInf(c1, psi)
        b2 = b2 / MrhoInf(c1, psi)

    rho1 = lambda t: Mchi(t, c1, psi)
    rho2 = lambda t: Mchi(t, c2, psi)


    if psi == "bisquare":
        rho_inv = lambda y: c1 * np.sqrt(1 - (1 - y)**(1/3))
    elif psi == "optimal":
        rho_inv = lambda y: np.sqrt(y/1.38) * c1 * 3

    M_scale = lambda sigma, u: np.sum( rho1(u/sigma) )/nobs - b1
    tau_scale2 = lambda u, sigma: sigma ** 2 * 1 / b2 * np.sum(rho2(u / sigma))/ nobs
    par = lower.keys()
    
    y = data[formula.split("~")[0].rstrip()].values
    right_hand_term = formula.split("~")[1]
    for i, pname in enumerate(par):
        right_hand_term = right_hand_term.replace(pname, "vector[%d]" % i)

    def objective(vector):
        fit = eval(right_hand_term)
        res = y - fit
        med_abs_res = np.median(np.abs(res))
        
        sigma = np.array(brentq(M_scale, constant[0] *  med_abs_res,
            constant[1] * med_abs_res, args=(res)))
        return tau_scale2(res, sigma)

    npar = len(par)
    nobs = y.size

    if npar > nobs:
        raise Exception("npar > nobs")
    if ctrl.fnscale:
        fnscale = ctrl.fnscale
    else:
        fnscale = np.mean((y - np.mean(y)) ** 2)

    if psi == "bisquare":
        constant = [1 / c1]
    elif psi == "optimal":
        constant = [1 / c1 * 1 / 3]

    constant.append(2  / rho_inv(2 / (nobs + 2)) if nobs % 2 else 1 / rho_inv( 1 /(nobs + 1)))

    if npar > upper.size:
        upper = np.repeat(upper, npar)
    optRes = JDEoptim(lower.values[0], upper, objective, tol=tol,
                       fnscale=fnscale)
    it = optRes.get("iter")
    status = "converged" if optRes.get("convergence") == 0 else "failed to convergence in %d steps" % it
    coef = dict(zip(lower.keys(),  optRes.get("par")))
    vector = optRes.get("par")
    fit = eval(right_hand_term)
    dictReturn = {"formula": formula,
                  "nobs": nobs,
                  "coefficients": coef,
                  "fitted_values": fit,
                  "residuals": y - fit,
                  "crit": optRes.get("value"),
                  "initial": optRes,
                  "Scale": np.sqrt(optRes.get("value")),
                  "status": status,
                  "iter": it,
                  "ctrl": ctrl}
    return dictReturn


def nlrob_CM(formula, data, lower, upper, tol=1e-6, psi="bisquare",
             ctrl=NlrobControl("CM")):
    """
    Compute an MM-estimator for nonlinear robust (constrained) regression.
    Returns a dictionary with all the variables of interest
    Parameters
    ----------

    formula: str
        A nonlinear formula including variables and parameters of the model,
        such as y ~ f(x, theta)

    data: pandas.core.frame.DataFrame
        Data frame containing the variables in the model. If not found in
        data, the variables are taken from environment(formula), typically
        the environment from which nlrob is called.
    lower: pandas.core.frame.DataFrame
        Dataframe with the initial guesses
    upper: array_like
        upper bound, the shape could be 1 or the the same of lower.
    psi: str
        A function (possibly by name) of the form g(x, 'tuning constant(s)',
        deriv)
        that for deriv=0 returns (x)/x and for deriv=1 returns  0 (x).
    init: str
    ctrl: object
        NlrobControl Class

    Returns
    -------

    coefficients: array_like
        Numeric vector of coefficient estimates.
    fitted_values : array_like
    residuals: array_like
        numeric vector of the residuals.
    hessian: array_like
        hessian matrix
    ctrl: object
        NlrobControl Class
    
    Examples
    --------
    >>> import pandas
    >>> from nlrob import *
    >>> 
    >>> formula = "density ~ Asym/(1 + np.exp(( xmid - np.log(conc) )/scal))"
    >>> lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
    ...         index=[0])
    >>> upper = np.array([1])
    >>> data = pandas.read_csv("Submodule=NLROBCM.Data=Input.csv")
    >>> Rfit_CM = nlrob_cm(formula, data, lower, lower, upper, method=method)

    See Also
    --------
    nlrob:
    nlrob_MM:
    nlrob_tau:
    nlrob_mtl:

    """
    if ctrl:
        psi = ctrl.psi
    if psi == "bisquare":
        t_chi = {"b": 0.5, "cc":1, "c":4.835}
    b = t_chi.get("b")
    c = t_chi.get("c")
    cc = t_chi.get("cc")

    rho = lambda t: Mchi(t, cc, psi)
    M_scale = lambda sigma, u: np.sum(rho(u / sigma)) / nobs - b

    # Defining local variables of functions
    data_names = data.keys()
    for var in data_names:
        if re.search("\s", var.strip()):
            continue
        globals().update({var:data[var].values})

    pnames = [name for name in lower.keys()  if name in lower.keys() or name == "sigma"]

    y = data[formula.split("~")[0].rstrip()].values

    right_hand_term = formula.split("~")[1]
    for i, pname in enumerate(pnames):
        right_hand_term = right_hand_term.replace(pname, "vector[%d]" % i)

    if "sigma" in pnames:
        if "sigma" in formula.split("~")[1] or "sigma" in data.keys():
            raise Exception("As \"sigma\" is in 'pnames', do not use it as variable or parameter name in 'formula'")
        def objective(vector):
            fit = eval(right_hand_term)
            sigma = vector[-1]
            return c * np.sum(rho( (y - fit)/sigma ))/nobs + np.log(sigma)
        def con(vector):
            fit = eval(right_hand_term)
            return M_scale(vector[-1], y - fit)
    else:
        def objective(vector):
            fit = eval(right_hand_term)
            res = y - fit
            sigma = np.median(np.abs(res - np.median(res)))
            return c * np.sum(rho(res / sigma)) / nobs + np.log(sigma)
        con = None

    npar = len(pnames)
    nobs = y.size

    if npar > nobs:
        raise Exception("npar > nobs")
    if ctrl.fnscale:
        fnscale = ctrl.fnscale
    else:
        fnscale = np.mean((y - np.mean(y)) ** 2)
    if npar > upper.size:
        upper = np.repeat(upper, lower.size)

    optRes = JDEoptim(lower.values[0], upper, objective, tol=tol,
                    fnscale=fnscale, constr=con)
    it = optRes.get("iter")
    status = "converged" if optRes.get("convergence") == 0 else "failed to convergence in %d steps" % it
    coef = dict(zip(lower.keys(),  optRes.get("par")))
    vector = optRes.get("par")
    fit = eval(right_hand_term)
    dictReturn = {"formula": formula,
                  "nobs": nobs,
                  "coefficients": coef,
                  "fitted_values": fit,
                  "residuals": y - fit,
                  "crit": optRes.get("value"),
                  "status": status,
                  "iter": it,
                  "ctrl": ctrl}
    return dictReturn



def nlrob_mtl(formula, data, lower, upper, tol=1e-6, psi="bisquare",
             ctrl=NlrobControl("mtl")):
    """
    Compute a mtl-estimator for nonlinear robust (constrained) regression

    Parameters
    ----------

    formula: str
        A nonlinear formula including variables and parameters of the model,
        such as y ~ f(x, theta)

    data: pandas.core.frame.DataFrame
        Data frame containing the variables in the model. If not found in
        data, the variables are taken from environment(formula), typically
        the environment from which nlrob is called.
    lower: pandas.core.frame.DataFrame
        Dataframe with the initial guesses
    upper: array_like
        upper bound, the shape could be 1 or the the same of lower.
    psi: str
        A function (possibly by name) of the form g(x, 'tuning constant(s)',
        deriv)
        that for deriv=0 returns (x)/x and for deriv=1 returns  0 (x).
    init: str
    ctrl: object
        NlrobControl Class

    Returns
    -------

    coefficients: array_like
        Numeric vector of coefficient estimates.
    fitted_values : array_like
    residuals: array_like
        numeric vector of the residuals.
    hessian: array_like
        hessian matrix
    ctrl: object
        NlrobControl Class


    Examples
    --------
    >>> import pandas
    >>> from nlrob import *
    >>> 
    >>> formula = "density ~ Asym/(1 + np.exp(( xmid - np.log(conc) )/scal))"
    >>> lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
    ...         index=[0])
    >>> upper = np.array([1])
    >>> data = pandas.read_csv("Submodule=NLROBMTL.Data=Input.csv")
    >>> Rfit_mtl = nlrob_mtl(formula, data, lower, lower, upper, method=method)

    See Also
    --------
    nlrob:
    nlrob_MM:
    nlrob_tau:
    nlrob_cm:
    """
    if ctrl:
        cutoff = ctrl.cutoff
    def trim(t):
        t = np.sort(t)
        i = np.where(t >= cutoff)[0]
        partial_h = np.min( (i - 1)/(2 * np.random.normal(t[i]) - 1))
        partial_h = np.max(np.floor(partial_h))
        h = np.max([hlow, partial_h]) if i.size else nobs
        return {"h": h, "t": t}


    rho = lambda t: Mchi(t, cc, psi)
    M_scale = lambda sigma, u: np.sum(rho(u / sigma)) / nobs - b

    # Defining local variables of functions
    data_names = data.keys()
    for var in data_names:
        if re.search("\s", var.strip()):
            continue
        globals().update({var:data[var].values})

    pnames = [name for name in lower.keys()  if name in lower.keys() or name == "sigma"]

    y = data[formula.split("~")[0].rstrip()].values
    right_hand_term = formula.split("~")[1]
    for i, pname in enumerate(pnames):
        right_hand_term = right_hand_term.replace(pname, "vector[%d]" % i)
    constant = np.log(2 * np.pi)
    if "sigma" in pnames:
        if "sigma" in formula.split("~")[1] or "sigma" in data.keys():
            raise Exception("As \"sigma\" is in 'pnames', do not use it as variable or parameter name in 'formula'")
        def objective(vector):
            fit = eval(right_hand_term)
            sigma = vector[-1]
            tp = trim(np.abs( (y - fit) / sigma))
            h = tp.get("h")

            return h * (constant + 2 * np.log(sigma)) + np.sum(tp.get("t")[:h] ** 2)
    else:
        def objective(vector):
            fit = eval(right_hand_term)
            res = y - fit
            sigma = np.median(np.abs(res - np.median(res)))
            tp = trim(np.abs(res / sigma))
            h = int(tp.get("h"))

            return h * (constant + 2 * np.log(sigma)) + np.sum(tp.get("t")[:h] ** 2)

    npar = len(pnames)
    nobs = y.size

    if npar > nobs:
        raise Exception("npar > nobs")
    if ctrl.fnscale:
        fnscale = ctrl.fnscale
    else:
        fnscale = np.mean((y - np.mean(y)) ** 2)
    if npar > upper.size:
        upper = np.repeat(upper, lower.size)
    hlow = (nobs + npar + 1) // 2

    optRes = JDEoptim(lower.values[0], upper, objective, tol=tol,
                    fnscale=fnscale)
    it = optRes.get("iter")
    status = "converged" if optRes.get("convergence") == 0 else "failed to convergence in %d steps" % it
    coef = dict(zip(lower.keys(),  optRes.get("par")))
    vector = optRes.get("par")
    fit = eval(right_hand_term)
    res = y - fit
    quan = trim( res/(coef["sigma"] if ("sigma" in pnames)  else np.median(np.abs(res - np.median(res))))).get("h") 
    dictReturn = {"formula": formula,
                  "nobs": nobs,
                  "coefficients": coef,
                  "fitted_values": fit,
                  "residuals": res,
                  "crit": optRes.get("value"),
                  "quan": quan,
                  "status": status,
                  "iter": it,
                  "ctrl": ctrl}
    return dictReturn


def nls(formula, data, start, algorithm="lm",
         weights=None, lower=np.array([-np.Inf]), upper=np.array([np.Inf])):

    """
    Determine the nonlinear (weighted) least-squares estimates of the parameters of a nonlinear model.

    Usage

    nls(formula, data, start, control, algorithm,
        trace, subset, weights, na.action, model,
        lower, upper, ...)
    
    Parameters
    ----------

    formula: str
        a nonlinear model formula including variables and parameters. Will be coerced to a formula if necessary.

    data: pandas.core.frame.DataFrame
        Data frame in which to evaluate the variables in formula and weights. Can also be a list or an environment, but not a matrix.

    start: pandas.core.frame.DataFrame
        A vector of starting estimates.

    algorithm: str
        Character string specifying the algorithm to use. The default algorithm is a "lm" algorithm. Other possible values are 'trf', 'dogbox'

    lower: scalar, array_like
        Lower bounds on Parameters. An array with the length equal
        to the number of parameters, or a scalar (in which case the bound is
        taken to be the same for all parameters.)

    upper: scalar, array_like
        Upper bounds on Parameters. An array with the length equal
        to the number of parameters, or a scalar (in which case the bound is
        taken to be the same for all parameters.)

    Returns
    -------

    coefficients: array_like
        Numeric vector of coefficient estimates.
    residuals: array_like
        numeric vector of the residuals.
    cov: array_like
        covariance matrix

    Examples
    --------
    >>> import pandas
    >>> from nlrob import *
    >>> 
    >>> formula = "density ~ Asym/(1 + np.exp(( xmid - np.log(conc) )/scal))"
    >>> lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
    ...         index=[0])
    >>> upper = np.array([1])
    >>> data = pandas.read_csv("Submodule=NLS.Data=Input.csv")
    >>> Rfit_nls = nls(formula, data, lower)

    See Also
    --------
    nlrob:
    nlrob_MM:
    nlrob_tau:
    nlrob_cm:

    """

    for var in data.keys():
        if re.search("\s", var.strip()):
            continue
        globals().update({var:data[var].values})

    y = data[formula.split("~")[0].rstrip()].values
    right_hand_term = formula.split("~")[1]
    for pname in data.keys():
        i = 0
        if pname in right_hand_term:
            break
    pnames = start.keys()
    p0 = start.values[0]
    x = pname
    if lower.any():
        
        if lower.size != p0.size:
            lower = np.repeat(lower, p0.size)
    if upper.any():
        if upper.size != p0.size:
            upper = np.repeat(upper, p0.size)
    bounds = (lower, upper)
    def get_func():
        env = {"np": np}
        code =  "def func(%s, %s):\n" % (x, ", ".join(pnames))
        code += "    return %s\n" % right_hand_term
        exec(code, env)
        return env.get("func")
    other_params = "bounds=bounds, method='%s'" % algorithm
    func = get_func()
    par, cov = eval("curve_fit(func, %s, y, p0, %s)" % (x, other_params),
                    globals(), locals())
    popt_str = ["par[%d]" % i for i, p in enumerate(par)]
    fit = eval("func(%s, %s)" % (x, ", ".join(popt_str)))
    res = y - fit
    dictReturn = {"coefficients": dict(zip(start.keys(), par)),
                  "residuals": res,
                  "cov":cov}
    return dictReturn



