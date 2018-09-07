################################################################################
# lmrob
#
# This is a translation from lmrob function of lmrob packages in R
#
# Includes the following functions:
#
#     - lmrob
#     - lmrob_control
#     - lmrob_fit
#     - lmrob_S
#     - lmrob_D_fit
#     - lmrob_M_fit
#     - lmrob_M_S
#     - lmrob_lar
#
#
################################################################################

from .utils import *
from .utils import _psi_conv_cc, _psi2ipsi, _regularize_Mpsi, _convSS, _Mpsi, _lmrob_hat, _vcov_avar1, _vcov_w

class LmrobControl():
    """
    Tuning parameters for lmrob, the MM-type regression estimator and the associated S-, M- and D-
    estimators. Using setting="KS2011" sets the defaults as suggested by Koller and Stahel (2011)
    and analogously for "KS2014".
    The .M*.default functions and .M*.defaults lists contain default tuning parameters for all
    the predefined ψ functions, see also Mpsi, etc.

    Parameters
    ----------

    setting: str
        a string specifying alternative default values. Leave empty for the defaults or use
        "KS2011" or "KS2014" for the defaults suggested by Koller and Stahel (2011,
        2017). See Notes.

    seed: int
        None or an integer vector compatible with _Random_seed: the seed to be used
        for random re-sampling used in obtaining candidates for the initial S-estimator.
        The current value of .Random.seed will be preserved if seed is set, i.e. non-
        None; otherwise, as by default, .Random.seed_will be used and modified as
        usual from calls to runif() etc.
    nResample: int
        number of re-sampling candidates to be used to find the initial S-estimator.
        Currently defaults to 500 which works well in most situations (see references).
    tuning_chi: int
        tuning constant vector for the S-estimator. If None, as by default, sensible de-
        faults are set (depending on psi) to yield a 50% breakdown estimator. See
        Notes.
    bb: int
        expected value under the normal model of the "chi" (rather ρ(rho)) function
        with tuning constant equal to tuning.chi. This is used to compute the S-
        estimator.
    tuning_psi: int
        tuning constant vector for the redescending M-estimator. If None, as by default,
        this is set (depending on psi) to yield an estimator with asymptotic efficiency
        of 95% for normal errors. See Notes.
    max_it: int
        integer specifying the maximum number of IRWLS iterations.
        groups (for the fast-S algorithm): Number of random subsets to use when the data set
        is large.
    n_group: int
        Size of each of the groups above. Note that this must
        be at least p.
    k_fast_s: int
        (for the fast-S algorithm): Number of local improvement steps ("I-steps") for
        each re-sampling candidate.
    k_m_s: int
        (for the M-S algorithm): specifies after how many unsucessful refinement steps
        the algorithm stops.
    best_r_s: int
        (for the fast-S algorithm): Number of of best candidates to be iterated further
        (i.e., "refined"); is denoted t in Salibian-Barrera & Yohai(2006).
    k_max: int
        (for the fast-S algorithm): maximal number of refinement steps for the "fully"
        iterated best candidates.
    maxit_scale: int
        integer specifying the maximum number of C level find_scale() iterations.
    refine_tol: float
        (for the fast-S algorithm): relative convergence tolerance for the fully iterated
        best candidates.
    rel_tol: float
        (for the RWLS iterations of the MM algorithm): relative convergence tolerance
        for the parameter vector
    scale_tol: float
        (for the scale estimation iterations of the S algorithm): relative convergence
        tolerance for the scale σ(.).
    solve_tol: float
        (for the S algorithm): relative tolerance for inversion. Hence, this corresponds
        to solve_default()’s tol.
    trace_lev: int
        integer indicating if the progress of the MM-algorithm should be traced (increas-
        ingly); default trace.lev = 0 does no tracing.
    mts: int
        maximum number of samples to try in subsampling algorithm.
    subsampling: str
        type of subsampling to be used, a string: "simple" for simple subsampling
        (default prior to version 0.9), "nonsingular" for nonsingular subsampling
    compute_rd: bool
        boolean indicating if robust distances (based on the MCD robust covariance es-
        timator covMcd) are to be computed for the robust diagnostic plots. This may
        take some time to finish, particularly for large data sets, and can lead to singu-
        larity problems when there are factor explanatory variables (with many levels,
        or levels with "few" observations). Hence, is False by default.
        method string specifying the estimator-chain. MM is interpreted as SM. See Details of
        lmrob for a description of the possible values.
    psi: str
        string specifying the type ψ-function used. See Details of lmrob. Defaults to
        "bisquare" for S and MM-estimates, otherwise "lqq".
    numpoints: int
        number of points used in Gauss quadrature.
    cov: function
        function name to be used to calculate covariance matrix
        estimate. The default is if(method in ['SM', 'MM']) "_vcov_avar1" else "_vcov_w".
    split_type: str
        determines how categorical and continuous variables are split. See splitFrame.
        fasts_large_n minimum number of observations required to switch from ordinary "fast S"
        algorithm to an efficient "large n" strategy.
    eps_outlier: float
        limit on the robustness weight below which an observation is considered to be an
        outlier.
    eps_x: float
        limit on the absolute value of the elements of the design matrix below which
        an element is considered zero.
    """
    def __init__(self,
            setting=None,
            seed=None,
            nResample=500,
            tuning_chi=None,
            bb=0.5,
            tuning_psi=None,
            max_it=50,
            groups=5,
            n_group=400,
            k_fast_s=1,
            best_r_s=2,
            k_max=200,
            maxit_scale=200,
            k_m_s=20,
            refine_tol=1e-7,
            rel_tol=1e-7,
            scale_tol=1e-10,
            solve_tol=1e-7,
            trace_lev=0,
            mts=1000,
            subsampling="nonsingular",
            compute_rd=False,
            method='MM',
            psi=None,
            numpoints=10,
            cov=None,
            split_type="f",
            fast_s_large_n=2000,
            eps_outlier=None,
            eps_x=None,
            compute_outlier_stats=None,
            warn_limit_reject=0.5,
            warn_limit_meanrw=0.5
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

        if not setting is None:
            if setting in ['KS2011', 'KS2014']:
                self.method = "SMDM" if method is None else method
                psi = 'lqq' if psi is None else _regularize_Mpsi(psi)
                self.max_it = 500 if max_it is None else max_it
                self.k_max = 2000 if k_max is None else k_max
                self.cov = _vcov_w if cov is None else cov
                if setting == 'KS2014':
                    self.best_r_s = 20 if best_r_s is None else None
                    self.k_fast_s = 2 if k_fast_s is None else None
                    self.nResample = 1000 if nResample is None else None
        else:
            if not psi is None:
                psi = psi
            elif psi is None and "D" in method:#re.search(r"D$", method):
                psi = 'lqq'
            elif not "D" in method:
                psi = 'bisquare'
            if cov is None:
                cov = _vcov_avar1 if method in ['SM', 'MM'] else _vcov_w
        if not psi is None:
            psi = _regularize_Mpsi(psi)
        
        # In ggw, lqq:  if tuning.{psi|chi}  are non-standard, calculate coefficients:
        compute_const = psi in ['ggw', 'lqq']
        if tuning_chi is None:
            tuning_chi = _Mchi_tuning_defaults.get(psi, None)
        else:
            if compute_const:
                tuning_chi = _psi_const(tuning_chi, psi)

        if tuning_psi is None:
            tuning_psi = _Mpsi_tuning_defaults.get(psi, None)
        else:
            if compute_const:
                tuning_psi = _psi_const(tuning_psi, psi)

        tuning_psi = np.array([tuning_psi], dtype="float") if np.isscalar(tuning_psi) else tuning_psi
        tuning_chi = np.array([tuning_chi], dtype="float") if np.isscalar(tuning_chi) else tuning_chi

        self.method = method
        self.max_it = max_it
        self.k_max = k_max
        self.cov = cov
        self.best_r_s = best_r_s
        self.k_fast_s = k_fast_s
        self.nResample = nResample
        self.seed = seed
        self.tuning_chi = tuning_chi
        self.bb = bb
        self.tuning_psi = tuning_psi
        self.groups = groups
        self.n_group = n_group
        self.maxit_scale = maxit_scale
        self.k_m_s = k_m_s
        self.refine_tol = refine_tol
        self.rel_tol = rel_tol
        self.scale_tol = scale_tol
        self.solve_tol = solve_tol
        self.trace_lev = trace_lev
        self.mts = mts
        self.subsampling = subsampling
        self.compute_rd = compute_rd
        self.psi = psi
        self.numpoints = numpoints
        self.split_type = split_type
        self.fast_s_large_n = fast_s_large_n
        self.eps_outlier = eps_outlier
        self.eps_x = eps_x
        self.compute_outlier_stats = compute_outlier_stats
        self.warn_limit_reject = warn_limit_reject
        self.warn_limit_meanrw = warn_limit_meanrw

    def copy(self):
        """
        Rerturns a copy of the object
        """
        return copy.copy(self)

    def __str__(self):
        string =  "setting = {:}\n".format(self.setting)
        string = string + "seed = {:}\n".format(self.seed)
        string = string + "nResample = {:}\n".format(self.nResample)
        string = string + "tuning_chi = {:}\n".format(self.tuning_chi)
        string = string + "bb = {:}\n".format(self.bb)
        string = string + "tuning_psi = {:}\n".format(self.tuning_psi)
        string = string + "max_it = {:}\n".format(self.max_it)
        string = string + "groups = {:}\n".format(self.groups)
        string = string + "n_group = {:}\n".format(self.n_group)
        string = string + "k_fast_s = {:}\n".format(self.k_fast_s)
        string = string + "best_r_s = {:}\n".format(self.best_r_s)
        string = string + "k_max = {:}\n".format(self.k_max)
        string = string + "maxit_scale = {:}\n".format(self.maxit_scale)
        string = string + "k_m_s = {:}\n".format(self.k_m_s)
        string = string + "refine_tol = {:}\n".format(self.refine_tol)
        string = string + "rel_tol = {:}\n".format(self.rel_tol)
        string = string + "scale_tol = {:}\n".format(self.scale_tol)
        string = string + "solve_tol = {:}\n".format(self.solve_tol)
        string = string + "trace_lev = {:}\n".format(self.trace_lev)
        string = string + "mts = {:}\n".format(self.mts)
        string = string + "subsampling = {:}\n".format(self.subsampling)
        string = string + "compute_rd = {:}\n".format(self.compute_rd)
        string = string + "method = {:}\n".format(self.method)
        string = string + "psi = {:}\n".format(self.psi)
        string = string + "numpoints = {:}\n".format(self.numpoints)
        string = string + "cov = {:}\n".format(self.cov)
        string = string + "split_type = {:}\n".format(self.split_type)
        string = string + "fast_s_large_n = {:}\n".format(self.fast_s_large_n)
        string = string + "eps_outlier = {:}\n".format(self.eps_outlier)
        string = string + "eps_x = {:}\n".format(self.eps_x)
        string = string + "compute_outlier_stats = {:}\n".format(self.compute_outlier_stats)
        string = string + "warn_limit_reject = {:}\n".format(self.warn_limit_reject)
        string = string + "warn_limit_meanrw = {:}\n".format(self.warn_limit_meanrw)
        return string




def lmrob(formula="", data={}, subset=None, weights=np.array([]), na_action='drop',
          method="", model=True, return_x=True, return_y=True,
          singular_ok=True, contrasts={}, offset=np.array([]), control=None, init={},
          **kwargs):
    """
    
    
    Computes fast MM-type estimators for linear (regression) models  and
    returns a dictionary object.
    

    Parameters
    ----------

    formula : str
        A symbolic description of the model to be fit.
    data : dict
        An optional data frame, dict containing the variables in the model.
    subset: array_like
        An optional vector specifying a subset of observations to be used in the fitting
        process.
    weights: array_like
        An optional vector of weights to be used in the fitting process (in addition to the
        robustness weights computed in the fitting process).
    na_action: str
        A function which indicates what should happen when the data contain nan. The
        default is set by the na_action setting of options, and is na_fail if that is
        unset. The "factory-fresh"   default is na_omit. Another possible value is None,
        no action. Value na_exclude can be useful.
    method: str
        string specifying the estimator-chain. "MM" is interpreted as SM.
        See Notes, notably the currently recommended setting = "KS2014".
    model, return_x, return_y: bool.
        If True the corresponding components of  the fit (the model frame, the model
        matrix, the response) are returned.
    singular_ok: bool
        If False (the default in S but not in R) a singular fit is an error.
    contrasts: dict
        An optional list. See the contrasts.arg of model.matrix.default.
    offset: array_like
        This can be used to specify an a priori known component to be included in the
        linear predictor during fitting. An offset term can be included in the formula
        instead or as well, and if both are specified their sum is used.
    control: object
        Instances of class LmrobControl.
    init: dict
        An optional argument to specify or supply the initial estimate. *See Notes*.
    kwargs: 
        additional arguments can be used to specify control parameters directly instead
        of (but not in addition to!) via control.

    Returns
    -------

    coefficients: array_like
        The estimate of the coefficient vector
    scale:
        The scale as used in the M estimator.
    residuals: array_like
        Residuals associated with the estimator.
    converged: bool
        True if the IRWLS iterations have converged.
    iter: int
        Number of IRWLS iterations
    rweights: array_like
        The "robustness weights" ψ(r i /S)/(r i /S).
    fitted_values: array_like
        Fitted values associated with the estimator.
    init: dict
        A similar dict that contains the results of intermediate estimates.
    rank: int
        the numeric rank of the fitted linear model.

    Examples
    --------
    >>> from lmrob import *
    >>> from rpy2.robjects.packages import importr
    >>> from rpy2.robjects import pandas2ri, r
    >>> import rpy2.robjects.numpy2ri as rpyn
    >>>
    >>> stackloss = r("stackloss")
    >>> data = {"AirFlow" : rpyn.ri2py(stackloss.rx2("Air.Flow")),
    ...         "WaterTemp" : rpyn.ri2py(stackloss.rx2("Water.Temp")),
    ...         "AcidConc" : rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
    ...         "stack_loss" : rpyn.ri2py(stackloss.rx2("stack.loss"))
    ...      }
    >>> formula = 'stack_loss ~ AirFlow + WaterTemp + AcidConc'
    >>> # S Method
    >>> m0 = lmrob(formula, data=data, method="S")

    Notes
    -----

    **Overview**:
        This function computes an MM-type regression estimator as described in
        Yohai (1987) and Koller and Stahel (2011). By default it uses a bi-square
        redescending score function, and it returns a highly robust and highly
        efficient estimator (with 50% breakdown point and 95% asymptotic
        efficiency for normal errors). The computation is carried out by a call
        to lmrob.fit(). The argument setting of lmrob.control is provided to set
        alternative defaults as suggested in Koller and Stahel (2011)
        (setting="KS2011"; now do use its extension setting="KS2014").
        For further details, see lmrob_control.

    **Initial Estimator** ``init``:
        The initial estimator may be specified using the argument init. This
        can either be a string, a function or a list. A string can be used to
        specify built in internal estimators (currently S and M-S, see See also
        below). A function taking arguments x, y, control, mf
        (where mf stands for model.frame) and returning a list containing at least
        the initial coefficients as coefficients and the initial scale estimate
        scale. Or a list giving the initial coefficients and scale as coefficients
        and scale. See also Examples. Note that if the init argument is a function
        or list, the method argument must not contain the initial estimator, e.g.,
        use MDM instead of SMDM.
        The default, equivalent to init = "S", uses as initial estimator an
        S-estimator (Rousseeuw and Yohai, 1984) which is computed using the
        Fast-S algorithm of Salibian-Barrera and Yohai (2006), calling lmrob.S().
        That function, since March 2012, by default uses nonsingular subsampling
        which makes the Fast-S algorithm feasible for categorical data as well,
        see Koller (2012). Note that convergence problems may still show up as
        warnings, e.g.,S refinements did not converge (to refine.tol=1e-07) in
        200 (= k.max) steps and often can simply be remedied by increasing
        (i.e. weakening) refine.tol or increasing the allowed number of iterations
        k.max, see lmrob.control.


    See Also
    --------
    LmrobControl:
    lmrob_S:
    lmrob_M_S:
    lmrob_fit:

    """
    if control == None:
        control = LmrobControl() if method is None else LmrobControl(method=method)
    #TODO: elif :

    #TODO: store the main call

    # Get model frame
    mf_dict = model_frame(formula, data, subset, weights, na_action, offset)

    mt = mf_dict.get("terms")
    y = mf_dict.get("response")
    w = mf_dict.get("weights")
    offset = mf_dict.get("offset")

    if len(w) == 0 and w.dtype.name != 'float64':
        raise Exception("'weights' must be a numeric vector")
    if len(offset) != 0 and len(offset) != y.shape[0]:
        raise Exception("number of offsets is %d, should equal %d (number of observations)" %(len(offset), y.shape[0]))

    if not mt:
        x = None
        singular_fit = False
        z = {
            "coefficients" : np.zeros((0, y.shape[1])) * np.nan if not np.isscalar(y) else np.array([]),
            "residuals" : y,
            "scale" : np.nan,
            "fitted_values" : np.zeros_like(y),
            "cov" : np.zeros((0, 0)) * np.nan,
            "weights" : w,
            "rank" : 0,
            "df_residuals" : np.sum(w[w != 0]) if not w is None else y.shape[0],
            "converged" : True,
            "iter" : 0
        }

        if not offset is None:
            z["fitted_values"] = offset
            z["residuals"] = y - offset
            z["offset"] = offset
    else:
        #mm_dict = model_matrix(mt, mf_dict)
        x = mf_dict.get("model_matrix")

        assign = mt.get("assign")
        p = x.shape[1]

        if offset.size > 0:
            y -= offset

        if w.size > 0:
            ny = y.shape[1]
            n = x.shape[0]
            if y.shape[1] != n or len(w) != n:
                Exception("incompatible dimensions")

            zero_weights = np.any(w == 0)
            if zero_weights:
                save_r  = copy.deepcopy(y)
                save_w  = w
                save_f  = y
                ok = w != 0
                nok = ok == False
                w = w[ok.ravel()]
                x0 = x[nok.ravel(), :]
                x = x[ok.ravel(), :]
                n = x.shape[0]
                y0 = y[nok.ravel(), :]
                y = y[ok.ravel(), :]
                ## add this information to model.frame as well
                ## need it in outlierStats()
                ## ?? could also add this to na.action, then
                ##    naresid() would pad these as well.
                mf_dict.update({"zero_weights" : [i for i, nokInf in enumerate(nok.ravel()) if nokInf]})

            wts = np.sqrt(w)
            save_y = y
            x = wts * x
            y = wts * y

        ## check for singular fit
        z0 = lm_fit(x, y, tol_=control.solve_tol)
        piv = z0.get("pivot")
        rankQR = z0.get("rank")
        rankQR = rankQR if np.isscalar(rankQR) else np.asscalar(rankQR)
        singular_fit = (rankQR < p)

        if rankQR > 0:
            if singular_fit:
                if not singular_ok:
                    raise Exception("singular fit encountered")
                pivot = piv
                p1 = pivot[:rankQR]
                p2 = pivot[rankQR:]
                #TODO: dn = dimnames(x)
                x = x[:, p1]
                mt["assign"] = assign[p1] ## needed for splitFrame to work

            if (inspect.isfunction(control.eps_x)):
                control.eps_x = lambda x: np.max(np.abs(x))#control.eps.x(max(abs(x)))

            ini = init
            if init:
                if init == "M-S": 
                    init = lmrob_M_S(x, y, control, mf_dict)
                elif init == "S":
                    init = lmrob_S(x, y, control)
                else:
                    raise Exception("init must be 'S', 'M-S', function or list")

                if ini == "M-S":
                    ini = init["control"].method

                ## modify (default) control$method, possibly dropping first letter:
                if control.method == "MM" or control.method[0] == "S":
                    control.method = control.method[1:2]
                ## check for control.cov argument
                if ini != "S" and control.cov == "_vcov_avar1":
                    control.cov = "_vcov_w"

            z = lmrob_fit(x, y, control, init=init) #-> ./lmrob.MM.R
            
            ##   ---------
            if ini.__class__ == str and not ini == control.method[:len(ini)]:
                control.method = ini + control.method
            if singular_fit:
                coefficients = np.zeros((p, 1))
                coefficients[p2, :] = np.nan
                coefficients[p1, :] = z.get("coefficients")
                z["coefficients"] = coefficients
                ## Update QR decomposition (z$qr)
                ## pad qr and qraux with zeroes (columns that were pivoted to the right in z0)
                d_p = p - rankQR
                n = y.shape[0]
                # TODO: z["qr"] = {"qr":f(z$qr), "qraux":[], "pivot":piv}

        else: ## rank 0
            z = {
                "coefficients" : np.zeros([p, y.shape[1]]) * np.nan if y.shape[0] > 1 and y.shape[1] > 1 else np.zeros(p) * np.nan,
                "residuals" : y,
                "scale" : np.nan,
                "fitted_values" : 0.0 * y,
                "cov" : np.zeros([0, 0]) * np.nan,
                "rweights" : np.zeros(y.shape[0]) * np.nan,
                "weights" : w,
                "rank" : 0,
                "df_residuals" : y.shape[0],
                "converged" : True,
                "iter" : 0,
                "control" : control
            }
            if not offset is None:
                z["residuals"] = y - offset

        if w.size > 0:
            z["residuals"] = z["residuals"] / wts
            z["fitted_values"] = save_y - z["residuals"]
            z["weights"] = w
            if zero_weights:
                coef = copy.deepcopy(z["coefficients"])
                coef[np.isnan(coef)] = 0.0
                f0 = np.dot(x0, coef)
                ##      above  ok := (w != 0);  nok := (w == 0)
                save_r[ok.ravel(), :] = z["residuals"]
                save_r[nok.ravel(), :] = y0 - f0
                save_f[ok.ravel(), :] = z["fitted_values"]
                save_f[nok.ravel(), :] = f0
                z["residuals"] = save_r
                z["fitted_values"] = save_f
                z["weights"] = save_w
                rw = z["rweights"]
                z["rweights"] = np.zeros((len(save_w), 1))
                z["rweights"][ok.ravel(), :] = rw

    if offset.size > 0:
        z["fitted_values"] += offset

    z["na_action"] = na_action
    z["offset"] = offset
    z["terms"] = mt
    z["assign"] = assign
    if control.compute_rd and x.size !=0:
        pass
    if model:
        z["model"] = mf_dict
    if return_x:
        if singular_fit or (w.size != 0 and zero_weights):
            z["x"] = mf_dict["model_matrix"]
        else:
            z["x"] = x
    
    if return_y:
        if singular_fit or (w.size != 0 and zero_weights):
            z["y"] = mf_dict["response"]
        else:
            z["y"] = y
    return z


def lm_fit(x, y, tol_=1e-07):
    """
    lm fit function

    
    Parameters
    ----------

    x: array_like
        array with the independets variables must be an fortran array
    y:array_like
        array with the dependent variable
    """
    x_copy = np.asfortranarray(copy.deepcopy(x))
    y_copy = np.asfortranarray(copy.deepcopy(y))
    n = np.array([x.shape[0]], dtype="int32")
    p = np.array([x.shape[1]], dtype="int32")
    ny = np.array([y.shape[1]], dtype="int32")
    rank = np.array([0],  dtype="int32")
    coefficients = np.zeros((p[0], 1), dtype="float")
    residuals = np.zeros((n[0],1), dtype="float")
    effects = np.zeros((n[0],1), dtype="float")
    pivot = np.arange(p[0], dtype="int32")
    qraux = np.zeros(p[0], dtype="float")
    work = np.zeros(2 * p[0], dtype="float")
    tol = np.array([tol_], dtype="float")

    x_p = x_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    n_p = n.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    p_p = p.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    y_p = y_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ny_p = ny.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    tol_p = tol.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    coefficients_p = coefficients.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    residuals_p = residuals.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    effects_p = effects.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    rank_p = rank.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    pivot_p = pivot.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    qraux_p = qraux.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    work_p = work.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    dqrls = lmroblib.dqrls_
    dqrls.argtypes = [ ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_int),
                       ctypes.POINTER(ctypes.c_int),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_int),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_int),
                       ctypes.POINTER(ctypes.c_int),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_double)]
    dqrls(
        x_p,
        n_p,
        p_p,
        y_p,
        ny_p,
        tol_p,
        coefficients_p,
        residuals_p,
        effects_p,
        rank_p,
        pivot_p,
        qraux_p,
        work_p
    )
    pivoted = 0
    resultsDict = {
        "coefficients": coefficients,
        "residuals": residuals,
        "rank": rank,
        "effects": effects,
        "pivot": pivot,
        "qraux": qraux,
        "pivoted": pivoted
    }

    return resultsDict

def lmrob_S(x_, y_, control, only_scale=False):
    """
    Computes an S-estimator for linear regression, using the "fast S" algorithm and
    returns a dictionary object.
    
    Parameters
    ----------

        x_: array_like
            Design matrix (n × p)

        y_: array_like
            Numeric vector of responses (or residuals for only_scale=True).

        control: object
            Instances of class LmrobControl.
        
        only_scale: bool
            Boolean indicating if only the scale of ``y`` should be computed. In this case, ``y``
            will typically contain residuals.

    Examples
    --------    
        >>> from lmrob import *
        >>> from rpy2.robjects.packages import importr
        >>> from rpy2.robjects import pandas2ri, r
        >>> import rpy2.robjects.numpy2ri as rpyn
        >>>
        >>> stackloss = r("stackloss")
        >>> data = {"AirFlow" : rpyn.ri2py(stackloss.rx2("Air.Flow")),
        ...         "WaterTemp" : rpyn.ri2py(stackloss.rx2("Water.Temp")),
        ...         "AcidConc" : rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
        ...         "stack_loss" : rpyn.ri2py(stackloss.rx2("stack.loss"))
        ...      }
        >>> formula = 'stack_loss ~ AirFlow + WaterTemp + AcidConc'
        >>> Y, X = patsy.dmatrices(formula, data, NA_action="drop", return_type='matrix')
        >>> control = LmrobControl()
        >>> m0 = lmrob_S(X, Y, control)
        >>> coefficients = m0.get("coefficients")
        >>> residuals = m0.get("residuals")
        >>> scale = m0.get("scale")
        >>> converged = m0.get("converged")
        >>> k_iter = m0.get("k_iter")
        >>> fitted_values = m0.get("fitted_values")
        >>> rweights = m0.get("rweights")

    Notes
    -----

        This function is used by ``lmrob_fit`` and typically not to be used on its own (because an S-estimator
        has too low efficiency ‘on its own’).
        By default, the subsampling algorithm uses a customized LU decomposition which ensures a non
        singular subsample (if this is at all possible). This makes the Fast-S algorithm also feasible for
        categorical and mixed continuous-categorical data.
        One can revert to the old subsampling scheme by setting the parameter subsampling in control
        to "simple".

    Returns
    -------

        coefficients: array_like
            numeric vector (length p) of S-regression coefficient estimates.
        scale: float
            The S-scale residual estimate
        fitted_values: array_like
            numeric vector (length n) of the fitted values.
        residuals: array_like
            numeric vector (length n) of the residuals.
        rweights: array_like
            numeric vector (length n) of the robustness weights.
        k_iter: int
            (maximal) number of refinement iterations used.
        converged: bool
            boolean indicating if all refinement iterations had converged.
        control: object
            the same object as the control argument.

    See Also
    --------
    
        lmrob
    """
    x = copy.deepcopy(x_)
    y = copy.deepcopy(y_)
    controlS = copy.deepcopy(control)
    nn = x.shape[0]
    pp = x.shape[1]

    n = ctypes.c_int(nn)
    p = ctypes.c_int(pp)
    nResample = ctypes.c_int(0) if only_scale else ctypes.c_int(controlS.nResample)
    large_n = nn > controlS.fast_s_large_n

    # R_lmrob_S params
    residuals = copy.deepcopy(y)
    xp = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    yp = residuals.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    coefficients = np.zeros((pp, 1), dtype="float")
    converged = np.zeros(1, dtype="int")
    scale = np.zeros(1, dtype="float")
    c_chi = _psi_conv_cc(controlS.psi, controlS.tuning_chi)
    #TODO scipy.stats.norm.ppf(3/4)
    c_constant = .6745
    mad = lambda arr,c,center: (np.median(arr) - center) / c
    #scale = ctypes.c_double(mad(x, c_constant, 0) if only_scale else np.zeros(1))
    scale_p = scale.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    coefficients_p = coefficients.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    #c_chi = ctypes.c_double(_psi_conv_cc(controlS.psi, controlS.tuning_chi))
    c_chi_p = c_chi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    psi2ipsi = ctypes.c_int(_psi2ipsi(controlS.psi))
    bb = ctypes.c_double(controlS.bb)
    best_r = ctypes.c_int(controlS.best_r_s)
    groups = ctypes.c_int(controlS.groups)
    n_group = ctypes.c_int(controlS.n_group)
    k_fast_s = ctypes.c_int(controlS.k_fast_s)
    k_iter = ctypes.c_int(controlS.k_max)
    maxit_scale = ctypes.c_int(controlS.maxit_scale)
    refine_tol = ctypes.c_double(controlS.refine_tol)
    inv_tol = ctypes.c_double(controlS.solve_tol)
    scale_tol = ctypes.c_double(controlS.scale_tol)
    converged_p = converged.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    trace_lev = ctypes.c_int(controlS.trace_lev)
    mts = ctypes.c_int(controlS.mts)
    ss = ctypes.c_int(_convSS(controlS.subsampling))
    fast_s_large_n = ctypes.c_int(controlS.fast_s_large_n if large_n else nn + 1)
    R_lmrob_S = lmroblib.R_lmrob_S
    R_lmrob_S.argtypes = [ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int)]
    R_lmrob_S( xp,
               yp,
               ctypes.byref(n),
               ctypes.byref(p),
               ctypes.byref(nResample), 
               scale_p,
               coefficients_p,
               c_chi_p,
               ctypes.byref(psi2ipsi),
               ctypes.byref(bb),
               ctypes.byref(best_r),
               ctypes.byref(groups),
               ctypes.byref(n_group),
               ctypes.byref(k_fast_s),
               ctypes.byref(k_iter),
               ctypes.byref(maxit_scale),
               ctypes.byref(refine_tol),
               ctypes.byref(inv_tol),
               ctypes.byref(scale_tol),
               converged_p,
               ctypes.byref(trace_lev),
               ctypes.byref(mts),
               ctypes.byref(ss),
               ctypes.byref(fast_s_large_n))

    if scale[0] < 0.0:
        raise Exception("C function R_lmrob_S() exited prematurely")
    elif scale[0] == 0.0:
        print("S-estimated scale == 0:  Probably exact fit; check your data")

    if trace_lev.value:
        if only_scale:
            print("lmrob.S(): scale = %g\n" % scale.value)
        else:
            print("lmrob.S(): scale = %g; coeff.=\n" % scale)

    if only_scale:
        return scale

    fitted_values = y - residuals
    rweights = lmrob_rweights(residuals, scale, controlS.tuning_chi, controlS.psi)
    controlS.method = "S"
    resultsDict = {
        "coefficients": coefficients,
        "residuals": residuals,
        "scale": scale,
        "k_iter": np.array([k_iter.value]),
        "converged": converged,
        "fitted_values": fitted_values,
        "rweights": rweights.reshape(( len(rweights), 1)),
        "nResample": controlS.nResample,
        "max_it": controlS.max_it,
        "groups": controlS.groups,
        "best_r_s": controlS.best_r_s,
        "k_fast_s": controlS.k_fast_s,
        "k_max": controlS.k_max,
        "max_it": controlS.max_it,
        "k_m_s": controlS.k_m_s,
        "trace_lev": controlS.trace_lev,
        "mts": controlS.mts,
        "compute_rd": controlS.compute_rd,
        "numpoints": controlS.numpoints,
        "fast_s_large_n": controlS.fast_s_large_n,
        "tuning_chi":controlS.tuning_chi,
        "bb": bb.value,
        "tuning_psi": controlS.tuning_psi,
        "refine_tol": controlS.refine_tol,
        "rel_tol": controlS.rel_tol,
        "scale_tol": controlS.scale_tol,
        "solve_tol": controlS.solve_tol,
        "warn_limit_reject": controlS.warn_limit_reject,
        "psi": controlS.psi,
        "subsampling": controlS.subsampling,
        "cov": controlS.cov,
        "split_type": controlS.split_type,
        "control": controlS

        }

    return resultsDict

def lmrob_fit(x, y, control, init):
    """
    Compute MM-type estimators of regression: An S-estimator is used as starting value, and an M-
    estimator with fixed scale and redescending psi-function is used from there. Optionally a D-step
    (Design Adaptive Scale estimate) as well as a second M-step is calculated. Returns a dictionary
    object.

    Parameters
    ----------
    x: array_like
        design matrix (n × p) typically including a column of 1s for the intercept.
    y: array_like
        numeric response vector (of length n).
    control: object 
        An instance of the LmrobControl Class
    init: dict
        Optional list of initial estimates. *See Notes.*

    
    Examples
    --------
    >>> from lmrob import *
    >>> from rpy2.robjects.packages import importr
    >>> from rpy2.robjects import pandas2ri, r
    >>> import rpy2.robjects.numpy2ri as rpyn
    
    >>> stackloss = r("stackloss")
    >>> data = {"AirFlow": rpyn.ri2py(stackloss.rx2("Air.Flow")),
    ...         "WaterTemp": rpyn.ri2py(stackloss.rx2("Water.Temp")),
    ...         "AcidConc": rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
    ...         "stack_loss": rpyn.ri2py(stackloss.rx2("stack.loss"))
    ...        }
    >>> formula = 'stack_loss ~ AirFlow + WaterTemp + AcidConc'
    >>> Y, X = patsy.dmatrices(formula, data, NA_action="drop",
    ...                        return_type='matrix')
    >>> control = LmrobControl()
    >>> m0 = lmrob_fit(X, Y, control=control, init=None)
    >>> coefficients = m0.get("coefficients")
    >>> residuals = m0.get("residuals")
    >>> scale = m0.get("scale")
    >>> converged = m0.get("converged")

    Notes
    -----

    This function is the basic fitting function for MM-type estimation, called by lmrob and typically
    not to be used on its own.
    If given, init must be a list of initial estimates containing at least the initial coefficients and scale
    as coefficients and scale. Otherwise it calls ``lmrob_S(..)`` and uses it as initial estimator.
    
    Returns
    --------

    fitted_values: array_like
        Xβ, i.e., ``X.dot(coefficients)``

    residuals: array_like
        the raw residuals, ``y - fitted_values``
    rweights: array_like
        robustness weights derived from the final M-estimator residuals (even when not
        converged).    
    rank: int
        degree_freedom n - rank
    coefficients: array_like
        estimated regression coefficient vector.
    scale: float
        the robustly estimated error standard deviation
    cov: array_like
        variance-covariance matrix of coefficients, if the RWLS iterations have con-
        verged (and control.cov is not "none").
    control: object
    iter: int
    converged: bool
        boolean indicating if the RWLS iterations have converged.
    init: dict
        A similar dict that contains the results of intermediate estimates (not for MM-
        estimates).
    

    See Also
    --------

    lmrob, lmrob__M__fit, lmrob__D__fit, lmrob_S

    """
    if control.method == "MM":
        control.method = "SM"
    ## Assumption:  if(is.null(init))  method = "S..."   else  method = "..."
    ## ---------    where "..." consists of letters {"M", "D"}
    if init:
        if init.__class__ != dict:
            raise Exception("init must be a dictionary")
        if init["converged"] is None:
            init["converged"] = True
        M = est = init["control"].method
        if init["control"] is None:
            init["control"] = control
            est = init["control"].method = "l"
        elif not len(M) or not M.__class__ == str:
            est = "l"
    else:
        M1 = control.method[0]
        if control.method[0] != "S":
            #TODO: warning using M1
            control.method = "S" + control.method[1:]
        init = lmrob_S(x, y, control)
        est = "S"
    if not "float" in init["coefficients"].dtype.name or not len(init["coefficients"]) == x.shape[1] or not "float" in init["scale"].dtype.name or not init["scale"] >= 0:
        raise Exception("some inconsistency")
    if est != "S" and control.cov == "_vcov_avar1":
        # TODO: warning
        control.cov = "_vcov_w"
    trace_lev = control.trace_lev
    if init["converged"]:
        ## --- loop through the other estimators; build up 'est' string
        #control.method[1:] if control.method[0] == est else control.method
        method = control.method[0:len(est)].replace(est, "") + control.method[len(est):]
        if trace_lev:
            print("init converged (remaining method = '%s') -> coef=\n" %(method))
            print(init["coefficients"])

        for step in method:
            ## now we have either M or D steps
            est += step
            ## 'control' may differ from 'init$control' when both (init, control) are spec. 
            ## D(AS)-Step
            if step == "D":
                init = lmrob__D__fit(init, x, control=control, method=init["control"].method)
            ## M-Step
            elif step == "M":
                init = lmrob__M__fit(x, y, obj=init, control=control, method=init["control"].method)
            else:
                raise Exception("only M and D are steps supported after 'init' computation")

            if trace_lev:
                print("step '%s' -> new coef=\n" %step)
                print(init["coefficients"])
            if not init["converged"]:
                #TODO: warning
                break

    control.method = est ## ~= original 'method', but only with the steps executed.
    init["control"] = control

    ## --- covariance estimate
    # TODO init["cov"] 
    # if init["scale"] == 0:
    #     init["cov"] = np.zeros_like(x)
    # elif not init["converged"] or x == None:
    #     init["cov"] = np.nan
    # else:
    #     if control.cov == None or control.cov == "none":
    #         init["cov"] = np.nan
        # TODO: else:
    # TODO: df = y.shape[0] - init["rank"] ## sum(init$r?weights)-init$rank
    # TODO: init$degree.freedom <- init$df.residual <- df
    return init 
    ## end{lmrob.fit}

def lmrob__D__fit(obj, x=np.array([]), control=None, method=None):
    """
    This function calculates a Design Adaptive Scale estimate for a given MM-estimate. This is sup-
    posed to be a part of a chain of estimates like SMD or SMDM. Returns a dictionart object
    
    Parameters
    ----------
    obj: dict
        based on which the estimate is to be calculated.
    x: array_like
        The design matrix; if missing, the method tries to get it from ``obj.get("x")``.
    control: object
        Instances of class LmrobControl.
    method: str
        (optional) the method used for obj computation.
    
    Returns
    -------

    scale: float
        The Design Adaptive Scale estimate

    converged: bool
        True if the scale calculation converged, False other.

    Examples
    --------
    >>> from lmrob import *
    >>> from rpy2.robjects.packages import importr
    >>> from rpy2.robjects import pandas2ri, r
    >>> import rpy2.robjects.numpy2ri as rpyn
    >>>
    >>> stackloss = r("stackloss")
    >>> data = {"AirFlow": rpyn.ri2py(stackloss.rx2("Air.Flow")),
    ...         "WaterTemp": rpyn.ri2py(stackloss.rx2("Water.Temp")),
    ...         "AcidConc": rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
    ...         "stack_loss": rpyn.ri2py(stackloss.rx2("stack.loss"))
    ...      }
    >>> formula = 'stack_loss ~ AirFlow + WaterTemp + AcidConc'
    >>> Y, X = patsy.dmatrices(formula, data, NA_action="drop",
    ...                        return_type='matrix')
    >>> control = LmrobControl()
    >>> init = lmrob_S(X, Y, control)
    >>> m0 = lmrob__D__fit(init, X, control)
    >>> # Getting the computed values
    >>> converged_py = m0.get("converged")
    >>> scale_py = m0.get("scale") 
    
    Notes
    -----
    This function is used by ``lmrob_fit`` and typically not to be used on its own. Note that ``lmrob_fit()``
    specifies control potentially differently than the default, but does use the default for method.

    See Also
    --------
    lmrob_fit, lmrob

    """
    objD = copy.deepcopy(obj)
    if control is None:
        control = objD.get("control")
        if control is None:
            Exception('lmrob__D::fit: control is missing')
    if method is None:
        method_ = objD.get("control").method
    else:
        method_ = copy.deepcopy(method)
    if not objD.get("converged"):
        Exception('lmrob..D..fit: prior estimator did not converge, stopping')
    controlD = copy.deepcopy(control)
    
    if x.size == 0:
        x_ = copy.deepcopy(objD.get("x"))
    else:
        x_ = copy.deepcopy(x)

    w_ = objD.get("rweights", np.array([])).ravel()
    if w_.size == 0:
        Exception('lmrob..D..fit: robustness weights undefined')

    r_ = objD.get("residuals", np.array([])).ravel()
    if r_.size == 0:
        Exception('lmrob..D..fit: residuals undefined')

    psi = controlD.psi
    if psi is None:
        Exception('lmrob..D..fit: parameter psi is not defined')

    c_psi_ = _psi_conv_cc(psi, controlD.tuning_chi if method_ in ['S', 'SD'] else controlD.tuning_psi)
    if not isinstance(c_psi_, (int,float)):
        Exception('lmrob..D..fit: parameter tuning_psi is not numeric')

    if objD.get("kappa") is None:
        objD.update({"kappa": lmrob_kappa(objD, controlD)})
    kappa_ = objD.get("kappa")

    if objD.get("tau") is None:
        objD.update({"tau": lmrob_tau(objD, x_, control)})
    tau_ = objD.get("tau")
    
    length_ = r_.size
    ipsi_ = _psi2ipsi(psi)
    type_ = 3
    rel_tol_ = controlD.rel_tol 
    k_max_ = controlD.k_max
    converged_ = 0
    
    scale_1_ = np.sqrt(np.sum(w_ * r_ ** 2) / kappa_ / np.sum(tau_**2 * w_))
    scale_1_ = np.array([scale_1_], dtype="float")
    R_find_D_scale = lmroblib.R_find_D_scale

    r_p = r_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    kappa_p = kappa_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    tau_p = tau_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    length = ctypes.c_int(length_)
    scale_1 = scale_1_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    #c = ctypes.c_double(c_psi_)
    #c_chi = _psi_conv_cc(psi, cc)
    c_p = c_psi_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ipsi = ctypes.c_int(ipsi_)
    Type = ctypes.c_int(type_)
    rel_tol = ctypes.c_double(rel_tol_)
    k_max = ctypes.c_int(k_max_)
    converged = ctypes.c_int(converged_)
    
    R_find_D_scale.argtypes = [ ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int)]
    R_find_D_scale(r_p,
                   kappa_p,
                   tau_p,
                   ctypes.byref(length),
                   scale_1,
                   c_p,
                   ctypes.byref(ipsi),
                   ctypes.byref(Type),
                   ctypes.byref(rel_tol),
                   ctypes.byref(k_max),
                   ctypes.byref(converged))

    objD.update({"scale": scale_1_ if converged else None})
    objD.update({"converged": np.bool(converged.value)})

    if not "D" == method_[-1]:
        if method_ == "MM":
            method_ = "SM"
        method_ = method_ + "D"


    objD.update({"control": controlD})
    if objD.get("cov"):
        if controlD.cov == _vcov_avar1:
            controlD.cov = _vcov_w
        lf_cov = control.cov
        #objD.update({"cov": lf_cov(objD, x)})

    return objD

def lmrob__M__fit(x=np.array([]), y=np.array([]),
                  beta_initial=np.array([]), scale=None,
                  control=None, method=None, obj=None):
    """
    This function performs RWLS iterations to find an M-estimator of regression. When started from
    an S-estimated beta.initial, this results in an MM-estimator. Returns a dictionary

    Parameters
    ----------
    x: array_like
        design matrix (n × p) typically including a column of 1s for the intercept.
    y: array_like
        numeric response vector (of length n).
    beta_initial: array_like
        numeric vector (of length p) of initial estimate. Usually the result of an S-
        regression estimator.
    scale: float
        robust residual scale estimate. Usually an S-scale estimator.
    control: object
        Instances of class LmrobControl.
    obj: dict
        If specified, this is typically used to set values for the
        other arguments.
    method: str
        (optional) the method used for obj computation.

    Returns
    -------
    
    coefficients: array_like
        the M-estimator (or MM-estim.) of regression
    control:
        the control obj input used
    scale: float
        The residual scale estimate
    converged: bool
        True if the RWLS iterations converged, False otherwise

    Examples
    --------
    >>> from lmrob import *
    >>> from rpy2.robjects.packages import importr
    >>> from rpy2.robjects import pandas2ri, r
    >>> import rpy2.robjects.numpy2ri as rpyn
    >>>
    >>> stackloss = r("stackloss")
    >>> data = {"AirFlow": rpyn.ri2py(stackloss.rx2("Air.Flow")),
    ...         "WaterTemp": rpyn.ri2py(stackloss.rx2("Water.Temp")),
    ...         "AcidConc": rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
    ...         "stack_loss": rpyn.ri2py(stackloss.rx2("stack.loss"))
    ...      }
    >>> formula = 'stack_loss ~ AirFlow + WaterTemp + AcidConc'
    >>> Y, X = patsy.dmatrices(formula, data, NA_action="drop",
    ...                        return_type='matrix')
    >>> control = LmrobControl()
    >>> init = lmrob_S(X, Y, control)
    >>> beta_inital = init.get("coefficients")
    >>> scale = init.get("scale")
    >>> m0 = lmrob__M__fit(X, Y, beta_inital, scale,
    ...                    control, method="MM")
    >>> # Getting the computed values
    >>> coefficients = m0.get("coefficients")
    >>> residuals = m0.get("residuals")
    >>> scale = m0.get("scale")
    >>> converged = m0.get("converged")
    >>> iter = m0.get("iter")

    Notes
    -----
        This function is used by lmrob_fit and typically
        not to be used on its own.

    See Also
    --------
    lmrob_fit, lmrob
    """
    if not obj is None:
        objM = copy.deepcopy(obj)
    else:
        objM = {}
    if x.size == 0:
        x_ = objM.get("x")
    else:
        x_ = copy.deepcopy(x)
    if y.size == 0:
        y_ = objM.get("y")
    else:
        y_ = copy.deepcopy(y)
    if beta_initial.size == 0:
        beta_initial_ = objM.get("coefficients")
    else:
        beta_initial_ = copy.deepcopy(beta_initial)
    if scale is None:
        scale_ = objM.get("scale")
    else:
        scale_ = copy.deepcopy(scale)
    if control is None:
        control = objM.get("control")
    else:
        controlM = copy.deepcopy(control)
    if method is None:
        method_ = copy.deepcopy(objM.get("control").method)
    else:
        method_ = copy.deepcopy(method)

    coefficients = np.zeros((x.shape[1], 1))
    residuals = np.zeros((x.shape[0], 1))
    # R_lmrob_MM params
    x_p = x_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    y_p = y_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    n = ctypes.c_int(x.shape[0])
    p = ctypes.c_int(x.shape[1])
    beta_initial_p = beta_initial_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    scale_p = scale_.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    coefficients_p = coefficients.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    residuals_p = residuals.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    iterations = ctypes.c_int(controlM.max_it)

    #c_psi = ctypes.c_double(_psi_conv_cc(controlM.psi, controlM.tuning_psi))
    controlM.tuning_psi = np.array([controlM.tuning_psi], dtype="float") if np.isscalar(controlM.tuning_psi) else controlM.tuning_psi
    c_psi = _psi_conv_cc(controlM.psi, controlM.tuning_psi)
    c_psi_p = c_psi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ipsi = ctypes.c_int(_psi2ipsi(controlM.psi))

    loss = ctypes.c_double(1.0)
    rel_tol = ctypes.c_double(controlM.rel_tol)
    converged = ctypes.c_int(1)
    trace_lev = ctypes.c_int(controlM.trace_lev)
    mts = ctypes.c_int(controlM.mts)
    ss = ctypes.c_int(_convSS(controlM.subsampling))
    R_lmrob_MM = lmroblib.R_lmrob_MM
    R_lmrob_MM.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    ]
    
    R_lmrob_MM( x_p,
                y_p,
                ctypes.byref(n),
                ctypes.byref(p),
                beta_initial_p,
                scale_p,
                coefficients_p,
                residuals_p,
                ctypes.byref(iterations),
                c_psi_p,
                ctypes.byref(ipsi),
                ctypes.byref(loss),
                ctypes.byref(rel_tol),
                ctypes.byref(converged),
                ctypes.byref(trace_lev),
                ctypes.byref(mts),
                ctypes.byref(ss)
               )
    objM.update({ "coefficients": coefficients,
                  "scale": scale_,
                  "residuals": residuals,
                  "loss": loss.value,
                  "converged": converged.value ,
                  "iter": iterations.value
                })
    fitted_values = x.dot(coefficients)
    rweights = lmrob_rweights(residuals, scale_, controlM.tuning_psi, controlM.psi)
    objM.update({"fitted_values": fitted_values,
                 "rweights": rweights.reshape(( len(rweights), 1))})
    objM.update({"control": controlM})
    if not method_[-1] == "M":
        method_ = method_ + "M"
    qr = np.linalg.qr(x * np.sqrt(rweights).reshape(x.shape[0], 1))
    rank = np.linalg.matrix_rank(qr[1])
    if not objM.get("cov") is None:
        if not method_ in ["SM","MM"] and objM.get("control").cov == _vcov_avar1:
            objM.get("control").cov = _vcov_w
        lf_cov = objM.get("control").cov
        #objM.update({"cov" : lf_cov(objM, x)})
    return objM


def lmrob_M_S(x, y, control, mf):
    """

    Computes an M-S-estimator for linear regression using the "M-S" algorithm.
    
    Parameters
    ----------

    x: array_like
        Design matrix (n × p)
    y: array_like
        Numeric vector of responses (or residuals for only_scale=True).
    control: object
        Instances of class LmrobControl.
    mf: array_like
        a model matrix as returned by model_frame.
    
    Returns
    coefficients: array_like
        Numeric vector (length p) of M-S-regression coefficient estimates.
    scale: float
        The M-S-scale residual estimate
    residuals: array_like
        numeric vector (legnth n) of the residuals.
    rweights: array_like
        numeric vector (length n) of the robustness weights.
    control: object
        the same as the control argument.
    converged: bool
        Convergence status (always True), needed for lmrob.fit.

    Examples
    --------
    >>> from lmrob import *
    >>> from rpy2.robjects.packages import importr
    >>> from rpy2.robjects import pandas2ri, r
    >>> import rpy2.robjects.numpy2ri as rpyn
    >>>
    >>> stackloss = r("stackloss")
    >>> data = {"AirFlow" : rpyn.ri2py(stackloss.rx2("Air.Flow")),
    ...         "WaterTemp" : rpyn.ri2py(stackloss.rx2("Water.Temp")),
    ...         "AcidConc" : rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
    ...         "stack_loss" : rpyn.ri2py(stackloss.rx2("stack.loss"))
    ...      }
    >>> formula = 'stack_loss ~ AirFlow + WaterTemp + AcidConc'
    >>> mf_dict = model_frame(formula, data, None, [], 'drop', [])
    >>> Y, X = patsy.dmatrices(formula, data, NA_action="drop", return_type='matrix')
    >>> control = LmrobControl()
    >>> m0 = lmrob_M_S(X, Y, control, mf_dict)
    >>> coefficients = m0.get("coefficients")
    >>> residuals = m0.get("residuals")
    >>> scale = m0.get("scale")

    Notes
    -----
    This function is used by lmrob and not intended to be used on its own (because an M-S-estimator
    has too low efficiency ‘on its own’).
    An M-S estimator is a combination of an S-estimator for the continuous variables and an L1-
    estimator (i.e. an M-estimator with ψ(t) = sign(t)) for the categorical variables.
    The S-estimator is estimated using a subsampling algorithm. If the model includes interactions
    between categorical (factor) and continuous variables, the subsampling algorithm might fail. In
    this case, one can choose to assign the interaction to the categorical side of variables rather than
    to the continuous side.

    Note that the return status converged does not refer to the actual convergence status. The algorithm
    used does not guarantee convergence and thus true convergence is almost never reached. This is,
    however, not a problem if the estimate is only used as initial estimate part of an MM or SMDM
    estimate.

    The algorithm sometimes produces the warning message "Skipping design matrix equilibration
    (dgeequ): row ?? is exactly zero.". This is just an artifact of the algorithm and can be ignored
    safely.

    See Also
    --------
    lmrob

    """
    # Make a copy of input variables -----
    x_copy = copy.deepcopy(x)
    mf_copy = copy.copy(mf)
    controlMS = copy.deepcopy(control)

    # Get a split frame of variables categoricals and continuos -----
    split = split_frame(mf_copy, x_copy, controlMS.split_type)
    if split.get("x1").shape[1] == 0:
        print("Warning: No categorical variables found in model. Reverting to S-estimator.")
        return lmrob_S(x, y, controlMS) 
    if split.get("x2").shape[1] == 0:
        print("No continuous variables found in model. Reverting to L1-estimator.")
        return lmrob_lar(x, y, controlMS) 

    # TODO: ## this is the same as in lmrob.S(): ...
    x1 = np.asfortranarray(split.get("x1"), dtype="float")
    x2 = np.asfortranarray(split.get("x2"), dtype="float")
    y_copy = np.asfortranarray(copy.deepcopy(y), dtype="float")
    n_ = y.shape[0]
    p1_ = x1.shape[1]
    p2_ = x2.shape[1]
    scale = np.zeros(1, dtype="float")
    residuals = np.zeros((n_, 1), dtype="float")
    coefficients1 = np.zeros((p1_, 1), dtype="float")
    coefficients2 = np.zeros((p2_, 1), dtype="float")
    converged = np.zeros(1, dtype="bool")
    orthogonalize = np.ones(1, dtype="bool")
    subsample = np.ones(1, dtype="bool")
    descent = np.ones(1, dtype="bool")    

    x1_p = x1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    x2_p = x2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    y_p = y_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    residuals_p = residuals.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    n = ctypes.c_int(n_)
    p1 = ctypes.c_int(p1_)
    p2 = ctypes.c_int(p2_)
    nResample = ctypes.c_int(controlMS.nResample)
    maxit_scale = ctypes.c_int(controlMS.maxit_scale)
    scale_p = scale.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    coefficients1_p = coefficients1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    coefficients2_p = coefficients2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_chi = _psi_conv_cc(controlMS.psi, controlMS.tuning_chi)
    c_chi_p = c_chi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ipsi = ctypes.c_int(_psi2ipsi(controlMS.psi))
    bb = ctypes.c_double(controlMS.bb)
    k_m_s = ctypes.c_int(controlMS.k_m_s)
    max_k = ctypes.c_int(controlMS.k_max)
    rel_tol = ctypes.c_double(controlMS.rel_tol)
    inv_tol = ctypes.c_double(controlMS.solve_tol)
    scale_tol = ctypes.c_double(controlMS.scale_tol)
    converged_p = converged.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    trace_lev = ctypes.c_int(0)
    orthogonalize_p = orthogonalize.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    subsample_p = subsample.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    descent_p = descent.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    mts = ctypes.c_int(controlMS.mts)
    ss = ctypes.c_int(_convSS(controlMS.subsampling))
    R_lmrob_M_S = lmroblib.R_lmrob_M_S
    R_lmrob_M_S.argtypes = [ ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int)
                ]
    R_lmrob_M_S(x1_p,
                x2_p,
                y_p,
                residuals_p,
                ctypes.byref(n),
                ctypes.byref(p1),
                ctypes.byref(p2),
                ctypes.byref(nResample),
                ctypes.byref(maxit_scale),
                scale_p,
                coefficients1_p,
                coefficients2_p,
                c_chi_p,
                ctypes.byref(ipsi),
                ctypes.byref(bb),
                ctypes.byref(k_m_s),
                ctypes.byref(max_k),
                ctypes.byref(rel_tol),
                ctypes.byref(inv_tol),
                ctypes.byref(scale_tol),
                converged_p,
                ctypes.byref(trace_lev),
                orthogonalize_p,
                subsample_p,
                descent_p,
                ctypes.byref(mts),
                ctypes.byref(ss))

    ## FIXME? warning if 'conv' is not ok ??
    ## coefficients :
    idx = split.get("x1_idx")
    not_idx = idx == False
    coef = np.zeros(shape=(len(idx), 1), dtype="float")
    coef[    idx, :] = coefficients1
    coef[not_idx, :] = coefficients2

    ## set method argument in control
    controlMS.method = 'M-S'

    ## Get rweights
    rweights = lmrob_rweights(residuals, scale, controlMS.tuning_chi, controlMS.psi)

    resultsDict = {
        "coefficients": coef,
        "scale": scale,
        "residuals": residuals,
        "rweights": rweights,
        "converged": True,
        "descent_conv": converged,
        "control": controlMS
    }

    # TODO: obj$ostats <- outlierStats(obj, x, control)

    return resultsDict

def lmrob_lar(x, y, control=None):
    """
    Compute MM-type estimators of regression: An S-estimator is used as starting value, and an M-
    estimator with fixed scale and redescending psi-function is used from there. Optionally a D-step
    (Design Adaptive Scale estimate) as well as a second M-step is calculated.
    Usage. Return a dictionary object.
    

    Parameters
    ----------

    x: array_like
        Design matrix (n × p)
    y: array_like
        Numeric vector of responses (or residuals for only_scale=True).
    control: object
        Instances of class LmrobControl.
    
    
    Returns
    -------

    fitted_values: array_like
        Xβ, i.e., ``X.dot(coefficients)``.
    residuals: array_like
        the raw residuals, ``y - fitted_values``
    rweights: array_like
        robustness weights derived from the final M-estimator residuals (even when not
        converged).
    rank: int
        ``degree_freedom n - rank``
    coefficients: array_like
        estimated regression coefficient vector.
    scale: float
        the robustly estimated error standard deviation.
    control: object
        the same as the control argument.
    iter: bool
        converged boolean indicating if the RWLS iterations have converged.
    init: dict
        A similar dict that contains the results of intermediate estimates (not for MM-
        estimates).    
    
    Examples
    --------
    >>> from lmrob import *
    >>> from rpy2.robjects.packages import importr
    >>> from rpy2.robjects import pandas2ri, r
    >>> import rpy2.robjects.numpy2ri as rpyn
    >>>
    >>> stackloss = r("stackloss")
    >>> data = {"AirFlow" : rpyn.ri2py(stackloss.rx2("Air.Flow")),
    ...         "WaterTemp" : rpyn.ri2py(stackloss.rx2("Water.Temp")),
    ...         "AcidConc" : rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
    ...         "stack_loss" : rpyn.ri2py(stackloss.rx2("stack.loss"))
    ...      }
    >>> formula = 'stack_loss ~ AirFlow + WaterTemp + AcidConc'
    >>> mf_dict = model_frame(formula, data, None, [], 'drop', [])
    >>> Y, X = patsy.dmatrices(formula, data, NA_action="drop", return_type='matrix')
    >>> control = LmrobControl()
    >>> m0 = lmrob_lar(X, Y, control)
    >>> coefficients = m0.get("coefficients")
    >>> residuals = m0.get("residuals")
    >>> scale = m0.get("scale")


    Notes
    -----

    This function is the basic fitting function for MM-type estimation, called by lmrob and typically
    not to be used on its own.
    If given, init must be a list of initial estimates containing at least the initial coefficients and scale
    as coefficients and scale. Otherwise it calls lmrob.S(..) and uses it as initial estimator.


    See Also
    --------
    lmrob, lmrob..M..fit, lmrob..D..fit, lmrob.S

    """
    if control is None:
        control = LmrobControl()
    controlLar = control.copy()
    x_copy = np.asfortranarray(x)
    y_copy = np.asfortranarray(y)
    p_ = x_copy.shape[1]
    n_ = x_copy.shape[0]
    rel_tol_ = controlLar.rel_tol
    NIT = np.zeros(1, dtype="int")
    K = np.zeros(1, dtype="int")
    KODE = np.zeros(1, dtype="int")
    SIGMA = np.zeros(1)
    THETA = np.zeros(n_)
    RS = np.zeros(n_)
    SC1 = np.zeros(n_)
    SC2 = np.zeros(p_)
    SC3 = np.zeros(p_)
    SC4 = np.zeros(p_)
    bet0_ = 0.773372647623

    x_p = x_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    y_p = y_copy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    p = ctypes.c_int(p_)
    n = ctypes.c_int(n_)
    MDX = ctypes.c_int(n_)
    MDT = ctypes.c_int(n_)
    rel_tol = ctypes.c_double(rel_tol_)
    NIT_p = NIT.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    K_p = K.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    KODE_p = KODE.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    SIGMA_p = SIGMA.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    THETA_p = THETA.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    RS_p = RS.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    SC1_p = SC1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    SC2_p = SC2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    SC3_p = SC3.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    SC4_p = SC4.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    bet0 = ctypes.c_double(bet0_)

    if not (p_ > 0 and n_ >= p_ and y.size == n_):
        Exception("")


    rllarsbi = lmroblib.rllarsbi_
    rllarsbi.argtypes = [ ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double)]
    rllarsbi( x_p,
              y_p,
              ctypes.byref(n),
              ctypes.byref(p),
              ctypes.byref(MDX),
              ctypes.byref(MDT),
              ctypes.byref(rel_tol),
              NIT_p,
              K_p,
              KODE_p,
              SIGMA_p,
              THETA_p,
              RS_p,
              SC1_p,
              SC2_p,
              SC3_p,
              SC4_p,
              ctypes.byref(bet0))

    obj = { "coefficients": THETA[:p_],
            "scale": SIGMA,
            "residuals": RS,
            "iter": NIT,
            "status": KODE,
            "converged": True,
            }

    return obj

if __name__ == '__main__':
    lmrob_control = LmrobControl(None)