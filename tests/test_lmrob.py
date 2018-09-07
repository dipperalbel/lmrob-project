import pandas
import sys
import os
import inspect
import patsy            
# R packages (provisional)
from rpy2.robjects.functions import SignatureTranslatedFunction
from rpy2 import robjects
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as rpyn
import unittest
path_of_script = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe() ))[0]))
sys.path.append(os.path.join(path_of_script, "../"))
from lmrob import *
from lmrob import _psi2ipsi, _lmrob_hat
from scipy.optimize import root

class LmrobTestCase(unittest.TestCase):


    def setUp(self):
        base = importr('base')
        self.stats = importr('stats')
        self.robustbase = importr('robustbase')
        STM = SignatureTranslatedFunction

        stackloss = r("stackloss")
        pandas2ri.activate()
        robjects.globalenv['dataframe'] = stackloss
        self.x = base.matrix(base.c(stackloss.rx2("Water.Temp"), stackloss.rx2("Air.Flow")), ncol=2 )
        self.y = base.matrix(stackloss.rx2("stack.loss"))
        
        # load stackloss data
        self.stackloss = stackloss
        self.dictStackloss = {
            "AirFlow"   : rpyn.ri2py(stackloss.rx2("Air.Flow")),
            "WaterTemp" : rpyn.ri2py(stackloss.rx2("Water.Temp")),
            "AcidConc"  : rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
            "stack_loss" : rpyn.ri2py(stackloss.rx2("stack.loss"))
        }

        # load categorical data
        pandas2ri.activate()
        r['load']('categorical_data.rda')
        self.categorical_data = r['data']

        self.robustbase.lmrob_control = STM(self.robustbase.lmrob_control, init_prm_translate = {'tuning.chi': 'tuning_chi', 'tuning.psi': 'tuning_psi'})
        #r['load']('categorical_control.rda')
        #self.categoricalal_control = r['ctrl']

    def assertNorm(self, m1, m2, messages):
        m1 = m1.ravel()
        m2 = m2.ravel()
        m1 = np.sort(m1[np.isnan(m1) == False])
        m2 = np.sort(m2[np.isnan(m2) == False])
        self.assertEqual(m1.size, m2.size, "Incorrect Size")
        size = m1.size
        error = np.linalg.norm(m1 - m2) / size
        tol = 1e-4
        self.assertLess(error, tol, messages)


    def test_lm_fit(self):
        X = self.x
        Y = self.y
        # R from python
        m0 = self.stats._lm_fit(X, Y)
        qr = rpyn.ri2py((m0.rx2('qr')))
        coefficients = rpyn.ri2py((m0.rx2('coefficients')))
        residuals = rpyn.ri2py((m0.rx2('residuals')))
        effects = rpyn.ri2py((m0.rx2('effects')))
        rank = rpyn.ri2py((m0.rx2('rank')))
        pivot = rpyn.ri2py((m0.rx2('pivot')))
        qraux = rpyn.ri2py((m0.rx2('qraux')))
        tol = rpyn.ri2py((m0.rx2('tol')))
        pivoted = rpyn.ri2py((m0.rx2('pivoted')))

         # Developed in python
        m0_py = lm_fit(rpyn.ri2py(X), rpyn.ri2py(Y))
        coefficients_py = m0_py['coefficients']
        residuals_py = m0_py['residuals']
        effects_py = m0_py['effects']
        rank_py = m0_py['rank']
        pivot_py = m0_py['pivot']
        qraux_py = m0_py['qraux']
        pivoted_py = m0_py['pivoted']

        self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients")
        self.assertNorm(residuals, residuals_py, "Incorrect residuals")
        self.assertNorm(effects, effects_py, "Incorrect effects")
        self.assertEqual(rank, rank_py, "Incorrect rank")
        self.assertNorm(pivot, pivot_py + 1, "Incorrect pivot")
        self.assertNorm(qraux, qraux_py, "Incorrect qraux")
        self.assertEqual(pivoted, pivoted_py, "Incorrect pivoted")

    def test_lmrob_S(self):
        X = self.x
        Y = self.y
        control = self.robustbase.lmrob_control()

        m0 = self.robustbase.lmrob_S(self.x, self.y, control)

        #print(control)
        coefficients = rpyn.ri2py((m0.rx2('coefficients')))
        residuals = rpyn.ri2py((m0.rx2('residuals')))
        scale = rpyn.ri2py((m0.rx2('scale')))
        converged = rpyn.ri2py((m0.rx2('converged')))
        k_iter = rpyn.ri2py((m0.rx2('k.iter')))
        fitted_values = rpyn.ri2py(m0.rx2('fitted.values'))
        rweights = rpyn.ri2py(m0.rx2('rweights'))

        control_py = LmrobControl()
        m0_py = lmrob_S(rpyn.ri2py(X), rpyn.ri2py(Y),control_py)
        coefficients_py = m0_py.get("coefficients")
        residuals_py = m0_py.get("residuals")
        scale_py = m0_py.get("scale")
        converged_py = m0_py.get("converged")
        k_iter_py = m0_py.get("k_iter")
        fitted_values_py = m0_py.get("fitted_values")
        rweights_py = m0_py.get("rweights")


        self.assertNorm(scale, scale_py, "Incorrect scale")
        self.assertEqual(converged, converged_py, "Incorrect converged")
        self.assertNorm(k_iter, k_iter_py, "Incorrect k_iter")
        self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients")
        self.assertNorm(residuals, residuals_py, "Incorrect residuals")
        self.assertNorm(fitted_values, fitted_values_py, "Incorrect fitted_values")
        self.assertNorm(rweights, rweights_py, "Incorrect rweights")

    def test_lmrob_M_fit(self):
        X = self.x
        Y = self.y
        control = self.robustbase.lmrob_control()
        init = self.robustbase.lmrob_S(self.x, self.y, control)
        m0 = self.robustbase.lmrob__M__fit(self.x, self.y, obj=init)
        
        coefficients = rpyn.ri2py(m0.rx2("coefficients"))
        scale = m0.rx2("scale")
        residuals = rpyn.ri2py(m0.rx2("residuals"))
        loss = rpyn.ri2py(m0.rx2("loss"))
        converged = rpyn.ri2py(m0.rx2("converged"))
        iterations = rpyn.ri2py(m0.rx2("iter"))

        control_py = LmrobControl()

        rpyn.ri2py(X), rpyn.ri2py(Y),control_py
        beta_inital = rpyn.ri2py(init.rx2("coefficients"))
        scale = rpyn.ri2py(init.rx2("scale"))
        m0_py = lmrob__M__fit(rpyn.ri2py(X), rpyn.ri2py(Y), beta_inital, scale, control_py, "MM")
        coefficients_py = m0_py.get("coefficients")
        residuals_py = m0_py.get("residuals")
        scale_py = m0_py.get("scale")
        converged_py = m0_py.get("converged")
        loss_py = m0_py.get("loss")
        iter_py = m0_py.get("iter")


        self.assertEqual(scale, scale_py, "Incorrect scale")
        self.assertEqual(converged, converged_py, "Incorrect converged")
        self.assertEqual(iterations, iter_py, "Incorrect iter")
        self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients")
        self.assertNorm(residuals, residuals_py, "Incorrect residuals")


    def test_lmrob_M_S(self):
        categorical_data = self.categorical_data
        categorical_data_py = pandas.DataFrame(self.categorical_data)

        # method = "S" with init = "M-S" makes to lmrob execute lmrob_M_S only
        method = "S"
        init   = "M-S"

        # R's lmrob
        for split_type in ["f", "fi", "fii"]:
            control = self.robustbase.lmrob_control(psi='optimal', tuning_chi=20.0, bb = 0.0003846154, tuning_psi=20.0, method=method , cov='.vcov.w', split_type=split_type)
            m0_R = self.robustbase.lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data, init=init, control=control)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            control_py = LmrobControl(psi="optimal", tuning_chi = 20, bb = 0.0003846154, tuning_psi=20, method=method, cov="_vcov_w", split_type=split_type)
            m0_py = lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data_py, init=init, control=control_py)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            
            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)

            # print("coefficients R:\n", coefficients)
            # print("coefficients py:\n", coefficients_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method) 

    def test_lmrob_fit(self):
        X = self.x
        Y = self.y
        control = self.robustbase.lmrob_control()
        control_py = LmrobControl()
        control_py.method = "MM"


        m0 = self.robustbase.lmrob_fit(X, Y, control)
        coefficients = rpyn.ri2py(m0.rx2("coefficients"))
        scale = rpyn.ri2py(m0.rx2("scale"))
        residuals = rpyn.ri2py(m0.rx2("residuals"))
        converged = rpyn.ri2py(m0.rx2("converged"))

        m0_py = lmrob_fit(rpyn.ri2py(X), rpyn.ri2py(Y), control_py, init=None)
        coefficients_py = m0_py.get("coefficients")
        residuals_py = m0_py.get("residuals")
        scale_py = m0_py.get("scale")
        converged_py = m0_py.get("converged")

        self.assertNorm(scale, scale_py, "Incorrect scale")
        self.assertEqual(converged, converged_py, "Incorrect converged")
        self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients")
        self.assertNorm(residuals, residuals_py, "Incorrect residuals")



    def test_lmrob(self):
        stackloss = self.stackloss
        dictStackloss = self.dictStackloss
        data = pandas.DataFrame(dictStackloss)

        # R's lmrob
        for method in ["S", "MM", "SMD", "SMDM"]:
            m0_R = self.robustbase.lmrob("stack.loss ~ Water.Temp + Air.Flow + Acid.Conc.", data=stackloss, method=method)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            m0_py = lmrob("stack_loss ~ WaterTemp + AirFlow + AcidConc", data=data, method=method)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            


            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method)

    def test_lmrob_offset(self):
        stackloss = self.stackloss
        dictStackloss = self.dictStackloss
        data = pandas.DataFrame(dictStackloss)

        # R's lmrob
        for method in ["S", "MM", "SMD", "SMDM"]:
            m0_R = self.robustbase.lmrob("stack.loss ~ Water.Temp + Air.Flow + offset(Acid.Conc.)", data=stackloss, method=method)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            m0_py = lmrob("stack_loss ~ WaterTemp + AirFlow", data=data, offset="AcidConc", method=method)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            
            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)

            # print("coefficients R:\n", coefficients)
            # print("coefficients py:\n", coefficients_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method)

    def test_lmrob_weights(self):
        stackloss = self.stackloss
        dictStackloss = self.dictStackloss
        data = pandas.DataFrame(dictStackloss)

        # R's lmrob
        for method in ["S", "MM", "SMD", "SMDM"]:
            m0_R = self.robustbase.lmrob("stack.loss ~ Water.Temp + Air.Flow", data=stackloss, weights=stackloss.rx2("Acid.Conc."), method=method)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            m0_py = lmrob("stack_loss ~ WaterTemp + AirFlow", data=data, weights="AcidConc", method=method)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            
            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)

            # print("coefficients R:\n", coefficients)
            # print("coefficients py:\n", coefficients_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method)

    def test_lmrob_offset_weights(self):
        stackloss = self.stackloss
        dictStackloss = self.dictStackloss
        data = pandas.DataFrame(dictStackloss)

        # R's lmrob
        for method in ["S", "MM", "SMD", "SMDM"]:
            m0_R = self.robustbase.lmrob("stack.loss ~ Water.Temp", data=stackloss, offset=stackloss.rx2("Air.Flow"), weights=stackloss.rx2("Acid.Conc."), method=method)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            m0_py = lmrob("stack_loss ~ WaterTemp", data=data, offset="AirFlow", weights="AcidConc", method=method)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            
            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)

            # print("coefficients R:\n", coefficients)
            # print("coefficients py:\n", coefficients_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method)



    def test_lmrob_categorical(self):
        categorical_data = self.categorical_data
        categorical_data_py = pandas.DataFrame(self.categorical_data)
        # dataFrame = patsy.demo_data("y", "x1", "x2", "x3", "cats", "dogs", nlevels=3, min_rows=10)
        # dataFrame = pandas.DataFrame(dataFrame)

        for method in ["S", "MM", "SMD", "SMDM"]:
            # R's lmrob
            control = self.robustbase.lmrob_control(psi='optimal', tuning_chi=20.0, bb = 0.0003846154, tuning_psi=20.0, method=method , cov='.vcov.w')
            m0_R = self.robustbase.lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data, control=control)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            control_py = LmrobControl(psi="optimal", tuning_chi = 20, bb = 0.0003846154, tuning_psi=20, method=method, cov="_vcov_w")
            m0_py = lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data_py, control=control_py)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            
            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)

            # print("coefficients R:\n", coefficients)
            # print("coefficients py:\n", coefficients_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method)    

    def test_lmrob_categorical_offset(self):
        categorical_data = self.categorical_data
        categorical_data_py = pandas.DataFrame(self.categorical_data)
        # dataFrame = patsy.demo_data("y", "x1", "x2", "x3", "cats", "dogs", nlevels=3, min_rows=10)
        # dataFrame = pandas.DataFrame(dataFrame)

        for method in ["S", "MM", "SMD", "SMDM"]:
            # R's lmrob
            control = self.robustbase.lmrob_control(psi='optimal', tuning_chi=20.0, bb = 0.0003846154, tuning_psi=20.0, method=method , cov='.vcov.w')
            m0_R = self.robustbase.lmrob("y ~ x1*x2 + x3 + x4 + x5 + offset(os)", data=categorical_data, control=control)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            control_py = LmrobControl(psi="optimal", tuning_chi = 20, bb = 0.0003846154, tuning_psi=20, method=method, cov="_vcov_w")
            m0_py = lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data_py, offset="os", control=control_py)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            
            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)

            # print("coefficients R:\n", coefficients)
            # print("coefficients py:\n", coefficients_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method)

    def test_lmrob_categorical_weights(self):
        categorical_data = self.categorical_data
        categorical_data_py = pandas.DataFrame(self.categorical_data)
        # dataFrame = patsy.demo_data("y", "x1", "x2", "x3", "cats", "dogs", nlevels=3, min_rows=10)
        # dataFrame = pandas.DataFrame(dataFrame)

        for method in ["S", "MM", "SMD", "SMDM"]:
            # R's lmrob
            control = self.robustbase.lmrob_control(psi='optimal', tuning_chi=20.0, bb = 0.0003846154, tuning_psi=20.0, method=method , cov='.vcov.w')
            m0_R = self.robustbase.lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data, weights=categorical_data["weights"], control=control)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            control_py = LmrobControl(psi="optimal", tuning_chi = 20, bb = 0.0003846154, tuning_psi=20, method=method, cov="_vcov_w")
            m0_py = lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data_py, weights="weights", control=control_py)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            
            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)

            # print("coefficients R:\n", coefficients)
            # print("coefficients py:\n", coefficients_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method)

    def test_lmrob_categorical_offset_weights(self):
        categorical_data = self.categorical_data
        categorical_data_py = pandas.DataFrame(self.categorical_data)
        # dataFrame = patsy.demo_data("y", "x1", "x2", "x3", "cats", "dogs", nlevels=3, min_rows=10)
        # dataFrame = pandas.DataFrame(dataFrame)

        for method in ["S", "MM", "SMD", "SMDM"]:
            # R's lmrob
            control = self.robustbase.lmrob_control(psi='optimal', tuning_chi=20.0, bb = 0.0003846154, tuning_psi=20.0, method=method , cov='.vcov.w')
            m0_R = self.robustbase.lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data, offset=categorical_data["os"], weights=categorical_data["weights"], control=control)
            coefficients = rpyn.ri2py(m0_R.rx2("coefficients"))
            scale = rpyn.ri2py(m0_R.rx2("scale"))
            residuals = rpyn.ri2py(m0_R.rx2("residuals"))
            
            # Python's lmrob
            control_py = LmrobControl(psi="optimal", tuning_chi = 20, bb = 0.0003846154, tuning_psi=20, method=method, cov="_vcov_w")
            m0_py = lmrob("y ~ x1*x2 + x3 + x4 + x5", data=categorical_data_py, offset="os", weights="weights", control=control_py)
            coefficients_py = m0_py.get("coefficients")
            residuals_py = m0_py.get("residuals")
            scale_py = m0_py.get("scale")
            
            # print("scale R:\n", scale)
            # print("scale py:\n", scale_py)
            
            # print("residuals R:\n", residuals)
            # print("residual py:\n", residuals_py)

            # print("coefficients R:\n", coefficients)
            # print("coefficients py:\n", coefficients_py)
            
            self.assertNorm(scale, scale_py, "Incorrect scale at method %s" % method)
            self.assertNorm(coefficients, coefficients_py, "Incorrect coefficients at method %s" % method)
            self.assertNorm(residuals, residuals_py, "Incorrect residuals at method %s" % method)



    def test_lmrob_tau(self):
        X = self.x
        Y = self.y
        control = self.robustbase.lmrob_control()
        init = self.robustbase.lmrob_S(self.x, self.y, control)
        tau = self.robustbase.lmrob_tau(init, X, control)
        control_py = LmrobControl()
        init_py = lmrob_S(rpyn.ri2py(X), rpyn.ri2py(Y),control_py)
        tau_py = lmrob_tau(init_py, rpyn.ri2py(X), control_py)

        self.assertNorm(rpyn.ri2py(tau), tau_py, "Incorrect tau")

    def test_lmrob_lar(self):
        X = self.x
        Y = self.y
        control = self.robustbase.lmrob_control()
        init = self.robustbase.lmrob_lar(self.x, self.y, control)
        
        coefficients = init.rx2("coefficients")
        scale = init.rx2("scale")
        residuals = init.rx2("residuals")
        k = init.rx2("iter")
        status = init.rx2("status")
        converged = init.rx2("converged")

        control_py = LmrobControl()
        init_py = lmrob_lar(rpyn.ri2py(X), rpyn.ri2py(Y), control_py)

        coefficients_py = init_py.get("coefficients")
        scale_py = init_py.get("scale")
        residuals_py = init_py.get("residuals")
        iter_py = init_py.get("iter")
        status_py = init_py.get("status")
        converged_py = init_py.get("converged")

        self.assertNorm(rpyn.ri2py(coefficients), coefficients_py, "Incorrect coefficients")
        self.assertNorm(rpyn.ri2py(scale), scale_py, "Incorrect scale")
        self.assertNorm(rpyn.ri2py(residuals), residuals_py, "Incorrect residuals")
        self.assertNorm(rpyn.ri2py(k), iter_py, "Incorrect iter")
        self.assertNorm(rpyn.ri2py(status), status_py, "Incorrect status")
        self.assertEqual(rpyn.ri2py(converged), converged_py, "Incorrect converged")

    def test_lmrob_D_fit(self):
        X = self.x
        Y = self.y
        control = self.robustbase.lmrob_control()
        init = self.robustbase.lmrob_S(self.x, self.y, control)
        
        m0 = self.robustbase.lmrob__D__fit(init, X, control)
        converged = m0.rx2("converged")
        scale = m0.rx2("scale")
        w_ = rpyn.ri2py(m0.rx2("rweights"))
        r_ = rpyn.ri2py(m0.rx2("residuals"))
        kappa = rpyn.ri2py(m0.rx2("kappa"))
        tau_ = rpyn.ri2py(m0.rx2("tau"))
        control_py = LmrobControl()
        init_py = lmrob_S(rpyn.ri2py(X), rpyn.ri2py(Y),control_py)
        m0_py = lmrob__D__fit(init_py, rpyn.ri2py(X), control_py)
        converged_py = m0_py.get("converged")
        scale_py = m0_py.get("scale")     

        self.assertEqual(rpyn.ri2py(converged), np.array([converged_py]), "Incorrect Converged")
        self.assertNorm(rpyn.ri2py(scale), np.array([scale_py]), "Incorrect Scale")

        return

    def tearDown(self):
        pandas2ri.deactivate()
        del(self.stats)
        del(self.robustbase)
        del(self.x)
        del(self.y)



class LmrobUtilities(unittest.TestCase):

    def setUp(self):
        self.robustbase = importr('robustbase')
        self.psiList = ['bisquare']#,'welsh','huber',  'optimal','hampel','ggw','lqq']

        pandas2ri.activate()
    

    def assertNorm(self, m1, m2, messages):
        self.assertEqual(m1.size, m2.size, "Incorrect Size")
        size = m1.size
        error = np.linalg.norm(m1 - m2) / size
        tol = 1e-10
        self.assertLess(error, tol, messages)


    def test_ghq(self):
        for i in range(1,20):
            vd = self.robustbase.ghq(i)
            nodes = rpyn.ri2py(vd.rx2("nodes"))
            weights = rpyn.ri2py(vd.rx2("weights"))
            vd_py = ghq(i)
            nodes_py = vd_py.get("nodes")
            weights_py = vd_py.get("weights")
            self.assertNorm(nodes, nodes_py, "Icorrect nodes at %d and %d" % (i,  i % 2))
            self.assertNorm(weights, weights_py, "Icorrect weights at %d and %d" % (i,  i % 2))


    def test_mpsi(self):
        x = np.random.rand(10)
        robjects.globalenv['dataframe'] = x
        control = self.robustbase.lmrob_control()   
        for psi in self.psiList:
            res = self.robustbase.Mpsi(x, control.rx2("tuning.psi"), control.rx2("psi"), deriv=0)
            res_py = Mpsi(x, rpyn.ri2py(control.rx2("tuning.psi")), psi, deriv=0)
            self.assertNorm(rpyn.ri2py(res), res_py, "Incorrect Operation for %s" % psi)


    def test_mchi(self):
        x = np.random.rand(10)
        robjects.globalenv['dataframe'] = x
        control = self.robustbase.lmrob_control()   
        for psi in self.psiList:
            res = self.robustbase.Mchi(x, control.rx2("tuning.chi"), control.rx2("psi"), deriv=0)
            res_py = Mchi(x, rpyn.ri2py(control.rx2("tuning.chi")), psi, deriv=0)
            self.assertNorm(rpyn.ri2py(res), res_py, "Incorrect Operation for %s" % psi)


    def test_mwgt(self):

        x = np.random.rand(10)
        robjects.globalenv['dataframe'] = x
        control = self.robustbase.lmrob_control()   
        for psi in self.psiList:
            res = self.robustbase.Mwgt(x, control.rx2("tuning.chi"), control.rx2("psi"))
            res_py = Mwgt(x, rpyn.ri2py(control.rx2("tuning.chi")), psi)
            self.assertNorm(rpyn.ri2py(res), res_py, "Incorrect Operation for %s" % psi)

    def test_lmrob_rweights(self):
        x = np.random.randn(10)
        scale = np.random.rand()
        robjects.globalenv['dataframe'] = x
        for psi in self.psiList:
            control = self.robustbase.lmrob_control()
            rw = self.robustbase.lmrob_rweights(x,
                                                scale,
                                                control.rx2("tuning.psi"),
                                                psi)
            rw_py = lmrob_rweights( x,
                                    scale,
                                    rpyn.ri2py(control.rx2("tuning.psi")),
                                    psi)
            self.assertNorm(rpyn.ri2py(rw), rw_py, "Incorrect Weights for %s" % psi)


    def test_lmrob_hat(self):
        x = np.random.randn(10,2)
        robjects.globalenv['dataframe'] = x
        h = self.robustbase._lmrob_hat(x)
        h_py = _lmrob_hat(x)
        self.assertNorm(rpyn.ri2py(h), h_py, "Incorrect Result")


    def test_lmrob_kappa(self):
        for method in ["S", "MM"]:
            control = self.robustbase.lmrob_control(method=method)
            control_py = LmrobControl(method=method)
            res = self.robustbase.lmrob_kappa(control=control)
            res_py = lmrob_kappa(control=control_py)
            self.assertNorm(rpyn.ri2py(res), res_py, "Incorrect residual for the method %s" % method)
    
    def tearDown(self):
        pandas2ri.deactivate()

if __name__ == '__main__':    
    unittest.main(verbose=2)
    


