import pandas
import sys
import os
import inspect
from time import time
import unittest
path_of_script = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe() ))[0]))
sys.path.append(os.path.join(path_of_script, "../"))
from nlrob import *


class NlrobTestCase(unittest.TestCase):

    def setUp(self):        
        self.MSE = lambda y1,y2: np.linalg.norm(y1 - y2) ** 2 / y1.size
        self.relative = lambda y1, y2: abs(y1 - y2) / y2

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


    def test_JDEoptim(self):
        lower = pandas.read_csv("Submodule=JDEOPTIM.Data=lower.csv")
        upper = pandas.read_csv("Submodule=JDEOPTIM.Data=upper.csv")
        par = pandas.read_csv("Submodule=JDEOPTIM.Data=par.csv")
        value = pandas.read_csv("Submodule=JDEOPTIM.Data=value.csv")
        lower = lower.set_index("Unnamed: 0")
        def griewank(x):
            return 1 + np.sum(x ** 2)/4000 - np.prod(np.cos(x/np.sqrt(np.arange(x.size) + 1)))
        init = JDEoptim(lower["x"].values, upper["x"].values, griewank, triter = 100,
                        trace=False)
        self.assertLess(init.get("value"), value["x"].values,
            "Value below threshold")


    def test_nlrob_mm(self):
        data = pandas.read_csv("Submodule=NLROBMM.Data=Input.csv")
        residuals = pandas.read_csv("Submodule=NLROBMM.Data=residuals.csv")
        coef = pandas.read_csv("Submodule=NLROBMM.Data=coefficients.csv").set_index("Unnamed: 0")
        formula = "density ~ Asym/(1 + np.exp(( xmid -np.log(conc) )/scal))"
        lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
            index=[0])
        upper = np.array([1])
        Rfit_MM = nlrob_MM(formula, data, lower, upper)
        self.assertLess(self.MSE(Rfit_MM.get("residuals"), 0),
            1.1 * self.MSE(residuals["x"].values, 0), "Residuals below the threshold")
        for pname in lower.keys():
            self.assertLess(self.relative(Rfit_MM["coefficients"][pname], coef["x"][pname]),
            5e-2,  "Coefficients below the threshold %s" % pname)

    def test_nlrob_tau(self):
        data = pandas.read_csv("Submodule=NLROBTAU.Data=Input.csv")
        residuals = pandas.read_csv("Submodule=NLROBTAU.Data=residuals.csv")
        coef = pandas.read_csv("Submodule=NLROBTAU.Data=coefficients.csv").set_index("Unnamed: 0")
        formula = "density ~ Asym/(1 + np.exp(( xmid -np.log(conc) )/scal))"
        lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
            index=[0])
        upper = np.array([1])
        Rfit_tau = nlrob_tau(formula, data, lower, upper)
        self.assertLess(self.MSE(Rfit_tau.get("residuals"), 0),
            1.1 * self.MSE(residuals["x"].values, 0), "Residuals below the threshold")
        for pname in lower.keys():
            self.assertLess(self.relative(Rfit_tau["coefficients"][pname], coef["x"][pname]),
            5e-2,  "Coefficients below the threshold %s" % pname)

    def test_nlrob_cm(self):
        data = pandas.read_csv("Submodule=NLROBCM.Data=Input.csv")
        residuals = pandas.read_csv("Submodule=NLROBCM.Data=residuals.csv")
        coef = pandas.read_csv("Submodule=NLROBCM.Data=coefficients.csv").set_index("Unnamed: 0")
        formula = "density ~ Asym/(1 + np.exp(( xmid -np.log(conc) )/scal))"
        lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
            index=[0])
        
        upper = np.array([1])
        Rfit_CM = nlrob_CM(formula, data, lower, upper)
        self.assertLess(self.MSE(Rfit_CM.get("residuals"), 0),
            1.1 * self.MSE(residuals["x"].values, 0), "Residuals below the threshold")
        for pname in lower.keys():
            self.assertLess(self.relative(Rfit_CM["coefficients"][pname], coef["x"][pname]),
            5e-2,  "Coefficients below the threshold %s" % pname)

    def test_nlrob_mtl(self):
        data = pandas.read_csv("Submodule=NLROBMTL.Data=Input.csv")
        residuals = pandas.read_csv("Submodule=NLROBMTL.Data=residuals.csv")
        coef = pandas.read_csv("Submodule=NLROBMTL.Data=coefficients.csv").set_index("Unnamed: 0")
        formula = "density ~ Asym/(1 + np.exp(( xmid -np.log(conc) )/scal))"
        lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
            index=[0])
        
        upper = np.array([1])
        Rfit_mtl = nlrob_mtl(formula, data, lower, upper)
        self.assertLess(self.MSE(Rfit_mtl.get("residuals"), 0),
            1.1 * self.MSE(residuals["x"].values, 0), "Residuals below the threshold")
        for pname in lower.keys():
            self.assertLess(self.relative(Rfit_mtl["coefficients"][pname], coef["x"][pname]),
            5e-2,  "Coefficients below the threshold %s" % pname)

    def test_nls(self):
        data = pandas.read_csv("Submodule=NLS.Data=Input.csv")
        residuals = pandas.read_csv("Submodule=NLS.Data=residuals.csv")
        coef = pandas.read_csv("Submodule=NLS.Data=coefficients.csv").set_index("Unnamed: 0")
        formula = "density ~ Asym/(1 + np.exp(( xmid -np.log(conc) )/scal))"
        lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], [1, 0.1, 0.4])),
            index=[0])
        
        upper = np.array([1])
        Rfit_nls = nls(formula, data, lower)
        self.assertLess(self.MSE(Rfit_nls.get("residuals"), 0),
            1.1 * self.MSE(residuals["x"].values, 0), "Residuals below the threshold")
        for pname in lower.keys():
            self.assertLess(self.relative(Rfit_nls["coefficients"][pname], coef["x"][pname]),
                5e-2,  "Coefficients below the threshold %s" % pname)


    def test_nlrob(self):
        methods = ["M", "MM", "CM", "tau", "mtl"]
        formula = "density ~ Asym/(1 + np.exp(( xmid -np.log(conc) )/scal))"
        lower = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.zeros(3))),
            index=[0])
        upper = np.array([1])
        data = pandas.read_csv("Function=NLROBInput.csv")
        for method in methods:
            residuals = pandas.read_csv("Function=NLROB%s.Data=residuals.csv" % method)
            coef = pandas.read_csv("Function=NLROB%s.Data=coefficients.csv" % method).set_index("Unnamed: 0")

            if method != "M":
                Rfit = nlrob(formula, data, lower, lower, upper, method=method)
            else:
                start = pandas.DataFrame(data=dict(zip(["Asym", "xmid", "scal"], np.ones(3))),
                                         index=[0])
                Rfit = nlrob(formula, data, start, method=method)
            self.assertLess(self.MSE(Rfit.get("residuals"), 0),
                1.3 * self.MSE(residuals["x"].values, 0), "Residuals below the threshold")
            for pname in lower.keys():
                self.assertLess(self.relative(Rfit["coefficients"][pname], coef["x"][pname]),
                5e-2,  "Coefficients below the threshold %s" % pname)

