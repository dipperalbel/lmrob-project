from lmrob import *

x = np.array([50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73])
y = np.array([0.44,0.47,0.47,0.59,0.66,0.73,0.81,0.88,1.06,1.20,1.35,1.49,1.61,2.12,11.9,12.4,14.2,15.9,18.2,21.2,4.30,2.40,2.70,2.90])
data = {"x" : x,
        "y" : y
}
formula = 'y ~ x'
    # S Method
m0 = lmrob(formula, data=data, method="S")
print("---------S-Residuals---------")
print(m0.get('residuals'))
m1 = lmrob(formula, data=data, method="MM")
print("---------MM-Residuals---------")
print(m1.get('residuals'))
