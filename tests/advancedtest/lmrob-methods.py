import csv
from lmrob import *
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, r
import rpy2.robjects.numpy2ri as rpyn

def main_test():
    stackloss = r("stackloss")
    #pandas2ri.activate()
    #robjects.globalenv['dataframe'] = stackloss
    data = {"AirFlow" : rpyn.ri2py(stackloss.rx2("Air.Flow")),
            "WaterTemp" : rpyn.ri2py(stackloss.rx2("Water.Temp")),
            "AcidConc" : rpyn.ri2py(stackloss.rx2("Acid.Conc.")),
            "stack_loss" : rpyn.ri2py(stackloss.rx2("stack.loss"))
            }
    formula = 'stack_loss ~ AirFlow + WaterTemp + AcidConc'
    # S Method
    m0 = lmrob(formula, data=data, method="S")
    write_csv(m0.get('coefficients'), 'coeff_Smethod_new.csv')
    write_csv(m0.get('residuals'), 'res_Smethod_new.csv')
    
    # MM Method
    m1 = lmrob(formula, data=data, method="MM")
    write_csv(m1.get('coefficients'), 'coeff_MMmethod_new.csv')
    write_csv(m1.get('residuals'), 'res_MMmethod_new.csv')
    
    # SMD Method
    m4 = lmrob(formula, data=data, method="SMD")
    write_csv(m4.get('coefficients'), 'coeff_SMDmethod_new.csv')
    write_csv(m4.get('residuals'), 'res_SMDmethod_new.csv')
    
    # SMDM Method
    m5 = lmrob(formula, data=data, method="SMDM")
    write_csv(m5.get('coefficients'), 'coeff_SMDMmethod_new.csv')
    write_csv(m5.get('residuals'), 'res_SMDMmethod_new.csv')

def write_csv(data, fname):
    with open(fname, 'w', newline='') as csvfile:
        spamwriter = csv.DictWriter(csvfile, delimiter=';',
                             quotechar='|', quoting=csv.QUOTE_MINIMAL,
                             fieldnames=['"x"'])
        spamwriter.writeheader()
        spamwriter.writerows([{'"x"': row[0]} for row in data])

if __name__ == '__main__':
    main_test()