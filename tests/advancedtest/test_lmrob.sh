#!/bin/bash

echo "------------------------------------------------------------------------------"
echo "                              Performing R calls                              "
echo "------------------------------------------------------------------------------"

Rscript lmrob-methods.R

echo "------------------------------------------------------------------------------"
echo "                           Performing Python calls                            "
echo "------------------------------------------------------------------------------"

python lmrob-methods.py


echo "METHOD: S"
echo "+----------------+---------------------+-------------+"
echo "| R Coefficients | Python Coefficients |    Diff     |"
echo "+----------------+---------------------+-------------+"
python  print_table.py coeff_Smethod_base.csv coeff_Smethod_new.csv
echo "+----------------+---------------------+-------------+"
echo "|  R residuals   |  Python residuals   |    Diff     |"
echo "+----------------+---------------------+-------------+"
python  print_table.py res_Smethod_base.csv res_Smethod_new.csv
echo "+----------------+---------------------+-------------+"
echo ""
echo "METHOD: MM"
echo "+----------------+---------------------+-------------+"
echo "| R Coefficients | Python Coefficients |    Diff     |"
echo "+----------------+---------------------+-------------+"
python  print_table.py coeff_MMmethod_base.csv coeff_MMmethod_new.csv
echo "+----------------+---------------------+-------------+"
echo "|  R residuals   |  Python residuals   |    Diff     |"
echo "+----------------+---------------------+-------------+"
python  print_table.py res_MMmethod_base.csv res_MMmethod_new.csv
echo "+----------------+---------------------+-------------+"
echo ""
echo "METHOD: SMD"
echo "+----------------+---------------------+-------------+"
echo "| R Coefficients | Python Coefficients |    Diff     |"
echo "+----------------+---------------------+-------------+"
python  print_table.py coeff_SMDmethod_base.csv coeff_SMDmethod_new.csv
echo "+----------------+---------------------+-------------+"
echo "|  R residuals   |  Python residuals   |    Diff     |"
echo "+----------------+---------------------+-------------+"
python  print_table.py res_SMDmethod_base.csv res_SMDmethod_new.csv
echo "+----------------+---------------------+-------------+"
echo ""
echo "METHOD: SMDM"
echo "+----------------+---------------------+-------------+"
echo "| R Coefficients | Python Coefficients |    Diff     |"
echo "+----------------+---------------------+-------------+"
python  print_table.py coeff_SMDMmethod_base.csv coeff_SMDMmethod_new.csv
echo "+----------------+---------------------+-------------+"
echo "|  R residuals   |  Python residuals   |    Diff     |"
echo "+----------------+---------------------+-------------+"
python  print_table.py res_SMDMmethod_base.csv res_SMDMmethod_new.csv
echo "+----------------+---------------------+-------------+"