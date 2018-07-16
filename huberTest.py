# -*- coding: utf-8 -*-


import imp
import numpy

huberM = imp.load_source('huberM', 'C:/Users/Alberto/Desktop/Tesi/huberM.py')
huberM.huberM()
aa = numpy.arange(1, 10)
aa = numpy.append(aa,1000)
huberM.huberM(aa)
w_ = [1]*10
w_[3] = 3
w_[8] = 2
huberM.huberM(aa,weights = w_)

