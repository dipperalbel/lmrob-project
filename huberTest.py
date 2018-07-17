# -*- coding: utf-8 -*-


import huberM

huberM.huberM([1,2,3,4,5,6,7,8,9,1000])

huberM.huberM([1,2,3,4,5,6,7,8,9,1000], 4)

huberM.huberM([11,20,14,15], 4, [2,2,4,5])

huberM.huberM([11,20,14,15], 4, [2,2,4,5], 0.00002)

huberM.huberM(x = [11,20,14,15], k = 4, tol = 0.00002, se = True)

huberM.tauHuber([1,2,3,4,5] , 1)