# -*- coding: utf-8 -*-

import numpy
import warnings

# high weighted median


def wgtHighMedian(x, weights=None):
    """Compute the weighted Hi-Median of x.
    :param x: numeric vector
    :param weights: numeric vector of weights; of the same length as x.
    """
    # make numpy array
    # numpy.array(x) create an array from x and it assaign to  the variable x
    x = numpy.array(x)

    # sort x
    # after the array is created, it is sorted
    x_sort = numpy.argsort(x)
    x = x[x_sort]
    # parse weights and sort
    if (weights != None):
        if(len(x) != len(weights)):
            raise ValueError('x and weights cannot have different length')
    if(weights is None):
        weights = numpy.array([1]*len(x))
    weights = numpy.array(weights)
    weights = weights[x_sort]

    # cumSum and interpolate
    # cumsum return an array with the cumulative sum of the elements. The array obtain it is substracted with  0.5 *  the array weights
    cumW = numpy.cumsum(weights) - 0.5 * weights
    cumW /= numpy.sum(weights)

    # function interp
    res = numpy.interp(0.5, cumW, x)

    # find closest value
    abs_diff = numpy.absolute(numpy.array(x-res))
    abs_idx = numpy.where(abs_diff == numpy.min(abs_diff))
    fin = numpy.max(x[abs_idx])
    return(fin)

# test
# wgtHighMedian(numpy.array([0,3,6,9]),numpy.array([1,3,3,1]))
# wgtHighMedian(numpy.array([0,3,3,3,6,6,6,9]),None)
# test cases wgt_himedian
#a = numpy.array([1, 3, 5, 7])
#b = numpy.array([1, 1, 0.5, 1])
#a * b
# wgtHighMedian(a)
#wgtHighMedian(a, b)


# function for mad (median absolute deviation)
def mad(x, center=None, constant=1.4826, na_rm=False, low=False, high=False):
    """Compute the median absolute deviation, i.e., the (lo-/hi-) median of the absolute deviations from the median, and (by default) adjust by a factor for asymptotically normal consistency.
    :param x: numeric vector.
    :param center: optionally, the centre: defaults to the median.
    :param constant: scale factor.
    :param na_rm: if True then NaN values are stripped from x before computation takes place.
    :param low: if True, compute the ‘lo-median’, i.e., for even sample size, do not average the two middle values, but take the smaller one.
    :param high: if True, compute the ‘hi-median’, i.e., take the larger of the two middle values for even sample size.
    """
    axis = None
    if (na_rm == True):
        i = numpy.isnan(x)
        x = numpy.array(x)[numpy.where(i == False)]
    if(center == None):
        center = numpy.median(x, axis=None)
    # Fin qua tutto apposto
    n = len(x)
    if ((low or high) and n % 2 == 0):
        if (low and high):
            raise ValueError('low and high cannot be both true')
        n2 = ((n//2) + high) - 1
        # return constant * (numpy.partition((numpy.absolute(numpy.array(x) - center)), n2)))[n2]
        return constant * numpy.partition((numpy.absolute(numpy.array(x) - center)), n2)[n2]
    return constant*numpy.median(numpy.absolute(numpy.array(x) - center), axis)

# test mad
#mad([-1, -0.5, 0, 1, 2], center = 3)
# mad([-1, -0.5, 0, 1, 2], center = 3, constant = 1.4826) # R default

# function to determine correction factor Tau for Huber-M-estimator


def tauHuber(x, mu, k=1.5, s=None, resid=None):
    """Calculate correction factor Tau for the variance of Huber-M-Estimators.
    :param x: numeric vector.
    :param mu: location estimator.
    :param k: tuning parameter of Huber Psi-function.
    :param s: scale estimator held constant through the iterations.
    :param resid: the value of (x - mu)/s .
    """
    # get needed parameters
    if(s == None):
        #s = mad(x)
        # mimic R default parameters: set constant to 1.4826
        s = mad(x, constant=1.4826)
    if(resid == None):
        resid = (numpy.array(x) - mu)/s
    # perform calculation
    inta = numpy.absolute(resid) <= k
    psi = numpy.copy(resid)
    psi[inta == False] = numpy.sign(resid[inta == False])*k
    psiP = numpy.bool_(inta)
    res = len(x) * numpy.sum(numpy.square(psi)) / numpy.square(numpy.sum(psiP))
    # return
    return(res)

# test tauHuber
#tauHuber([-1, -0.5, 0, 1, 2], mu = 3, resid = None)

# helper function sum


def sumU(x, weights):
    """Calculate weighted sum
    :param x: numeric vector.
    :param weights: numeric vector of weights; of the same length as x.
    """
    if(len(x) != len(weights)):
        raise ValueError('x and weights not of same length')
    sums = numpy.sum(numpy.array(x) * numpy.array(weights))
    return(sums)

# test sumU
# load a and b above
#sumU(a, b)
# sumU(a, [1,2,3]) # cause error

# huberM function


def huberM(x,
           k=1.5,
           weights=None,
           tol=0.000001,
           mu=None,
           s=None,
           se=False,
           warn0scale=True):
    """(Generalized) Huber M-estimator of location with MAD scale.
    :param x: numeric vector.
    :param k: positive factor; the algorithm winsorizes at k standard deviations.
    :param weights: numeric vector of non-negative weights of same length as x, or None.
    :param tol: convergence tolerance.
    :param mu: initial location estimator.
    :param s: scale estimator held constant through the iterations.
    :param se: logical indicating if the standard error should be computed and returned (as SE component). Currently only available when weights is None.
    :param warn0scale: logical; if True, and s is 0 and len(x) > 1, this will be warned about.
    """
    # parse x
    x = numpy.array(x)
    # parse mu
    if(mu == None):
        if(weights == None):
            mu = numpy.median(x)
        else:
            mu = wgtHighMedian(x, weights)
    # parse s
    if(s == None):
        if(weights == None):
            #s = mad(x, center = mu)
            # mimic R default parameters: set constant to 1.4826
            s = mad(x, center=mu, constant=1.4826)
        else:
            s = wgtHighMedian(numpy.abs(x - mu), weights)

    # get only valid data
    if(numpy.sum(numpy.isnan(x)) > 0):
        i = numpy.isnan(x)
        x = numpy.array(x)[numpy.where(i == False)]
        if(weights != None):
            weights = numpy.array(weights)[numpy.where(i == False)]
    # get length
    n = len(x)
    # check weights
    if(weights == None):
        sumW = numpy.copy(n)
    else:
        if(numpy.sum(numpy.array(weights) <= 0) or len(weights) != n):
            raise ValueError('Something wrong with weights')
        sumW = numpy.sum(weights)
    it = 0
    # if sum of weights 0
    if(sumW == 0):
        # return(list(mu = NA., s = NA., it = it, se = NA.)) #R
        return((numpy.NaN, numpy.NaN, it, numpy.NaN))
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

            y = numpy.minimum(numpy.maximum(mu - k * s, x), mu + k * s)
            if(weights == None):
                mu1 = numpy.sum(y)/sumW
            else:
                mu1 = sumU(y, weights)/sumW
            if(numpy.absolute(mu - mu1) < (tol * s)):
                break
            mu = numpy.copy(mu1)

    if(se):
        se2 = s * numpy.sqrt(tauHuber(x, mu, k, s)/n)
    else:
        se2 = numpy.NaN
    # return
    mu2 = mu.item(0)
    return((mu2, s, it, se2))

# test huberM
#aa = numpy.arange(1, 10)
#aa = numpy.append(aa,1000)
# huberM(aa)
#w_ = [1]*10
#w_[3] = 3
#w_[8] = 2
#huberM(aa,weights = w_)
#huberM([11,20,14,15], 4, [2,2,4, 5], 0.00002, 20)
#huberM(x = [11,20,14,15], k = 4, tol = 0.2, se = True)
