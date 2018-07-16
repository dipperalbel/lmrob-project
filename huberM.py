# -*- coding: utf-8 -*-

import numpy
import warnings

# high weighted median 
def wgtHighMedian(values, sample_weight=None, quantiles=[0.5], values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = numpy.array(values)
    quantiles = numpy.array(quantiles)
    if sample_weight is None:
        sample_weight = numpy.ones(len(values))
    sample_weight = numpy.array(sample_weight)
    assert numpy.all(quantiles >= 0) and numpy.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = numpy.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = numpy.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= numpy.sum(sample_weight)
    res = numpy.interp(quantiles, weighted_quantiles, values)
    # find closest value
    abs_diff = numpy.absolute(numpy.array(values-res))
    abs_idx = numpy.where(abs_diff == numpy.min(abs_diff))
    fin = numpy.max(values[abs_idx])
    return(fin)
    
# test
#wgtHighMedian(numpy.array([0,3,6,9]),numpy.array([1,3,3,1])) 
#wgtHighMedian(numpy.array([0,3,3,3,6,6,6,9]),None)
# test cases wgt_himedian
#a = numpy.array([1, 3, 5, 7])
#b = numpy.array([1, 1, 0.5, 1])
#a * b
#wgtHighMedian(a)
#wgtHighMedian(a, b)


# function for mad (median absolute deviation)
def mad(data, center = None, constant = 1, axis = None):
  "mad function mimicing R's mad function"
  # if no center given, take median)
  if(center == None):
    center = numpy.median(data, axis = None)  
  # ! to mimic R mad function constant must be put to 1.4826 
  # (I dont't see a obvious sense of the constant)
  return constant*numpy.median(numpy.absolute(numpy.array(data) - center), axis)

# test mad
#mad([-1, -0.5, 0, 1, 2], center = 3)
#mad([-1, -0.5, 0, 1, 2], center = 3, constant = 1.4826) # R default

# function to determine correction factor Tau for Huber-M-estimator
def tauHuber(x, mu, k = 1.5, s = None, resid = None):
  "Calculate correction factor Tau for the variance of Huber-M-Estimators"
  # get needed parameters
  if(s == None):
    #s = mad(x)  
    s = mad(x, constant = 1.4826) # mimic R default parameters: set constant to 1.4826
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
  "Calculate weighted sum"
  if(len(x) != len(weights)):
    raise ValueError('x and weights not of same length')
  sums = numpy.sum(numpy.array(x) * numpy.array(weights))  
  return(sums)
# helper function pmax
def pmin(mmax, vec):
  "Parallel minimum, mimicing R's pmin function"
  vecMin = numpy.array(vec)
  vecMin[vecMin > mmax] = mmax
  return(vecMin)
def pmax(mmin, vec):
  "Parallel maximum, mimicing R's pmax function"
  vecMax = numpy.array(vec)
  vecMax[vecMax < mmin] = mmin
  return(vecMax)

# test function
# load a and b above
#sumU(a, b)
#sumU(a, [1,2,3])
#pmin(4,[1,2,5])
#pmax(4,[1,2,5])

# huberM function
def huberM(x, 
          k = 1.5, 
          weights = None, 
          tol = 0.000001, 
          se = False,
		  warn0scale = True):
  # get needed parameters
  if(weights == None):
    mu = numpy.median(x)
    #s = mad(x, center = mu) 
    s = mad(x, center = mu, constant = 1.4826) # mimic R default parameters: set constant to 1.4826
  else:
    mu = wgtHighMedian(x, weights)
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
    #return(list(mu = NA., s = NA., it = it, se = NA.)) #R
    return((numpy.NaN, numpy.NaN, it, numpy.NaN))
  if(se and weights != None):
    raise ValueError('Std.error computation not yet available for the case of weights')
  if(s <= 0):
    if(s < 0): 
      raise ValueError('negative scale s')
    if(warn0scale == True):
      warnings.warn("scale 's' is zero -- returning initial 'mu'")
  else:
    while(1):
      it = it + 1
      # import pmax and pmin here
      y = numpy.minimum(numpy.maximum(mu - k * s, x), mu + k * s)
      if(weights == None):
        mu1 = numpy.sum(y)/sumW
      else:
        mu1 = sumU(y,weights)/sumW
      if(numpy.absolute(mu - mu1) < (tol * s)):
        break
      mu = numpy.copy(mu1)
  
  if(se):
    se2 = s * numpy.sqrt(tauHuber(x, mu, s, k)/n)
  else:
    se2 = numpy.NaN
  # return
  return((mu, s, it, se2))
 
# test huberM
#aa = numpy.arange(1, 10)
#aa = numpy.append(aa,1000)
#huberM(aa)
#w_ = [1]*10
#w_[3] = 3
#w_[8] = 2
#huberM(aa,weights = w_)