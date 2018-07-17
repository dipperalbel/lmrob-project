# -*- coding: utf-8 -*-

import numpy
import warnings

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
  return (mu, s, it, se2)

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

# function for mad (median absolute deviation)
def mad(data, center = None, constant = 1, axis = None):
  "mad function mimicing R's mad function"
  # if no center given, take median)
  if(center == None):
    center = numpy.median(data, axis = None)  
  # ! to mimic R mad function constant must be put to 1.4826 
  # (I dont't see a obvious sense of the constant)
  return constant*numpy.median(numpy.absolute(numpy.array(data) - center), axis)


# helper function sum
def sumU(x, weights):
  "Calculate weighted sum"
  if(len(x) != len(weights)):
    raise ValueError('x and weights not of same length')
  sums = numpy.sum(numpy.array(x) * numpy.array(weights))  
  return(sums)

# high weighted median 
def wgtHighMedian(values, sample_weight=None, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param sample_weight: array-like of the same length as `array`
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    quantiles=[0.5]
    values = numpy.array(values)
    quantiles = numpy.array(quantiles)
    if sample_weight is None:
        sample_weight = numpy.ones(len(values))
    sample_weight = numpy.array(sample_weight)
    assert numpy.all(quantiles >= 0) and numpy.all(quantiles <= 1), 'quantiles should be in [0, 1]'
    sorter = numpy.argsort(values)
    values = values[sorter]
    sample_weight = sample_weight[sorter]
    weighted_quantiles = numpy.cumsum(sample_weight) - 0.5 * sample_weight
  
    weighted_quantiles /= numpy.sum(sample_weight)
    res = numpy.interp(quantiles, weighted_quantiles, values)
    # find closest value
    abs_diff = numpy.absolute(numpy.array(values-res))
    abs_idx = numpy.where(abs_diff == numpy.min(abs_diff))
    fin = numpy.max(values[abs_idx])
    return(fin)