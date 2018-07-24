import numpy
import warnings
"""
A Sphinx example and cheat sheet.

Documenting Python Code
#######################

For documenting Python code, see ``module/__init__.py`` or the code example below:

.. code-block:: python

    class AClass:
        \"\"\"
        Class docstring, with reference to the :mod:`module`, or another class
        :class:`module.AnotherClass` and its function :func:`module.AnotherClass.foo`.
        \"\"\"
    
    class AnotherClass:
        \"\"\"
        Another class' docstring.
        \"\"\"
        
        def foo(arg1, arg2):
            \"\"\"
            A method's docstring with parameters and return value.
            
            Use all the cool Sphinx capabilities in this description, e.g. to give
            usage examples ...
            
            :Example:
    
            >>> another_class.foo('', AClass())        
            
            :param arg1: first argument
            :type arg1: string
            :param arg2: second argument
            :type arg2: :class:`module.AClass`
            :return: something
            :rtype: string
            :raises: TypeError
            \"\"\"
            
            return '' + 1

Functions
#########

Robustbase functions

.. rubric:: Footnotes

"""

def wgtHighMedian(x, weights=None):
    """
    Compute the weighted Hi-Median of x.

    :Example:

    >>> module.wgtHighMedian([1,2,4,5,6,8,10,12])        

    :param x: numeric vector
    :type x: list
    :param weights: numeric vector of weights; of the same length as x.
    :type weights: list
    :return: The weighted Hi-Median of x.
    :rtype: float
    :raises: ValueError
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

# function for mad (median absolute deviation)
def mad(x, center=None, constant=1.4826, na_rm=False, low=False, high=False):
    """
    Compute the median absolute deviation, i.e., the (lo-/hi-) median of the absolute deviations from the median, and (by default) adjust by a factor for asymptotically normal consistency.
    
    :Example:

    >>> module.mad([1,2,4,5,6,8,10,12])

    :param x: numeric vector.
    :type x: list
    :param center: optionally, the centre: defaults to the median.
    :type center: float
    :param constant: scale factor.
    :type constant: float
    :param na_rm: if True then NaN values are stripped from x before computation takes place.
    :type na_rm: bool
    :param low: if True, compute the ‘lo-median’, i.e., for even sample size, do not average the two middle values, but take the smaller one.
    :type low: bool
    :param high: if True, compute the ‘hi-median’, i.e., take the larger of the two middle values for even sample size.
    :type high: bool
    :return: the median absolute deviation
    :rtype: float
    :raises: ValueError
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
        return constant * numpy.partition((numpy.absolute(numpy.array(x) - center)), n2)[n2]
    return constant*numpy.median(numpy.absolute(numpy.array(x) - center), axis)

def tauHuber(x, mu, k=1.5, s=None, resid=None):
    """
    Calculate correction factor Tau for the variance of Huber-M-Estimators.
    
    :Example:

    >>> module.tauHuber([1,2,4,5,6,8,10,12])

    :param x: numeric vector.
    :type x: list
    :param mu: location estimator.
    :type mu: float
    :param k: tuning parameter of Huber Psi-function.
    :type k: float
    :param s: scale estimator held constant through the iterations.
    :type s: float
    :param resid: the value of (x - mu)/s .
    :type resid: float
    :return: The correction factor Tau for the variance of Huber-M-Estimators.
    :rtype: float
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

def sumU(x, weights):
    """
    Calculate weighted sum
    
    :Example:

    >>> module.sumU([1,2,4,5,6,8,10,12], [52,44,82,24,68,42,82,20])
    
    :param x: numeric vector.
    :type x: list
    :param weights: numeric vector of weights; of the same length as x.
    :type weights: list
    :return: The weighted sum.
    :rtype: float
    :raises: ValueError
    """
    if(len(x) != len(weights)):
        raise ValueError('x and weights not of same length')
    sums = numpy.sum(numpy.array(x) * numpy.array(weights))
    return(sums)

def huberM(x,
           k=1.5,
           weights=None,
           tol=0.000001,
           mu=None,
           s=None,
           se=False,
           warn0scale=True):
    """
    (Generalized) Huber M-estimator of location with MAD scale.

    :Example:

    
    >>> module.huberM([1,2,4,5,6,8,10,12])



    :param x: numeric vector.
    :type x: list
    :param k: positive factor; the algorithm winsorizes at k standard deviations.
    :type k: float
    :param weights: numeric vector of non-negative weights of same length as x, or None.
    :type weights: list
    :param tol: convergence tolerance.
    :type tol: float
    :param mu: initial location estimator.
    :type mu: float
    :param s: scale estimator held constant through the iterations.
    :type s: float
    :param se: logical indicating if the standard error should be computed and returned (as SE component). Currently only available when weights is None.
    :type se: bool
    :param warn0scale: logical; if True, and s is 0 and len(x) > 1, this will be warned about.
    :type warn0scale: bool
    :return: The tuple (mu, s , it , se ) containing the location, scale parameters, number of iterations used and the se component.
    :rtype: tuple
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