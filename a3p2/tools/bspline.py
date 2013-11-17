#===========================================================================
# Imports

import logging

#===========================================================================
# Functions & classes

#---------------------------------------------------------------------------
class BSpline1D :
    """
    1D b-spline object.

    WARNING: The data (knots, weights) are not checked for consistency.

    For more information in b-splines see:
    http://en.wikipedia.org/wiki/B-spline
    http://mathworld.wolfram.com/B-Spline.html
    """
    
    knots, weights = [], []

    def __init__(self, knots, weights) :
        self.knots, self.weights = knots, weights

    def eval_base_spline(self, t, j, n) :

        if j + n + 1 > len(self.knots) - 1 :
            logging.warning("Request knot is out of bounds => returning 0.")
            return 0.

        if n == 0 :
            if (t >= self.knots[j]) and (t < self.knots[j + 1]) :
                return 1.
            else :
                return 0.

        result = 0.

        if self.knots[j + n] - self.knots[j] > 0. :
            result = (t - self.knots[j]) * self.eval_base_spline(t, j, n - 1) / (self.knots[j + n] - self.knots[j])
        if self.knots[j + n + 1] - self.knots[j + 1] > 0. :
            result += (self.knots[j + n + 1] - t) * self.eval_base_spline(t, j + 1, n - 1) / (self.knots[j + n + 1] - self.knots[j + 1])
        
        return result

    def eval(self, t, n) :

        if t > max(self.knots) or t < min(self.knots) :
            return 0.

        val = 0.
        for i in range(len(self.knots) - n - 1) :
            val += self.weights[i] * self.eval_base_spline(t, i, n)

        return val

#---------------------------------------------------------------------------
def data_to_bspline1d(data_x, data_y, order=2) :
    """
    Returns a BSpline1D object matching the data set.
    """

    if order < 2 or order > 3 :
        raise ValueError("Only spline of order 2 or 3 are supported")

    if len(data_x) < 2 :
        raise ValueError("Data lists need to have at least two entries.")

    if len(data_x) != len(data_y) :
        raise ValueError("Data lists (x,y) need to have the same length.")

    k, w = [], []

    for i in range(order + 1) :
        k.append(data_x[0])
        w.append(data_y[0])
    w.pop()
    w.pop()
    for i in range(len(data_x) - 1) :
        if order % 2:
            k.append(data_x[i])
        else :
            k.append((data_x[i] + data_x[i+1])/2.)
    w.extend(data_y)
    for i in range(order + 1) :
        k.append(data_x[len(data_x) - 1])
        w.append(data_y[len(data_y) - 1])

    return BSpline1D(k, w)

#===========================================================================
