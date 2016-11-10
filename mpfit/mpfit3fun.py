from mpfit.mpfit3 import mpfit

def mpfitfun(func, x, y, err, start_params, full_output=False, **kw):
	"""Fit the used defined function to the data
	Input:
	- func: the function
	- x: x vector
	- y: y vector
	- err: vector with the errors of y
	- start_params: the starting parameters for the fit
	Output:
	- The tuple (params, yfit) with best-fit params and the values of func
	  evaluated at x
	Keywords:
	- full_output: boolean parameter. If True(default is False) then instead of
	 best-fit parameters the mpfit object is returned

	# Example:
	# params,yfit=mpfitfun(function, x, y, err, [0, 10, 1])

	"""
	def myfunct(p, fjac=None, x=None, y=None, err=None):
	    model = func(x, p)
	    status = 0
	    return [status, (y-model)/err]

	fa={'x' : x, 'y' : y, 'err' : err}
	res = mpfit(myfunct, start_params, functkw=fa, **kw)
	yfit = func(x, res.params)
	if full_output:
		return (res, yfit)
	else:
		return (res.params, yfit)
