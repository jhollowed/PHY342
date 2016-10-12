% Joe Hollowed
% PHY342
%
% Function to execute the secant root finding algorithm
%
% Last edited 9/18/16

function [root, iter] = secant_mod(xl, xr, tol, xall, fall, dfdx1='none', dfdxN='none', i=0, err=0)

	% param func: function as a valid matlab string
	% param xl: left endpoint of secant line
	% param xr: right endpoint of secant line
	% param tol: desired tolerance
	% return: [root, iterations to find root]
	
	x = [xl, xr];
	f = zeros(1,2);

	f(1) = cubic_spline_func(xl, xall, fall, dfdx1, dfdxN);
	f(2) = cubic_spline_func(xr, xall, fall, dfdx1, dfdxN);
	dfdx = diff(f) / (xr-xl);
	x_new = xr - (f(2) / dfdx);
	
	
	if(err <= tol && i > 0)
		root = x_new;
		iter = i;
		return;
	end

	err = relativeErr(xr, x_new);
	i = i + 1;
	[root, iter] = secant_mod(xr, x_new, tol, xall, fall, dfdx1, dfdxN, i, err);
	return;
