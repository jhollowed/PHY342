% Joe Hollowed
% PHY342
%
% Function to execute the secant root finding algorithm
%
% Last edited 9/18/16

function [root, iter] = secant(func, xl, xr, tol, i=0, err=0)

	% param func: function as a valid matlab string
	% param xl: left endpoint of secant line
	% param xr: right endpoint of secant line
	% param tol: desired tolerance
	% return: [root, iterations to find root]
	
	x = [xl, xr];
	f = eval(func);
	dfdx = diff(f) / (xr-xl);
	
	x_new = xr - (eval(strrep(func,'x','xr')) / dfdx);
	
	
	if(err <= tol && i > 0)
		root = x_new;
		iter = i;
		return;
	end

	err = relativeErr(xr, x_new);
	i = i + 1;
	[root, iter] = secant(func, xr, x_new, tol, i, err);
	return;
