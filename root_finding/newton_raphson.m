% Joe Hollowed
% PHY342
%
% Function to execute the newton-raphson root finding algorithm
%
% Last edited 9/12/16

function [root, iter] = newton_raphson(func, xi, tol, makePlot=false, i=0, err=0)

	% param func: function as a valid matlab string
	% param a: point to evaluate derivative
	% param xd: desired tolerance
	% param makePlot: whether or not to animate the algorithm
	% return: [root, iterations to find root]
	
	x = [xi, xi+tol];
	f = eval(func);
	dfdx = diff(f) / tol;
	
	x_new = xi - (eval(strrep(func,'x','xi')) / dfdx);
	
	
	if(err <= tol && i > 0)
		root = x_new;
		iter = i;
		return;
	end
	
	err = relativeErr(xi, x_new);
	i = i + 1;
	[root, iter] = newton_raphson(func, x_new, tol, makePlot, i, err);
	return;
