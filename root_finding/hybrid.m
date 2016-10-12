% Joe Hollowed
% PHY342
%
% Function to execute the newton-raphson root finding algorithm
%
% Last edited 9/12/16

function [root, iter] = hybrid(func, xi, xl, xr, tol, i=0, err=0)

	% param func: function as a valid matlab string
	% param a: point to evaluate derivative
	% param xd: desired tolerance
	% param makePlot: whether or not to animate the algorithm
	% return: [root, iterations to find root]
	
	dbstop;
	
	dx = 0.0001;	
	x = [xi, xi+dx];
	f = eval(func);
	dfdx = diff(f) / 0.0001;
	
	x_new = xi - (eval(strrep(func,'x','xi')) / dfdx);
	
	if(x_new >= xl && x_new <= xr)
		
		if(err <= tol && i > 0)
			root = x_new;
			iter = i;
			return;
		end
	
		err = relativeErr(xi, x_new);
		i = i + 1;
		[root, iter] = hybrid(func, x_new, tol, i, err);
		return;
	else
		xc = xl + (abs(xr-xl))/2;
		if(err <= tol && i > 0)
			return;
		end
		
		%----------if root not found, find function values at root bounds----------
		funcl = strrep(func, 'x', 'xl');
		funcr = strrep(func, 'x', 'xr');
		funcc = strrep(func, 'x', 'xc');
		fl = eval(funcl);
		fr = eval(funcr);
		fc = eval(funcc);

		%----------recursive call to bisect----------
		if(fl * fc < 0)
			i = i+1;
			err = relativeErr(xc, xl);
			[xc, i] = hybrid(func, xl, xl, xc, tol, i, err);
			root = xc;
			iter = i;
			return;
		elseif(fc * fr < 0)
			i = i+1;
			err = relativeErr(xc, xr);
			[xc, i] = hybrid(func, xc, xc, xr, tol, i, err);
			root = xc;
			iter = i;
			return;
		end
	end
end
				
		
