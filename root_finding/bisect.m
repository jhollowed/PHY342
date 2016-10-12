% Joe Hollowed
% PHY342
%
% Function to preform the recursive bisection root finding algorithm
%
% Last edited 9/12/16

function [xc, i] = bisect(func, xl, xr, dx, makePlot=false, handles=0, i=0, err=0)
		
	% param func: a function of x in the form of a string
	% param xl: left bound for the root range
	% param xr: right bound for the root range
	% param dx: the desired tolerance
	% param makePlot: whether or not to animate the algorithm
	% param handles**: figure and plot handles to share across recursive calls.
	% param i**: current iteration to share across reacursive calls.
	% param err**: relative error between last and current iteration
	% return: [root, iterations required to find root]
	% ** denotes optional args that the user probably shouldn't pass		
	
	%----------Check for root----------
	xc = xl + (abs(xr-xl))/2;
	if(err <= dx && i > 0)
		return;
	end
		
	%----------if root not found, find function values at root bounds----------
	funcl = strrep(func, 'x', 'xl');
	funcr = strrep(func, 'x', 'xr');
	funcc = strrep(func, 'x', 'xc');
	fl = eval(funcl);
	fr = eval(funcr);
	fc = eval(funcc);

		
	%-----------Plot bisection of root range if valid handles were passed----------
	if(makePlot == true)
		handles = bisect_plot(func, dx, err, i, handles, xl, xr, xc);
	end
	%----------recursive call to bisect----------
	if(fl * fc < 0)
		i = i+1;
		err = relativeErr(xc, xl);
		[xc, i] = bisect(func, xl, xc, dx, makePlot, handles, i, err);
		return;
	elseif(fc * fr < 0)
		i = i+1;
		err = relativeErr(xc, xr);
		[xc, i] = bisect(func, xc, xr, dx, makePlot, handles, i, err);
		return;
	end
end
