% Joe Hollowed
% PHY 342
%
% Function to find root of some function f(x), given initial root bounds
%
% Last Edited 9/18/16

function [] = rootSolver3(func, xl, xr, tols)

	% param func: function to pass to root finder, as a string 
	% param xl: left endpoint of secant line
	% param xr: right endpoint of secant line
	% return: None

	%----------Check validity of function----------
 	x = linspace(1,10,10);
	try
		test = eval(func);
	catch
		error('Function f(x) uses invalid operators');
	end
	
	roots = zeros(1, length(tols));
	iters = zeros(1, length(tols));

	for n = 1:length(tols)
		%----------Preform root finding----------
		[root, iter] = secant(func, xl, xr,  tols(n));
		roots(n) = root;
		iters(n) = iter;
	end
	disp(roots)
	disp(iters)
	
	
	if(length(tols) > 1)
		figure;
		semilogx(1./tols, iters, 'o-', 'linewidth', 2, 'Color', 'blue');
		xlabel('1 / tolerance', 'fontsize', 16);
		ylabel('iterations', 'fontsize', 16);
		title(sprintf('Tolerance vs. Iterations\nSecant'), 'fontsize', 18);
	end
end	
