% Joe Hollowed
% PHY342
%
% This function solves the general ODE dy/dx
% by repeatedly approximating the solution, y(x), at 
% various x-values, via numerous integration methods. 
%
% Last edited 10/10/16

function [x, y] = odeSolver(y0, x0, xN, steps, ODE, method)
	
	if(method > 3 || method < 1)
		error('method argument must be 1, 2, or 3');
	end
	
	use = {@euler, @modEuler, @rk4};
	x = linspace(x0, xN, steps);
	y = zeros(1, length(x));
	y(1) = y0;
	h = diff(x)(1);

	for i = 1:length(x)-1;
		yn = use{method}(x(i), y(i), h, ODE);
		y(i+1) = yn;
	end
end
