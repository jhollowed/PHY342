% Joe Hollowed
% PHY342
% 
% Function to find coefficients fo a system of nonlinear equation 
% via crude variation of each coefficient, until the least squares
% sum, or chi squared, is minimized. 
%
% Last edited 10/2/16

function [a] = minimize(x, y, a, h, Y, method='crude', maxIter = 15)

	% param x: experimental x-values
	% param y: experimental y-values
	% param a: vector of initial guesses for all coefficients
	% param h: vector of initial step sizes for each coefficeint
	% param Y: a function handle of the expected function to fit to
	% param method: method to use to minimize the least squares sum 
	%		(default is 'crude')
	% param maxIter: maximum iterations to run before accepting value 
	%		 and returning
	% return: the final coefficients

	% Dependence: crude.m, newton.m
	
	
	if(method == 'crude')

		LSE = 1;
		LSE_new = 0;
		i = 0;

		while(LSE_new < LSE & i < maxIter)
	
			LSE = sum((y - Y(a, x, y)).^2);
			[a, h] = crude(a, h, x, y, Y);			
			LSE_new = sum((y - Y(a, x, y)).^2);
			i = i+1;
		end
 
	elseif(method == 'newton')
	end
end

