
% Joe Hollowed
% PHY342
%
% This function preforms the cubic spline interpolation algorithm on a dataset. Some of this code
% is taken from the text (Devries & Hasbun 2011)
%
% Last edited 6/24/16

function [x_interp, p] = cubic_spline(x, f, dfdx1 = 'none', dfdxN = 'none', resolution = 100)

	% param x: the x-values of the data
	% param f: the function values of the data
	% param resolution: desired function resolution increase. 
	%    		    For example, if x was N values long, the 
	%		    interpolation will contain N*resolution 
	%		    values. Default is 100. 

	% find second derivative of function at all x values using gaussian elimination
	second = splineInit(x, f, dfdx1, dfdxN);
	
	%Now all that's left to do is preform eq. 
	N = length(x);
	N_interp = N*resolution;
	x_interp = linspace(x(1), x(N), N_interp);
	p = zeros(1, length(x_interp));

	for i = 1:N_interp
		xval = x_interp(i);

		%find which two original data points the current x value lies between
		j = 1;
		while(xval > x(j+1))
			j = j + 1;
		end

		%now calculate p(x)
		hj = (x(j+1) - x(j));
		p(i) = f(j) + ( (f(j+1) - f(j))/hj - (hj*second(j+1))/6 - (hj*second(j))/3 ) *...
		       (xval - x(j)) + (second(j)/2)*((xval - x(j))^2) + ...
		       ( (second(j+1) - second(j))/(6*hj) ) * ((xval - x(j))^3);
	end
end
		       
	
		
