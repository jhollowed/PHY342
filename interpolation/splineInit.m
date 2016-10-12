% Joe Hollowed
% PHY342
% 
% Function to solve a tridiagonal system of linear equations, in the context of solving for initial
% values (second derivatives) to use in a cubic spline interpolation. This is simply a more
% specialized version of trisolve.m, which uses Gaussian elimination to solve tridiagonal matrices.
% This function was mostly taken from the text (Devries & Hasbun 2011)
%
% Last edited 9/24/16

function [second] = splineInit(x, f, dpdx1='none', dpdxN='none')

	% param x: vector of data x-coordinates
	% param f: vector of data function values
	% param dpdx1: the first derivative of f(x1)
	% param dpdxN: the first derivative of f(xN)
	% return: vector of second deriavatives of f(x) evaluated at each x
	
	% initialize some vars	
	N = length(x);
	beta = zeros(1, N);
	second = zeros(1, N);
	b1 = 2.0 * (x(2) - x(1));
	beta(1) = b1;
	% error check
	if(beta(1) == 0), error('Zero diagonal element in SPLINEINIT'); end
	% use natural spline if first derivative at x1 not given
	if(dpdx1 ~= 'none'), r1 = 6.0 * ((f(2)-f(1))/(x(2)-x(1)) - dpdx1);
	else, r1 = 0; end
	second(1) = r1;

	% loop through each index, finding current beta and r values (r given by spline algorithm),
	% Eq. 3.38
	for j=2:N
		if(j==N) % last element
			bj = 2.0 * (x(N) - x(N-1));
			% use natural spline if first deriavtive at xN not given
			if(dpdxN ~= 'none'), rj = -6.0 * ((f(N)-f(N-1))/(x(N)-x(N-1)) - dpdxN);
			else, rj = 0; end
		else
			bj = 2.0 * (x(j+1) - x(j-1));
			rj = 6.0 * ((f(j+1) - f(j)) / (x(j+1) - x(j)) - (f(j) - f(j-1))/(x(j) - x(j-1)));
		end
		
		% evaluate the off-diagonal elements. Since the matrix is symmetric, a and c are
		% equivalent. This is Eq. 3.38 and 3.49
		aj = x(j) - x(j-1);
		c = aj;
		beta(j) = bj - aj * c / beta(j-1);
		second(j) = rj - aj * second(j-1) / beta(j-1);
		if(beta(j) == 0), error('Zero diagonal element in SPLINEINIT'); end
	end
	
	% finish with back substitution, using Eq. 3.50, 2.52
	second(N) = second(N) / beta(N);
	for j=1:N-1
		c = x(N-j+1)-x(N-j);
		second(N-j) = (second(N-j) - c*second(N-j+1))/beta(N-j);
	end	
end
