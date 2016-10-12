% Joe Hollowed
% PHY 342
%
% Function to solve tridiagonal set of equations. Adapted from trisolve.m (Devries & Hasbun 2011)
%
% Last edited 9/24/16

function [x] = trisolve(a, b, c, r)

	% param a: lower diagonal
	% param b: middle diagonal
	% param c: upper diagonal
	% param r: LHS column vector
	% param N: dimension of square matrix
	% return: solution column vector x

	if(b(1) == 0) error('Zero diagonal element in TRISOLVE'); end
		
	N = length(r);
	beta = zeros(1,N);
	x = zeros(1,N);
	beta(1) = b(1);
	x(1) = r(1);

	for j=2:N
		% Eq. 3.49
		beta(j) = b(j) - a(j-1) * c(j-1) / beta(j-1);
		x(j) = r(j) - a(j-1) * x(j-1) / beta(j-1);
		if(b(j) == 0)
			error('Zero diagonal element in TRISOLVE')
		end
	end
	
	% Back substitution (Remove the dependence on x(j-1) from every x(j))
	% Eq. 3.50
	x(N) = x(N) / beta(N);
	for j=1:N-1
		% Eq. 3.51
		x(N-j) = (x(N-j) - c(N-j)*x(N-j+1)) / beta(N-j);
	end
end
