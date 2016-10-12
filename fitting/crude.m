% Joe Hollowed
% PHY342
% 
% Function to preform one step of least squares minimization by varying each 
% parameter one by one
%
% Last edited 10/2/16

function [a, h] = crude(a, h, x, y, Y)

	% param a: vector of coefficients
	% param h: vector of step sizes
	% param y: vector of experimental y-values
	% param Y: function handle for fitting function
	% return: new coefficients and step sizes
	
	for i = 1:length(a)
	%loop through each parameter
		for k = 1:length(a)
			if(k == i)
				%vary only one parameter at a time
				a_plus(i) = a(i) + h(i);
				a_minus(i) = a(i) - h(i);
			else
				a_plus(k) = a(k);
				a_minus(k) = a(k);
			end
		end
	
		%find least-squares at plus and minus h, and h=0
		s0 = sum((y - Y(a, x, y)).^2); 
		sp = sum((y - Y(a_plus, x, y)).^2); 
		sm = sum((y - Y(a_minus, x, y)).^2);

		%Minimum of quadtratic polynomial of the LSE
		a(i) = a(i) - (0.5*h(i)*((sp-sm)/(sp-2.*s0+sm)));
		h(i) = 0.5*h(i);
	end
end
