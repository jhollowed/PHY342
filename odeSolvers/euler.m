% Joe Hollowed
% PHY342
%
% This function preforms one step of Euler's method 
% to solve the general ODE dy/dx. This is just a Taylor
% series for y(x), truncated to two terms. 
%
% Last edited 10/10/16

function [yf] = euler(xi, yi, h, ODE)

	% param xi: initial x-value (x at t = 0)
	% param yi: value of y(xi)  (y at t = 0)
	% param h: stepsize
	% param ODE: function handle to the ODE dy/dx
	%	     (expected inputs of (x, y) )

	yf = yi + (h*ODE(xi, yi));
end
