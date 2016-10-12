% Joe Hollowed
% PHY342
%
% This function preforms one step of the modified Euler's method 
% to solve the general ODE dy/dx. This is just a Taylor series 
% for y(x), truncated to two terms, just like Euler's method, except
% we now use the value of the derivative MIDWAY through each step 
% pf size h to advance the solution 
%
% Last edited 10/10/16

function [yf] = modEuler(xi, yi, h, ODE)

	% param xi: initial x-value (x at t = 0)
	% param yi: value of y(xi)  (y at t = 0)
	% param h: stepsize
	% param ODE: function handle to the ODE dy/dx
	%	     (expected inputs of (x, y) )
	
	xmid = xi + (h/2);
	ymid = yi + ((h/2) * ODE(xi, yi));
	
	yf = yi + (h*ODE(xmid, ymid));
end
