% Joe Hollowed
% PHY342
%
% This function preforms one step of the fourth order Runge-Kutta
% method to solve the general ODE dy/dx. This is just a Taylor
% series for y(x), truncated to four terms. 
%
% Last edited 10/10/16

function [yf] = rk4(xi, yi, h, ODE)

	% param xi: initial x-value (x at t = 0)
	% param yi: value of y(xi)  (y at t = 0)
	% param h: stepsize
	% param ODE: function handle to ODE dy/dx
	%	     (expected inputs of (x, y) )

	xmid = xi + (h/2);
	xf = xi + h;

	%approximations of derivative dy/dx
	f0 = ODE(xi, yi);
	f1 = ODE(xmid, yi + (h/2)*f0);
	f2 = ODE(xmid, yi + (h/2)*f1);
	f3 = ODE(xf, yi + h*f2);

	yf = yi + (h/6)*(f0 + 2*f1 + 2*f2 + f3);

end
