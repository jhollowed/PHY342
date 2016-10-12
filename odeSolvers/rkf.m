% Joe Hollowed
% PHY342
%
% This function solves the general ODE dy/dt via the Runge-Kutta
% Fehlberg method, which uses adaptive stepsizes based on the 
% 4th and 5th order Runge-Kutta approximations.  
%
% Last edited 10/11/16

function [t, y] = rkf(t0, tf, y0, ODE, h_min = 0.00001, h_max = 0.1,...
		      epsilon = 1e-5, maxFactor = 4.0, minFactor = 0.1)

	%----------------------------------------------------------------------
	
	% param t0: initial time
	% param tf: final time
	% param y0: initial conditions (can be a vector for systems of ODE's. 
	%	    Vector can be either row or column, doesn't matter)
	% param ODE: function handle for ODE, expecting inputs (t, y)
	%	     (for example, if your ODE is, dp/dt = mg -kv^2 (prob 5.6),
	%	      then the ODE argument, say 'dpdt', should be defined as 
	%	      'dpdt = @(t, p) (100*9.8 - 0.1*(p/100)^2)', where I assumed
	%	      m = 100, g = 9.8, k = 0.1. It does not matter if the 
	%             variable names match those of this function, as long as the 
	%	      dependent and independent vars are in the right order. This 
	%	      is basically just turning the derivs() function from MyRKF.m
	%	      into an input parameter, so you never need to edit any source 
	%	      code per problem.
	% param h_min: minimum step size at which to halt execution
	%	       (algorithm should not reach this step size on
	%		a successful approximation).
	% param h_max: maximum allowed (and starting) value for h
	% param epsilon: desired accuracy (maximum accepted error for each step)	
	% param maxFactor: the most that h_new can increase by, relative to h, in 
	%		   a single step
	% param minFactor: the least that h_new can decrease by, relative to 
	%		   h, in a single step
	% return: the discretized function y(t), as the arrays t and y
	%	  where each row of y is the function values of each 
	%	  dependent variable (for systems of ODES)
	%
	%
	% Note that most of the parameters have a value already assigned to them. 
	% This makes them optional, and you only need to pass them if you want to
	% change the default value. You only need to pass the first 4 args
	
	%------------------------------------------------------------------------

	% initialize som stuff
	if(size(yi)(1) ~= 1), y = y'; end;
	dimensions = length(y0);
	t(1) = t0;
	y(1, 1:dimensions) = y0;
	h = h_max;
	i = 0;
	
	% read coefficients from file
	a_file = textread('coeffs.dat', '%s');
	c_file = textread('y_coeffs.dat', '%s');
	a_coeffs = cellfun(@eval, a_file);
	c = cellfun(@eval, c_file);	

	% loop until we reach 'tf', increasing 't' by some optimal step size 'h' each time
	while(t(end) < tf)
		h_found = false;
		i = i+1;
		yi = y(i, :);
			
		% loop until we find an optimal value for 'h', which satisfies the
		% error tolerance 'epsilon'
		while(h_found == false)
			% find each intermediate function value, by evaluating the ODE using 
			% RKF coefficients	
			a = a_coeffs*h;
			f0 = ODE(t(i), yi);
			f1 = ODE(t(i)+a(1), yi+a(2)*f0);
			f2 = ODE(t(i)+a(3), yi+a(4)*f0 + a(5)*f1);
			f3 = ODE(t(i)+a(6), yi+a(7)*f0 + a(8)*f1 + a(9)*f2);
			f4 = ODE(t(i)+a(10), yi+a(11)*f0 + a(12)*f1 + a(13)*f2 + a(14)*f3);
			f5 = ODE(t(i)+a(15), yi+a(16)*f0 + a(17)*f1 + a(18)*f2 + a(19)*f3 + a(20)*f4);
			
			% find 'y' and 'yhat', the 4th and 5th order RK approximations. 
			% compare them to find error (depending on number of dimensions) 
			% in the solution.
			y_n = yi + h*(c(1)*f0 + c(2)*f2 + c(3)*f3 + c(4)*f4);
			yhat_n = yi + h*(c(5)*f0 + c(6)*f2 + c(7)*f3 + c(8)*f4 + c(9)*f5);
			
			if(dimensions == 1)	
				relErr = abs((yhat_n - y_n) / yhat_n);
			else
				relErr = sqrt( sum(((yhat_n - y_n)./yhat_n).^2)) / dimensions;
			end
			h_new = (0.9*h) * ((epsilon/relErr) ^(0.2));
			
		
			% check quality of new stepsize. If too small, assume there is a 
			% problem, return. If too large, use max stepsize
			% if 'h_new' is larger than the previous stepsize 'h', then we were
			% already below our error tolerance with 'h', use it and advance 
			% solution, and storing next values of 'y' and 't'
			if(h_new > maxFactor*h), h_new = maxFactor*h;
			elseif(h_new < minFactor*h), h_new = minFactor*h; end;
			if(h_new < h_min), error('h is getting too small!');
			elseif(h_new > h_max), h_new = h_max; end;
			if(h_new >= h)
				h_found = true;	
				t(i+1) = t(i) + h;	
				y(i+1, 1:dimensions) = yhat_n;
			end
			% 'h_new' was not larger than 'h', meaning we can do better. 
			% Redefine 'h' and loop again
			h = h_new;
		end
	end
end
