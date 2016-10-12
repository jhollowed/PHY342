% Joe Hollowed
% PHY 342
%
% Function to visualize the bisection root finding algorithm
% (expected to be called from bisect.m
%
% Last edited 9/9/16

function handles = bisect_plot(func, dx, err, i, handles, xl, xr, xc) 
	
	% param func: function being handled in root finding
	% param dx: desired tolerance
	% param err: approximate relative error in root estiamtion
	% param i: current iteration
	% param xl: left x-bound of root region
	% param xr: right x-bound of root region
	% param xc: bisection point of root region
	% return: vector of figure and plot handles
	
	if(i == 0)
		handles = [1,0,0,0];
	end
	
	fig = handles(1);
	lbound = handles(2);
	rbound = handles(3);
	cbound = handles(4);
	h = get(0, 'screensize')(4);
	w = get(0, 'screensize')(3);

	if(i==0)
		figure(fig, 'Position', [w/4, h/4, w/2, h/2]);
		hold
		x = linspace(-3,3,100);
		plot(x, eval(func), 'Color', 'black', 'linewidth', 2);
		xbounds = xlim;
		ybounds = ylim;
		xlabel('x', 'fontsize', 19)
		ylabel('f(x)', 'fontsize', 19)
		plot(x, 0*x, 'Color', 'black');
		lbound = plot([0,0], ylim, 'Color', 'red', 'linewidth', 1.5);
		rbound = plot([0,0], ylim, 'Color', 'red', 'linewidth', 1.5);
		cbound = plot([0,0], ylim, '--', 'Color', 'blue', 'linewidth', 1.5);
		handles = [fig, lbound, rbound, cbound];	

	else
		figure(fig, 'Position', [w/4, h/4, w/2, h/2])
		hold
	end
		
	funcstr = strrep(func, '*', '');
	funcstr = strrep(funcstr, '.','');
	stattitle = sprintf('f(x) = %s\ntolerance = %.5f', funcstr, dx);
	uptitle = sprintf('relative err = %.5f\ni = %d', err, i);	
	title(sprintf('%s\n%s', stattitle, uptitle), 'fontsize', 16);
	set(lbound, 'XData', [xl, xl])
	set(rbound, 'XData', [xr, xr])
	set(cbound, 'YData', [0,0])
	if(i==0)
		pause
	end
	pause(0.5)
	set(cbound, 'XData', [xc, xc])
	set(cbound, 'YData', ylim)	
	pause(0.5)		
end
