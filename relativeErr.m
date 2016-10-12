% Joe Hollowed
% PHY 342
%
% Function to calculate relative error of two values
%
% Last Edited 9/9/16

function err = relativeErr(true, approx)
	
	% param true: "true" or theoretical value
	% param approx: approximate or experimental value
	% return: the relative error

	err = abs((true - approx) / true);
	return
