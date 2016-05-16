% Brent Hyman
% Actual Solution
function[output] = exactSol(x,t)
	output = exp(-pi^(2)*t)*sin(pi*x) + x*(1-x);
end