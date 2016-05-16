% Brent Hyman
% the initial condition in time
function[output] = IC(x)
	output = sin(pi*x) + x*(1-x);
end