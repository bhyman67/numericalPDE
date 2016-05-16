% Brent Hyman
% Tridiagonal Matrix Algorithm
function [output] = TMA( nx,a,d,c,u )
	for i = 2:nx
		d(i) = d(i)-( a(i-1)/d(i-1) )*c(i-1);
		u(i) = u(i)-( a(i-1)/d(i-1) )*u(i-1);
	end
	v(nx) = u(nx)/d(nx);
	for i = nx-1:-1:1
		v(i) = ( u(i)-c(i)*v(i+1) )/d(i);
	end
	output = v;
end