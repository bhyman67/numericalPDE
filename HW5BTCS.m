% Brent Hyman
% BTCS approximation code for PDE given in homework 5.

				%Algorithm input/setup
%===========================================================
clear all;
%step sizes:
% -> space variable step size
h = 0.1;
% -> time variable step size
k = 0.01; 
%Constants:
s = k/(h^2);
tMax = 0.25;
nx = (1-h)/h;
nt = tMax/k;
					% Algorithm
%===========================================================
for i = 1:nx
	u(i) = IC(i*h);
end
BTCSapprox = u;
for i = 1:nx-1
	a(i) = -s;
	c(i) = -s;
end
for j = 1:nt
	for i = 1:nx
		d(i) = 1+2*s;
	end
	%solve the system 
	% -> calls the TMA funct to implement the Tridiagonal Matrix Algorithm
	newU = TMA( nx,a,d,c,(BTCSapprox+2*k) );
	BTCSapprox = newU;
end
					  %Graph
%===========================================================
% -> The exact solution vector
exactSolVect(1) = exactSol(0,0.25);
exactSolVect(nx+2) = exactSol(1,0.25);
xVect(1) = 0;
xVect(nx+2) = 1;
for x = 1:nx
	exactSolVect(x+1) = exactSol(x*h,0.25);
	xVect(x+1) = x*h;
end
hold on	
plot( xVect,exactSolVect,'R' )
% -> The BTCS approximation
for a = 2:nx+1
	BTCSapproxActual(a) = BTCSapprox(a-1);
end
BTCSapproxActual(1) = 0;
BTCSapproxActual(nx+2) = 0;
plot( xVect,BTCSapproxActual )
title('BTCS approximation plotted with the true solution (@t = 0.25)');
xlabel('x');
ylabel('u(x,0.25)');
legend('exactSoln','approxSoln')
hold off
					%Chart
%===========================================================
disp('BTCS chart for each spacial step at time t = 0.25')
fprintf('%6s %8s %8s \n','trueSol','approx','absError')
for l = 1:nx+2
    fprintf('%6.6f %6.6f %6.6f\n',exactSolVect(l),BTCSapproxActual(l),abs(exactSolVect(l)-BTCSapproxActual(l)))
end