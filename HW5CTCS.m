% Brent Hyman
% Crank-Nicolson (CTCS) approximation code for PDE given in homework 5.


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
% matrix that obtains b vector per time step
matrix = zeros(nx,nx);
matrix(1,1) = 2-2*s;
matrix(1,2) = s;
j = 1;
for i = 2:nx-1
	matrix(i,j) = s;
	matrix(i,j+1) = (2-2*s);
	matrix(i,j+2) = s;
	j = j+1;
end
matrix(nx,nx-1) = s;
matrix(nx,nx) = 2-2*s;
					% Algorithm
%===========================================================
for i = 1:nx
	u(i) = IC(i*h);
end
CTCSapprox = u;
for i = 1:nx-1
	a(i) = -s;
	c(i) = -s;
end
for j = 1:nt
	% set up the b vector
	%===	
	bVector = (matrix*CTCSapprox') + 4*k;
	%===
	for i = 1:nx
		d(i) = 2+2*s;
	end
	%solve the system 
	% -> calls the TMA funct to implement the Tridiagonal Matrix Algorithm
	newU = TMA( nx,a,d,c,bVector );
	CTCSapprox = newU;
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
	CTCSapproxActual(a) = CTCSapprox(a-1);
end
CTCSapproxActual(1) = 0;
CTCSapproxActual(nx+2) = 0;
plot( xVect,CTCSapproxActual )
title('CTCS approximation plotted with the true solution (@t = 0.25)');
xlabel('x');
ylabel('u(x,0.25)');
legend('exactSoln','approxSoln')
hold off
					%Chart
%===========================================================
disp('CTCS chart for each spacial step at time t = 0.25')
fprintf('%6s %8s %8s \n','trueSol','approx','absError')
for l = 1:nx+2
    fprintf('%6.6f %6.6f %6.6f\n',exactSolVect(l),CTCSapproxActual(l),abs(exactSolVect(l)-CTCSapproxActual(l)))
end