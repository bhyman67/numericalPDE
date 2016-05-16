% Brent Hyman
% FTCS approximation code for PDE given in homework 5.


				%Algorithm input/setup
%===========================================================
clear all;
%step sizes:
% -> space variable step size
h = 0.1;
% -> time variable step size
k = 0.001; 
%Constants:
s = k/(h^2);
tMax = 0.25;
nx = (1-h)/h;
nt = tMax/k;
w = zeros(nx+2,1);
u = zeros(nx+2,1);
					% Algorithm
%===========================================================
%set up initial condition
for i = 2:nx+1
	w(i) = IC(i*h);
end 
%FTCS approximation
for j = 1:nt
	u(1) = 0;
	u(nx+2) = 0;
	for i = 2:nx+1
		u(i) = s*w(i+1) + (1-2*s)*w(i) + s*w(i-1) + 2*k;
	end
	w = u;
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
% -> The FTCS approximation
plot( xVect,u )
title('FTCS approximation plotted with the true solution (@t = 0.25)');
xlabel('x');
ylabel('u(x,0.25)');
legend('exactSoln','approxSoln')
hold off
					%Chart
%===========================================================
disp('FTCS chart for each spacial step at time t = 0.25')
fprintf('%6s %8s %8s \n','trueSol','approx','absError')
for l = 1:nx+2
    fprintf('%6.6f %6.6f %6.6f\n',exactSolVect(l),u(l),abs(exactSolVect(l)-u(l)))
end