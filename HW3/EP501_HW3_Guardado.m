%%
% EP 501 
% HW 3
% Julio Guardado

clear, clc, close all

%% Problem 1
%load test problem
load test_lsq.mat

%create line and quadratic function
yline = 5+6*x;
yquad = 2*x.^2+3*x+6;

%get fit using lsqfit.m
yfitl1 = lsqfit(x,yline);
yfitq1 = lsqfit(x,yquad);

%get fit using polyfit and polyval
yfitl2 = polyfit(x,yline,1);
yfitl2 = polyval(yfitl2,x);
yfitq2 = polyfit(x,yquad,1);
yfitq2 = polyval(yfitq2,x);

%plot line
figure(1);
%Least squares fit
subplot(2,1,1)
plot(x,yline,'b-',x,yfitl1,'r--')
legend('Actual','Fit','location','northwest')
xlabel('x'); ylabel('y(x)')
title('lsqfit.m fit')

%matlab fit
subplot(2,1,2)
plot(x,yline,'b-',x,yfitl2,'r--')
legend('Actual','Fit','location','northwest')
xlabel('x'); ylabel('y(x)')
title('Polyval fit')

%plot quadratic
figure(2);
title('Fit of a quadratic function')
%Least squares fit
subplot(2,1,1)
plot(x,yquad,'b-',x,yfitq1,'r--')
legend('Actual','Fit','location','northwest')
xlabel('x'); ylabel('y(x)')
title('lsqfit.m fit')

%matlab fit
subplot(2,1,2)
plot(x,yquad,'b-',x,yfitq2,'r--')
legend('Actual','Fit','location','northwest')
xlabel('x'); ylabel('y(x)')
title('Polyval fit')




