%%
% EP 501 
% HW 3
% Julio Guardado

clear, clc, close all

%% Problem 1.a + 1.b
%load test problem
load test_lsq.mat

%get first and second order fit using lsqfit.m
yfitlsq1 = lsqfit(x,ynoisy,1);
yfitlsq2 = lsqfit(x,ynoisy,2);

%get fit using polyfit and polyval
ypolyfit1 = polyval(polyfit(x,ynoisy,1),x);
ypolyfit2 = polyval(polyfit(x,ynoisy,2),x);

%plot results of first order fit from lsqfit.m vs polyfit and polyval
figure(1)
%lsqfit.m
subplot(1,2,1)
plot(x,ynoisy,'b.',x,yfitlsq1,'r--')
xlabel('x'); ylabel('y(x)');
title('First order fit using lsqfit.m')
legend('Data','Fit','location','northwest')

%polyfit and polyval
subplot(1,2,2)
plot(x,ynoisy,'b.',x,ypolyfit1,'r--')
xlabel('x'); ylabel('y(x)');
title('First order fit using polyfit')
legend('Data','Fit','location','northwest')

%plot results of second order fit from lsqfit.m vs polyfit and polyval
figure(2)
%lsqfit.m
subplot(1,2,1)
plot(x,ynoisy,'b.',x,yfitlsq2,'r--')
xlabel('x'); ylabel('y(x)');
title('Second order fit using lsqfit.m')
legend('Data','Fit','location','northwest')

%polyfit and polyval
subplot(1,2,2)
plot(x,ynoisy,'b.',x,ypolyfit2,'r--')
xlabel('x'); ylabel('y(x)');
title('Second order fit using polyfit')
legend('Data','Fit','location','northwest')


%% Problem 1.c + 1.d
%compute fits
yfitlsq3 = lsqfit(x,ynoisy,3);

%compute chi squared values
chi1 = gof(ynoisy,yfitlsq1,sigmay);
chi2 = gof(ynoisy,yfitlsq2,sigmay);
chi3 = gof(ynoisy,yfitlsq3,sigmay);

%determine best order of fit from chi squared values
disp('Chi Squared values for linear, quadratic, and cubic fits:')
fprintf('\t Linear:     %f\n',chi1)
fprintf('\t Quardratic: %f\n',chi2)
fprintf('\t Cubic:      %f\n\n',chi3)
disp('Cubic fit is the best fit')

figure(3)
plot(x,ynoisy,'b.',x,yfitlsq3,'r--')
xlabel('x'); ylabel('y(x)');
title('Fit of best order according to chi squared values')
legend('Data','Fit','location','northwest')


%% Problem 2
%load test problem
load test_interp.mat

%perform interpolation using bilin_interp.m
f_biinterp = bilin_interp(xg,yg,f2D,xgi,ygi);

%interpolate with built in function
[Xg,Yg] = meshgrid(xg,yg);
[Xgi,Ygi] = meshgrid(xgi,ygi);
f_interp2 = interp2(Xg,Yg,f2D,Xgi,Ygi);

%plot original function
figure(4)
imagesc(xg,yg,f2D)
colorbar;
title('Original Function')

%plot bilin_interp.m interpolated function
figure(5)
imagesc(xgi,ygi,f_biinterp)
colorbar;
title('bilin\_interp.m')

figure(6)
imagesc(xgi,ygi,f_interp2)
colorbar;
title('interp2.m')


























