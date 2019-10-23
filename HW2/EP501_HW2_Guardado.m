%EP 501
%Homework 2
%Julio Guardado

clc; clear; close all;

%% Problem 1, Part a: first root of 0 order bessel function

%define function and setup
rho = linspace(0,20,200);
eta = .01;
maxit = 100;
tol = 1e-9;
f = @besselj;
f_rho = besselj(0,rho);

%plot function to find initial guess
figure(1);
plot(rho,f_rho)
xlabel('\rho');
ylabel('J_0(\rho)');
title('Bessel Function of Order 0')

%find first root using approximate newton method
x0 = 2.5;
[xNewtApprox, niter1,~] = newton_approx_bess(x0,eta,maxit,tol);


%% Problem 1,Part C: first six roots of order 0 bessel function
bess0_roots = zeros(1,7);       %allocate space for roots and number of iterations
niter = zeros(1,7);             %index 1 is only used for checking
x0 = 0;
j = 2;
while (j<8)                      %iterate using approx newton function
    x0 = x0 + 3;
    [bess0_roots(j),niter(j),~] = newton_approx_bess(x0,eta,maxit,tol);
    %check for convergence on the same root or out of order root
    if(bess0_roots(j) == bess0_roots(j-1))
        j = j-1;
        x0 = x0 + x0;
    else
        j = j+1;
    end
    
end

%display results
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 1:')
disp('First six roots of 0th order Bessel function:')
disp(bess0_roots)
fprintf('\t(first 5 roots checked against Wolfram Mathworlds Bessel Function Zeros)\n')
fprintf('\t http://mathworld.wolfram.com/BesselFunctionZeros.html\n\n')


%% Problem 2.a
%get function handles
f = @funcHW2a;
fder = @funcderHW2a;

%define terms for root finding
maxit = 100;
tol = 1e-9;
n = 5;                      %order of poly
x0 = linspace(.75,4.75,5);
rootsa = zeros(1,n);

%find roots
for j = 1:n
    rootsa(j) = newton_exact(f,fder,x0(j),maxit,tol);
end


%% Problem 2.b
%get function handles
f = @funcHW2b;
fder = @funcderHW2b;

%define terms for root finding
n = 3;                      %order of poly
x0 = [.5,.5+.5*1i,.5-1i];
rootsb = zeros(1,n);

for j = 1:n
    rootsb(j) = newton_exact(f,fder,x0(j),maxit,tol);
end

%display results
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 1a:')
disp('Roots of x^5-15x^4+95x^3-225x^2+274x-120: ')
disp(rootsa)
disp('Problem 1b:')
disp('Roots of x^3-3x^2+4x-2: ')
disp(rootsb)


%% Problem 4a
%define function handles
fm=@HW2fun2Df3;
gm=@HW2fun2Dg3;
gradfm=@grad_HW2fun2Df3;
gradgm=@grad_HW2fun2Dg3;

%define terms for root finding
nx = 2; %order of x
ny = 2; %order of y
n = nx+ny;  %total order

x0 = [2 4];
y0 = [4 9];
rootx = zeros(1,n);
rooty = zeros(1,n);

%find roots
for i =1:nx
    for j = 1:ny
        [rootx(i),rooty(i),~,~] = newton2D_exact(fm,gradfm,gm,gradgm,x0(i),y0(j));
    end
end


disp(rootx)










