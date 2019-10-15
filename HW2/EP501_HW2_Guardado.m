%EP 501
%Homework 2
%Julio Guardado

clc;clear;close all;

%% Problem 1, Part a: first root of 0 order bessel function

%define function and setup
rho = linspace(0,30,200);
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
i = 2;
while (i<8)                      %iterate using approx newton function
    x0 = x0 + 3;
    [bess0_roots(i),niter(i),~] = newton_approx_bess(x0,eta,maxit,tol);
    disp(bess0_roots(i) == bess0_roots(i-1))
    %check for convergence on the same root or out of order root
    if(bess0_roots(i) == bess0_roots(i-1))% || bess0_roots(i) < bess0_roots(i-1))
        i = i-1;
        x0 = x0 + x0;
        disp('ahhhhh')
    else
        i = i+1;
    end
    
end
disp(bess0_roots)

