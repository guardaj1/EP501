% EP 501 Homework 1 Main
% Julio Guardado

clear;
clc;
close all;

%% Problem 1.a + 1.b
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 1 Part a+b:')

%load test problem
load testproblem.mat

%used simple elimination function
Amod = forward_elim(A,b);

%test function with built in matlab functions and provided backsub.m
test_ans = A\b;
answer = backsub(Amod);

%display solution
fprintf('\tMATLAB:\t  forward_elim.m:\n')
disp(cat(2,test_ans,answer))

%% Problem 1.c
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 1 Part c+d:')

%load test problem
load lowertriang_testproblem.mat

%solve using matlab built in func
test_ans = L\bL;

%solve using forwardsub func
answer = forwardsub_lt(L,bL);

%display solution
fprintf('\tMATLAB:\t  forward_sub.lt.m:\n')
disp(cat(2,test_ans,answer))

%% Problem 2
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 2:')

%load test problem
load testproblem.mat

%solve using matlab built in func
test_ans = inv(A);

%solve using forwardsub func
[~, invA] = GJ_elim(A,eye(8));

%display solution
disp('Solution from built in MATLAB function:')
disp(test_ans)
disp('Solution from GJ_elim.m:')
disp(invA)

%% Problem 3.a-c
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 3.a-c:')

%load test problem
load testproblem.mat

%perform Doolittle LU factorization
[L, U,X] = DLU_fact(A,cat(2,b,b2,b3));

%display solution of first b matrix
disp('Solution of test problem from L and U:')
disp(X(:,1))

%display solution of multiple right hand sides
disp('Solution of multiple right hand sides from L and U:')
disp(X)

%% Problem 3.d
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 3.d:')

%solve using matlab built in func
test_ans = inv(A);

%solve using LU Factorization
[~,~,X] = DLU_fact(A,eye(size(A,1)));

disp('Solution from built in MATLAB function:')
disp(test_ans)
disp('Solution from Doolittle LU factorization:')
disp(X)

%% Problem 4
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 4:')

%load test problem
load iterative_testproblem.mat

%set up
nit = size(Ait,1);
x0=randn(nit,1);
tol=1e-9;
omega = 1;

%solve using matlab built in func
test_ans = Ait\bit;

%solve using Successive over relaxation
[xit,nit] = SOR(x0,Ait,bit,tol,omega);

%display solution
fprintf('\tMATLAB:\t  SOR.m:\n')
disp(cat(2,test_ans,xit))

%% Problem 5
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Problem 5:')

%load test problem
load testproblem.mat

%solve using matlab built in func
test_ans = det(A);

%solve using determintant.m
[~,determ] = Gauss_elim_det(A,b);

%display solution
fprintf('Determinant calculated by MATLAB:\t\t\t %f\n',test_ans)
fprintf('Determinant calculated by Gauss_elim_det.m,: %f\n',determ)










