%%
% EP 501 Homework 5
% Julio Guardado

clc;clear;close all;

%% Problem 1 Setup
%define parameters
a = 0.01;           %m
l = a/5;            %m
x_prime = -9*a/10;  %m
x_dprime = 9*a/10;  %m
eps0 = 8.854e-12;

%define grid
lx = 20;
x = linspace(-a,a,lx);
dx = x(2)-x(1);
dp_dx = zeros(lx,1);
phi = zeros(lx,1);

%define boundary conditions
dp_dx(1) = 1000;        %V/m
phi(lx) = 100;          %V

%% Problem 1a
%calculate dielectric function
eps = eps0*(10*tanh((x-x_prime)/l)-10*tanh((x-x_dprime)/l));

%plot
figure(1)
plot(x,eps)
title('Dielectric Function')
xlabel('x'); ylabel('\epsilon(x)');









