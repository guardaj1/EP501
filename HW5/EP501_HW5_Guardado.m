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

%% Problem 1b
%calculate derivative of epsilon
deps_dx = zeros(size(eps));

%forward diff for beginning 
deps_dx(1) = (eps(2)-eps(1))/dx;

%take centered diff for middle
for i = 2:lx-1
    deps_dx(i) = (eps(i+1) - eps(i-1))/(2*dx);
end %for

%backward difference at end
deps_dx(lx) = (eps(lx)-eps(lx-1))/dx;

%make matrix
M = zeros(lx,lx);
b = zeros(lx,1);

%build matrix
M(1,1) = -1/dx;
M(1,2) = 1/dx; 
M(lx,lx) = 1;
for ix = 2:lx-1
    %alpha
    M(ix,ix-1) = eps(ix)/dx^2 + deps_dx(ix)/(2*dx);
    
    %beta
    M(ix,ix) = -2*eps(ix)/dx^2;
    
    %gamma
    M(ix,ix+1) = eps(ix)/dx^2 - deps_dx(ix)/(2*dx);
    
end %for

%apply BCs
b(1) = 1000;
b(lx) = 100;

%solve
soln = M\b;

figure(2)
plot(x,soln)
title('Electric Potential')
xlabel('x(m)'); ylabel('Potential (V)')
%% Problem 2
%set up
m = 1.67e-27;   %kg
q = 1.6e-19;    %C
B = 50000;      %nT

%particle oscillation period
T = 2*pi*m/q/B;
omega = 2*pi/T;

%time grid
lt = 100;
t = linspace(0,T,lt);
dt = t(2)-t(1);

%rk2
vx_rk2=zeros(1,lt);
vy_rk2=zeros(1,lt);
vy_rk2(1)=1e3;
vx_rk2(1)=1e3;
for n=2:lt
    %step x and y components together, this is the half update
    vxhalf=vx_rk2(n-1)+dt/2*(omega*vy_rk2(n-1));
    vyhalf=vy_rk2(n-1)-dt/2*(omega*vx_rk2(n-1));
    
    %now the full update
    vx_rk2(n)=vx_rk2(n-1)+dt*(omega*vyhalf);
    vy_rk2(n)=vy_rk2(n-1)-dt*(omega*vxhalf);    
end %for

%rk4
vx_rk4=zeros(1,lt);
vy_rk4=zeros(1,lt);
vy_rk4(1)=1e3;
vx_rk4(1)=1e3;
for n=2:lt
    %calculate vx and vy at several time steps
    vx1=dt*fRK(vy_rk4(n-1),omega);
    vx2=dt*fRK(vy_rk4(n-1)+vx1/2,omega);
    vx3=dt*fRK(vy_rk4(n-1)+vx2/2,omega);
    vx4=dt*fRK(vy_rk4(n-1)+vx3,omega);
    
    vy1=dt*fRK(vx_rk4(n-1),omega);
    vy2=dt*fRK(vx_rk4(n-1)+vy1/2,omega);
    vy3=dt*fRK(vx_rk4(n-1)+vy2/2,omega);
    vy4=dt*fRK(vx_rk4(n-1)+vy3,omega);
    
    %now the full update
    vx_rk4(n)=vx_rk4(n-1)+1/6*(vx1+2*vx2+2*vx3+vx4);
    vy_rk4(n)=vy_rk4(n-1)+1/6*(vy1+2*vy2+2*vy3+vy4);    
end %for

%plot rk2 vs rk4
figure(3)
ax = plotyy(t,vx_rk2,t,vy_rk2);
title('RK2')
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');

figure(4)
ax = plotyy(t,vx_rk4,t,vy_rk4);
title('RK4')
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');

%% Problem 3
tt = linspace(0,.5,100);









