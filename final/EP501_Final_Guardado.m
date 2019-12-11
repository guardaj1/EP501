%% EP 501 Final
% Julio Guardado

clc;clear;close all
%% Problem 3
%% Define constants
% Define space grid
lx = 100;
x = linspace(0,1,lx);
dx = x(2)-x(1);

%constants
lambda = .25;
v=20;            %velocity of wave propagation

%initial condition
finitial=exp(-(x-.5).^8/(2*.1^8));

%% Part b: implicit solution
%set stable time step for implicit solution
dt=1/lambda/2*dx^2;
t=0:dt:.05;
lt=numel(t);

%initialize
fimpli=zeros(lx,lt);
fimpli(:,1)=finitial;

%calculate diffusion
for n=1:lt-1
    fimpli(1,n+1)=0;   %assume temperature goes to some small number on left boundary
    for i=2:lx-1      %interior grid points
        fimpli(i,n+1)=fimpli(i,n)+dt/dx^2*lambda*(fimpli(i+1,n)-2*fimpli(i,n)+fimpli(i-1,n));
    end %for
    fimpli(lx,n+1)=0;  %assume temperature goes to some small number on right boundary
end %for

%calculate advection
for n = 1:lt-1
    fimpli(:,n+1)=Godunov(dt,dx,v,fimpli(:,n));
end

%plot
figure(3)
for n = 1:lt
    imagesc(x,t,fimpli(:,n)' )
    axis xy
    colorbar
    caxis([0 1])
    shading flat
    xlabel('x');
    title(sprintf('Implicit Solution to Advection-Diffusion Eq, t=%5.3f',t(n)));
    
    if n == 1
        pause
    else
    pause(.1)
    end %if
end %for

%% Part c
%set unstable time step
dt2= 1/lambda/2*dx^2;
t2=0:2.5*dt2:.05;
lt2=numel(t2);

funstable=zeros(lx,lt2);
funstable(:,1)=finitial;

for n=1:lt2-1
    funstable(1,n+1)=0;   %assume temperature goes to some small number on left boundary
    for i=2:lx-1      %interior grid points
        funstable(i,n+1)=funstable(i,n)+dt2/dx^2*lambda*(funstable(i+1,n)-2*funstable(i,n)+funstable(i-1,n));
    end %for
    funstable(lx,n+1)=0;  %assume temperature goes to some small number on right boundary
end %for

for n = 1:lt2-1
    funstable(:,n+1)=Godunov(dt2,dx,v,funstable(:,n));
end

%plot instability at t = 0.05
fprintf('Max stable timestep: %f s\n',2.4*dt2) 
figure(4)
plot(x,funstable(:,end))

