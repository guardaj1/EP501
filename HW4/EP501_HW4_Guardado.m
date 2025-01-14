
clc;clear;close all;
%% Problem 1a + 1b
%define terms
I = 10;                        %A
mu = 4*pi*10^(-7);      %H/m
a = .005;                      %m

%define grid
lx = 50;
ly = 50;
x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
[X,Y] = meshgrid(x,y);

% Allocate space and calculate B
Bx = 0.*X;
By = 0.*Y;
Bmag = 0.*X;
for ix = 1:lx
    for iy = 1:ly
        %components
        if sqrt(X(ix,iy)^2 + Y(ix,iy)^2) < a
            Bx(ix,iy) = (mu*I/(2*pi*a^2))* sqrt(X(ix,iy)^2 + Y(ix,iy)^2)*(-Y(ix,iy)/sqrt(X(ix,iy)^2 + Y(ix,iy)^2));
            By(ix,iy) = (mu*I/(2*pi*a^2))* sqrt(X(ix,iy)^2 + Y(ix,iy)^2)*(X(ix,iy)/sqrt(X(ix,iy)^2 + Y(ix,iy)^2));
        elseif sqrt(X(ix,iy)^2 + Y(ix,iy)^2) >= a
            Bx(ix,iy) = (mu*I/(2*pi*sqrt(X(ix,iy)^2 + Y(ix,iy)^2)))*(-Y(ix,iy)/sqrt(X(ix,iy)^2 + Y(ix,iy)^2));
            By(ix,iy) = (mu*I/(2*pi*sqrt(X(ix,iy)^2 + Y(ix,iy)^2)))*(X(ix,iy)/sqrt(X(ix,iy)^2 + Y(ix,iy)^2));
        end %if
        
        %magnitude
        Bmag(ix,iy) = sqrt(Bx(ix,iy)^2+By(ix,iy)^2);
        
    end %for y
end %for x

% plot
figure(1)
%bx
subplot(1,3,1)
imagesc(x,y,Bx)
title('Bx')
colorbar;
axis xy

%by
subplot(1,3,2)
imagesc(x,y,By)
title('By')
colorbar;
axis xy

%b
subplot(1,3,3)
imagesc(x,y,Bmag)
title('|B|')
colorbar;
axis xy

%B
figure(2)
quiver(x,y,Bx,By,'AutoScale','off')
title('B')

%% Problem 1c+1d
%calculate curl of B numerically
dx = x(2)-x(1);
dy = y(2)-y(1);
curlb = 0.*X;

%x part of curl
curlx = zeros(size(By));
for ix = 2:lx-1
    curlx(:,ix) = (By(:,ix+1)-By(:,ix-1))/2/dx;
end %for
curlx(:,1) = (By(:,2)-By(:,1))/dx;
curlx(:,lx)=(By(:,lx)-By(:,lx-1))/dx;

%y part of curl
curly = zeros(size(Bx));
for iy = 2:ly-1
    curly(iy,:) = (Bx(iy+1,:)-Bx(iy-1,:))/2/dy;
end %for
curly(1,:) = (Bx(2,:)-Bx(1,:))/dy;
curly(ly,:)=(Bx(ly,:)-Bx(ly-1,:))/dy;

numeric_curl = curlx - curly;

%calculate curl of b analytically
analytic_curl = zeros(size(Bx));
for ix = 1:lx
    for iy = 1:ly
        if sqrt(x(ix)^2 + y(iy)^2) < a
            analytic_curl(ix,iy) = mu*I/(2*pi*a^2);
        elseif sqrt(x(ix)^2 + y(iy)^2) >= a
            analytic_curl(ix,iy) = 0;
        end
    end
end

figure(3)
subplot(1,2,1)
imagesc(x,y,numeric_curl)
axis xy
title('Numerical Curl')

subplot(1,2,2)
imagesc(x,y,analytic_curl)
axis xy
title('Analytic Curl')

%% Problem 1e+1f
%define constants
Q = 1;
a = 1;
eps = 8.854e-12;

%define grid
lx = 100;
ly = 100;
lz = 100;

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
z = linspace(-3*a,3*a,lz);

[X,Y,Z] = meshgrid(x,y,z);

%allocate space and calculate potential
phi = 0.*X;
for ix = 1:lx
    for iy = 1:ly
        for iz = 1:lz
            if sqrt(x(ix)^2+y(iy)^2+z(iz)^2) < a
                phi(iy,ix,iz) = Q/(4*pi*eps*a) - (Q/(8*pi*eps*a^3))*(x(ix)^2+y(iy)^2+z(iz)^2-a^2);
            elseif sqrt(x(ix)^2+y(iy)^2+z(iz)^2) >= a
                phi(iy,ix,iz) = Q/(4*pi*eps*sqrt(x(ix)^2+y(iy)^2+z(iz)^2));
            end %if
        end %for z
    end %for y
end %for x

%plot
figure(4);
imagesc(x,y,phi(:,:,50))
title('Scalar Potential')
axis xy
colorbar;

%calculate laplacian
figure(5);
lapl_phi = delsqr(phi,x(2)-x(1));
imagesc(lapl_phi(:,:,50))
title('Laplacian of Potential')
colorbar;
axis xy;


%% Problem 2a
%define grid
lx = 100;
ly = 100;
lz = 100;

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
z = linspace(-3*a,3*a,lz);

%define function to be integrated
func = -.5*eps*lapl_phi.*phi;

%integrate in x
W_x = zeros(ly,lz);
for iy = 1:ly-1
    for iz = 1:lz-1
        W_x(iy,iz) = num_int(func(iy,:,iz),x);
    end %for z
end %for y

%integrate in y
W_y = zeros(lz,1);
for iz = 1:lz-2
    W_y(iz) = num_int(W_x(:,iz),y);
end %for z

%integrate in z
W_E = num_int(W_y,z(1:end-3));

%display result
disp('Electrostatic Energy in the region R in Joules is: ')
disp(W_E)

%% Problem 2b
a = .005;         %m
r0 = 2*a;

%calculate parametric path
lp = 100;
phi = linspace(0,2*pi,lp);
dp = phi(2)-phi(1);

x_phi = r0*cos(phi);
y_phi = r0*sin(phi);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);

%plot path with 
figure(6)
imagesc(x,y,Bmag)
hold on
plot(x_phi,y_phi)
hold off
title('|B| with Parametric Path')

%% Problem 2c
%calculate B at r
B_phix = zeros(size(x_phi));
B_phiy = zeros(size(y_phi));
for ip = 1:lp
    B_phix(ip) = (mu*I/(2*pi*sqrt(x_phi(ip)^2 + y_phi(ip)^2)))*(-y_phi(ip)/sqrt(x_phi(ip)^2 + y_phi(ip)^2));
    B_phiy(ip) = (mu*I/(2*pi*sqrt(x_phi(ip)^2 + y_phi(ip)^2)))*(x_phi(ip)/sqrt(x_phi(ip)^2 + y_phi(ip)^2));
end %for

B_phi = sqrt(B_phix.^2 + B_phiy.^2);

%plot components
figure(7)
plotyy(x_phi,B_phix,y_phi,B_phiy)
title('Components of B along the Parametric Path')
%% Problem 2d
%allocate space for derivatives
drdp_x = zeros(size(x_phi));
drdp_y = drdp_x;

%forward difference at beginning
drdp_x(1) = (x_phi(2)-x_phi(1))/dp;
drdp_y(1) = (y_phi(2)-y_phi(1))/dp;

%take centered diff for middle
for i = 2:lp-1
    drdp_x(i) = (x_phi(i+1) - x_phi(i-1))/(2*dp);
    drdp_y(i) = (y_phi(i+1) - y_phi(i-1))/(2*dp);
end %for

%backward difference at end
drdp_x(lp) = (x_phi(lp)-x_phi(lp-1))/dp;
drdp_y(lp) = (y_phi(lp)-y_phi(lp-1))/dp;

%calculate analytical derivative for comparison
drdp_xreal = -r0*sin(phi);
drdp_yreal = r0*cos(phi);

%plot derivatives
figure(8);
subplot(1,2,1)
plot(drdp_x,drdp_y)
title('Numerically Computed Derivative')
subplot(1,2,2)
plot(drdp_xreal,drdp_yreal)
title('Analytical Derivative')


%% Problem 2e
%calculate integrand by taking dot product
bxdl = B_phix.*drdp_x;
bydl = B_phiy.*drdp_y;
Bdl = bxdl+bydl;

%perform numerical integral
current = num_int(Bdl,phi)/mu;

%display result

disp('Current in the loop in Amps is: ')
disp(current)











