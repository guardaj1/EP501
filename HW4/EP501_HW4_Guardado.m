
clc;clear;close all;
%% Problem 1a + 1b
%define terms
I = 10;                        %A
mu = 4*pi*10^(-7);      %H/m
a = .005;                      %m

%define grid
lx = 100;
ly = 100;
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
subplot(1,4,1)
imagesc(x,y,Bx)
title('Bx')
colorbar;
axis xy

%by
subplot(1,4,2)
imagesc(x,y,By)
title('By')
colorbar;
axis xy

%b
subplot(1,4,3)
imagesc(x,y,Bmag)
title('|B|')
colorbar;
axis xy

%bx
subplot(1,4,4)
quiver(x,y,Bx,By,0)
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

figure(2)
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
                phi(iy,ix,iz) = Q/(4*pi*eps*a) + (Q/(8*pi*eps*a^3))*(x(ix)^2+y(iy)^2+z(iz)^2-a^2);
            elseif sqrt(x(ix)^2+y(iy)^2+z(iz)^2) >= a
                phi(iy,ix,iz) = Q/(4*pi*eps*sqrt(x(ix)^2+y(iy)^2+z(iz)^2));
            end %if
        end %for z
    end %for y
end %for x

%plot
figure(3);
imagesc(x,y,phi(:,:,25))
title('Scalar Potential')
axis xy
colorbar;

%calculate laplacian
figure(4);
lapl_phi = delsqr(phi,x(2)-x(1));
imagesc(lapl_phi(:,:,25))
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

%% Problem 2b
a = .005;         %m
r0 = 2*a;

%calculate parametric path
phi = linspace(0,2*pi,100);
x_phi = r0*cos(phi);
y_phi = r0*sin(phi);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);

%plot path with 
figure(1)
subplot(1,4,3)
hold on
plot(x_phi,y_phi)
hold off
























