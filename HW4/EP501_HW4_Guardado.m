
clc;clear;close all;
%% Problem 1a + 1b
%define terms
I = 10;                        %A
mu_naught = 4*pi*10^(-7);      %H/m
a = .005;                      %m

%define grid
lx = 50;
ly = 50;
x = linspace(-3*a,3*a,50);
y = linspace(-3*a,3*a,50);
[X,Y] = meshgrid(x,y);

% Allocate space and calculate B
Bx = 0.*X;
By = 0.*Y;
Bmag = 0.*X;
for ix = 1:lx
    for iy = 1:ly
        %components
        if sqrt(X(ix,iy)^2 + Y(ix,iy)^2) < a
            Bx(ix,iy) = (mu_naught*I/(2*pi*a^2))* sqrt(X(ix,iy)^2 + Y(ix,iy)^2)*(-Y(ix,iy)/sqrt(X(ix,iy)^2 + Y(ix,iy)^2));
            By(ix,iy) = (mu_naught*I/(2*pi*a^2))* sqrt(X(ix,iy)^2 + Y(ix,iy)^2)*(X(ix,iy)/sqrt(X(ix,iy)^2 + Y(ix,iy)^2));
        elseif sqrt(X(ix,iy)^2 + Y(ix,iy)^2) >= a
            Bx(ix,iy) = (mu_naught*I/(2*pi*sqrt(X(ix,iy)^2 + Y(ix,iy)^2)))*(-Y(ix,iy)/sqrt(X(ix,iy)^2 + Y(ix,iy)^2));
            By(ix,iy) = (mu_naught*I/(2*pi*sqrt(X(ix,iy)^2 + Y(ix,iy)^2)))*(X(ix,iy)/sqrt(X(ix,iy)^2 + Y(ix,iy)^2));
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

%% Problem 1c
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

curl = curlx - curly;

figure(2)
imagesc(x,y,curl)












