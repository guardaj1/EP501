function laplacian = delsqr(f,dx)
% Calculates the 3D laplacian of the scalar field given by f assuming equal
% grid spacing dx and a cube domain

%allocate space
gradx=zeros(size(f));
grady=zeros(size(f));
gradz=zeros(size(f));

divx=zeros(size(f));
divy=zeros(size(f));
divz=zeros(size(f));

%get dimensions
lx = size(f,1);
ly = size(f,2);
lz = size(f,3);

%calculate gradient
%x component of gradient
for ix=2:lx-1
    gradx(:,ix,:)=(f(:,ix+1,:)-f(:,ix-1,:))/2/dx;    %\partial/\partial x
end %for
gradx(:,1,:)=(f(:,2,:)-f(:,1,:))/dx;
gradx(:,lx,:)=(f(:,lx,:)-f(:,lx-1,:))/dx;

%y component of gradient
for iy=2:ly-1
    grady(iy,:,:)=(f(iy+1,:,:)-f(iy-1,:,:))/2/dx;    %\partial/\partial y
end %for
grady(1,:,:)=(f(2,:,:)-f(1,:,:))/dx;
grady(ly,:,:)=(f(ly,:,:)-f(ly-1,:,:))/dx;

%z component of gradient
for iz=2:lz-1
    gradz(:,:,iz)=(f(:,:,iz+1)-f(:,:,iz-1))/2/dx;    %\partial/\partial y
end %for
gradz(:,:,1)=(f(:,:,2)-f(:,:,1))/dx;
gradz(:,:,lz)=(f(:,:,lz)-f(:,:,lz-1))/dx;

%calculate laplacian
%x component of div
for ix=2:lx-1
    divx(:,ix,:)=(gradx(:,ix+1,:)-gradx(:,ix-1,:))/2/dx;    %\partial/\partial x
end %for
divx(:,1,:)=(gradx(:,2,:)-gradx(:,1,:))/dx;
divx(:,lx,:)=(gradx(:,lx,:)-gradx(:,lx-1,:))/dx;

%y component of div
for iy=2:ly-1
    divy(iy,:,:)=(grady(iy+1,:,:)-grady(iy-1,:,:))/2/dx;    %\partial/\partial y
end %for
divy(1,:,:)=(grady(2,:,:)-grady(1,:,:))/dx;
divy(ly,:,:)=(grady(ly,:,:)-grady(ly-1,:,:))/dx;

%z component of div
for iy=2:lz-1
    divz(:,:,iz)=(gradz(:,:,iz+1)-gradz(:,:,iz-1))/2/dx;    %\partial/\partial y
end %for
divz(:,:,1)=(gradz(:,:,2)-gradz(:,:,1))/dx;
divz(:,:,lz)=(gradz(:,:,lz)-gradz(:,:,lz-1))/dx;

laplacian = divx+divy+divz;


end




