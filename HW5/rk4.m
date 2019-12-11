function yRK4 = rk4(t,y0,f)

lt = length(t);
ly = length(y0);
dt = t(2)-t(1);

yRK4=zeros(ly,lt);
yRK4(:,1)=y0;

for n=2:lt
    dy1=dt*f(t(n-1),yRK4(:,n-1));
    dy2=dt*f(t(n-1)+dt/2,yRK4(:,n-1)+dy1/2);
    dy3=dt*f(t(n-1)+dt/2,yRK4(:,n-1)+dy2/2);
    dy4=dt*f(t(n-1)+dt,yRK4(:,n-1)+dy3);
    
    yRK4(:,n)=yRK4(:,n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for