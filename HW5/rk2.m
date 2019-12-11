function [vx,vy] = rk2(vx0,vy0,t,omega)
lt = length(t);
dt = t(2)-t(1);
vx=zeros(1,lt);
vy=zeros(1,lt);
vy(1)=vy0;
vx(1)=vx0;
for n=2:lt
    %step x and y components together, this is the half update
    vxhalf=vx(n-1)+dt/2*(omega*vy(n-1));
    vyhalf=vy(n-1)-dt/2*(omega*vx(n-1));
    
    %now the full update
    vx(n)=vx(n-1)+dt*(omega*vyhalf);
    vy(n)=vy(n-1)-dt*(omega*vxhalf);    
end %for