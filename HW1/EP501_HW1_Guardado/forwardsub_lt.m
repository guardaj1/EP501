function x = forwardsub_lt(A,b)

% This function performs back substitution on an lower triangular matrix that has
% been modified by concatenating the RHS of the system.  
% Note that B is assumed to be lower triangular at this point.

Amod = cat(2,A,b);

n=size(Amod,1);               %number of unknowns in the system
x=zeros(n,1);                 %space in which to store our solution vector
x(1)=Amod(1,n+1)/Amod(1,1);   %finalized solution for last variable, resulting from lower triangular conversion

for ir1=2:n
    x(ir1)=Amod(ir1,n+1);     %assume we're only dealing with a single right-hand side here.
    fact=Amod(ir1,ir1);       %diagonal element to be divided through doing subs for the ir2 row
    for ic=1:ir1-1
        x(ir1) = x(ir1)-Amod(ir1,ic)*x(ic);
    end %for
    x(ir1)=x(ir1)/fact;       %divide once at the end to minimize number of ops
end %for 

end %function