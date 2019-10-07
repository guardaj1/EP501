function Amod = forward_elim_mRHS(A,B)

% This function perfoms simple forward elimination for systems with
% multiple right hand sides

% B = matrix of solution vectors

%Allocation of space and setup
Amod=cat(2,A,B);          %make a copy of A and modify with RHS of system
n=size(A,1);              %number of unknowns
m = size(B,2);

%perform simle elimination
for b = 1:m
    for ir1 = 1:n-1
        for ir2 = ir1+1:n
            elim_factor = Amod(ir2,ir1); %factor used to eliminate the next row
            for ir3 = ir1:n+m
                Amod(ir2,ir3) = Amod(ir2,ir3)-(elim_factor/Amod(ir1,ir1))*Amod(ir1,ir3);
            end %for
        end %for
    end %for
end %for
end %function