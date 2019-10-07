function [L, U, X] = DLU_fact(A,B)

% This function perfoms Doolittle LU Factorization

% B = matrix of solution vectors

%Allocation of space and setup
U = A;                    %allocate space for U
X = B;                    %allocate space for X

n=size(A,1);              %number of unknowns
m = size(B,2);            %number of b vectors
b_prime = zeros(n,1);
L = ones(1,n);            %set up lower triang matrix
L = diag(L);

%perform simple elimination to obtain L and U
for ir1 = 1:n-1
    for ir2 = ir1+1:n
        elim_factor = U(ir2,ir1); %factor used to eliminate the next row
        for ir3 = ir1:n
            U(ir2,ir3) = U(ir2,ir3)-(elim_factor/U(ir1,ir1))*U(ir1,ir3);
            L(ir2,ir1) =elim_factor/U(ir1,ir1);
        end %for
    end %for
end %for

%Use L and U to calculate X
for i = 1:m
    b_prime = forwardsub_lt(L,B(:,i));
    X(:,i) = backsub(cat(2,U,b_prime));
end %for

end %function






