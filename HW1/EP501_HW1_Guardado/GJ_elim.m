function [Amod, B] = GJ_elim(A,B)

%this function performs gauss jordan elimination

%Allocation of space and setup
n=size(A,1);              %number of unknowns
m=size(B,2);              %number of solution vectors

%perform forward elimination using made function
Amod = forward_elim_mRHS(A,B);

%perfom backwards elimination to get get diagonal matrix
for ir1 = 1:n
    for ir2 = ir1+1:n
        elim_factor = Amod(ir1,ir2)/Amod(ir2,ir2);
        for ir3 = ir2:n+m
            Amod(ir1,ir3) = Amod(ir1,ir3) - elim_factor*Amod(ir2,ir3);
        end
    end
end

% %perfom backwards elimination to get get diagonal matrix
% for ir1 = n:-1:1
%     for ir2 = n-1:-1:n
%         elimf = Amod(ir2,ir1);
%         for ir3 = n+m:-1:ir1
%             Amod(ir1,ir3) = Amod(ir1,ir3) - (elimf/Amod(ir1,ir1))*Amod(ir1,ir3);
%         end
%     end
% end

%scale matrix to get identity
scalef = diag(Amod);
for i = 1:n
    for j = 1:n+m
        Amod(i,j) = Amod(i,j)/scalef(i);
    end %for
end %for

A = Amod(:,1:n);
B = Amod(:,n+1:n+m);

end










