function x=MultiRHSbacksub(A)

% This function performs back substitution on an upper triangular matrix that has
% been modified by concatenating the RHS of the system.
% Note that B is assumed to be upper triangular at this point.


[n,nRHS]=size(A);                   %number of unknowns in the system
x=zeros(n,1);                  %space in which to store our solution vector
x(n)=A(n,n+1)/A(n,n);          %finalized solution for last variable, resulting from upper triangular conversion
for ir1 =1:nRHS
    for ir2=n-1:-1:1
        x(ir2,ir1)=A(ir2,n+1);       %assume we're only dealing with a single right-hand side here.
        fact=A(ir2,ir2);         %diagonal element to be divided through doing subs for the ir2 row
        for ic=ir2+1:n
            x(ir2,ir1)=x(ir2)-A(ir2,ic)*x(ic,ir1);
        end %for
        x(ir2,ir1)=x(ir2,ir1)/fact;      %divide once at the end to minimize number of ops
    end %for
end %for
end %function