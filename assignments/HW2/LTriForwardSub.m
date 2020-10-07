function  x = LTriForwardSub(A,b)
% This function performs a forward elimination on a lower triangular
% matrix and its RHS
%
% A= LHS matrix
% b= RHS matrix
%The code used to write this function was written by Dr. Zettergen. Find it
%at https://github.com/Zettergren-Courses/EP501_matlab/blob/master/linear_algebra/simple_elim_example.m
nref=length(b);             %size of system 
Awork=cat(2,A,b);           %Catenating A and B matrices
x=zeros(nref,1);
x(1)= Awork(1,nref+1)/Awork(1,1); %Finalised Solution for the first variable 
for ir=2:nref
    fact=Awork(ir,ir);          %diagonal element
    x(ir)=Awork(ir,nref+1);     %Assuming we're only dealing with 1 RHS
    for ic=1:ir-1       
        x(ir)=x(ir)-(Awork(ir,ic)*x(ic));
    end %for
    x(ir) = x(ir)/fact;         %Solution of respective variable
end %for
end %function