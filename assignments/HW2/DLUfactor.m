%Doolittle LU Factorization
function [L,U] = DLUfactor(A)
%SForwardElim Simple Forward Elimination
% This function performs a simple forward elimination on a given square
% matrix and its RHS
%
% A= LHS matrix
% b= RHS matrix
% c= 0 if matrix needs to be catenated, and 1 if matrix is already catenated.
%The code used to write this function was written by Dr. Zettergen. Find it
%at https://github.com/Zettergren-Courses/EP501_matlab/blob/master/linear_algebra/simple_elim_example.m

nref=length(A);             %size of system 
L=zeros(nref);
L=eye(nref,nref);
for ir1=2:nref              %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from (i.e. using) for this particular column
    for ir2=ir1:nref        %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
        fact=A(ir2,ir1-1); %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
        A(ir2,:)=A(ir2,:)-fact/A(ir1-1,ir1-1).*A(ir1-1,:);    %subtract off previous row modified by a factor that eliminates the ir-1 column term in this row (so it has only super-diagonal elements)
    L(ir2,ir1-1)=fact/A(ir1-1,ir1-1);
    end %for
end %for
U=A;
end %function