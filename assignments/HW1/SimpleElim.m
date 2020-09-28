%Simple Elimination
function [x] = SimpleElim(A,b)
%SimpleElim 
% This function performs both a simple back elimination and a simple 
%forward eliminatio on a given square matrix and its RHS
%
% A= LHS matrix
% b= RHS matrix
%The code used to write this function was modified from a code written by Dr. Zettergen. Find it
%at https://github.com/Zettergren-Courses/EP501_matlab/blob/master/linear_algebra/simple_elim_example.m
nref=length(b);
x = SForwardElim(A,b,0);
x = SBackElim(x,b,1);
for ir=1:nref
    x(ir,:)=x(ir,:)/(x(ir,ir));
end
end %function