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

x = SForwardElim(A,b);
x = SBackElim(x,b,1);
end %function