%% 
%Aldous George
%EP 501
%Project 1
%This code contains excerpts from the simple_elim_example provided by Dr.
%Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/linear_algebra/simple_elim_example.m
clc
clearvars
close all

%%Problem 1

%a 
load testproblem.mat    
ForwardElim= SForwardElim(A,b);

%b
x=backsub(ForwardElim);
disp('Elimination/back sub solution:  ');
disp(x);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(A\b);

%c
load lowertriang_testproblem.mat
TriForElim=LTriFowardrSub(L,bL);
disp(TriForElim);

%d verification
disp('Matlab,GNU/Octave built-in solution:  ');
disp(L\bL);

%% Problem 2
load testproblem.mat    
BackElim= SBackElim(A,b);