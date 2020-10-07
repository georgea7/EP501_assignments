%% Introduction
%Aldous George
%EP 501
%Project 2
%This code contains excerpts from the simple_elim_example provided by Dr.
%Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/...
%linear_algebra/simple_elim_example.m
clc
clearvars
close all
%% Problem 1
%Example 1.7
% A=[80, 20, -20; -20, 40,-20; -20, -20, 130];
% b=[20; 20; 20];

%a Kindly refer to DLUfactor function

%b
disp('b)');
load 'testproblem.mat'
[L,U] = DLUfactor(A);

%Forward sub for b'
bprime=LTriForwardSub(L,b);

%Back substitution for x
x=backsub(cat(2,U,bprime));

disp('Doolittle LU factorisation, b:  ');
disp(x);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(A\b);

%c
disp('c)');
%forward sub for b2
%Forward sub for b'
b2prime=LTriForwardSub(L,b2);

%Back substitution for x
x=backsub(cat(2,U,b2prime));

disp('Doolittle LU factorisation, b2:  ');
disp(x);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(A\b2);

%forward sub for b3
%Forward sub for b'
b3prime=LTriForwardSub(L,b3);

%Back substitution for x
x=backsub(cat(2,U,b3prime));

disp('Doolittle LU factorisation, b3:  ');
disp(x);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(A\b3);

%d
disp('d)');

nref=length(b);
InvA=[];
for ir=1:nref
    B=zeros(nref,1);
    B(ir) = 1;
    %Forward sub for b'
    Bprime=LTriForwardSub(L,B);
    %Back substitution for x
    x=backsub(cat(2,U,Bprime));
    InvA=cat(2,InvA,x);
end
disp('Inverse of A:  ');
disp(InvA);

%% Problem 2


