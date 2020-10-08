%% Introduction
%Aldous George
%EP 501
%Project 2
%This code contains excerpts from codes provided by Dr. Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/...
%linear_algebra
clc
clearvars
close all
%% Problem 1

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
%% 1 C)

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
%% 1 D)
%d
%Calculating Inverse one column at a time
nref=length(b);
InvA=[];
for ir=1:nref
    B=zeros(nref,1);
    B(ir) = 1;
    %Forward sub for b'
    Bprime=LTriForwardSub(L,B);
    %Back substitution for x
    x=backsub(cat(2,U,Bprime)); %ith column of the Inverse
    InvA=cat(2,InvA,x);
end %for
disp('Inverse of A:  ');
disp(InvA);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(inv(A));
%% Problem 2

%a Kindly refer to SoR.m function

%b
disp('2 B)');
%Initialisation
load 'iterative_testproblem.mat'
nref=size(Ait,1);
nit=10;
x0=zeros(nref,1);
tol=1e-10;
w=1.1;
%Testing Successive over-Relaxation
[xit,iter]=SoR(x0,Ait,bit,tol,false,w);
disp('Solution with Successive over-Relaxation iteration: ')
disp(xit);
disp('Number of Iterations required: ')
disp(iter);
disp('Tolerance: ')
disp(tol);
disp('MATLAB built-in solution: ')
disp(Ait\bit);

%% 2 C)
i=1;
for w=0.38:0.025:1.73
    [xit,iter]=SoR(x0,Ait,bit,tol,false,w);
    Iteration(i)=iter;
    Omega(i)=w;
    i=i+1;
end
[Min, Index] = min(Iteration);
disp('Relaxation parameter (Omega) that minimizes iterations: ')
disp(Omega(Index));
disp('Number of iterations performed with that Relaxation parameter (Omega): ')
disp(Iteration(Index));
%% 2 D)
[xit,iter]=SoR(x0,Ait,bit,tol,false,1);
disp('Number of iterations performed by Gauss-Seidel: ')
disp(iter);
disp('Difference in number of iterations performed by Gauss-Seidel vs Optimal case: ')
disp(abs(iter-Iteration(Index)));
j=1;
for i=1:length(Iteration)
    if iter==Iteration(i)
        ind(j,1)=i;
        ind(j,2)=Omega(i);
        j=j+1;
    end %if
end %for
disp('The lowest value of the Relaxation parameter (Omega) that performs atleast as well as the Gauss-Seidel method: ')
disp(ind(1,2));
disp('The highest value of the Relaxation parameter (Omega) that performs atleast as well as the Gauss-Seidel method: ')
disp(ind(2,2));
