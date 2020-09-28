%% Introduction
%Aldous George
%EP 501
%Project 1
%This code contains excerpts from the simple_elim_example provided by Dr.
%Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/...
%linear_algebra/simple_elim_example.m
clc
clearvars
close all

%% Problem 1

%a 
load testproblem.mat    
ForwardElim= SForwardElim(A,b,0);

%b
x=backsub(ForwardElim);
disp('Elimination/back sub solution:  ');
disp(x);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(A\b);

%c
load lowertriang_testproblem.mat
disp('Lower Triangular Matrix/forward substitution:  ');
TriForElim=LTriFowardrSub(L,bL);
disp(TriForElim);

%d verification
disp('Matlab,GNU/Octave built-in solution:  ');
disp(L\bL);

%Due to some unknown reason the solutions to Problem 1 appears at the very
%end of the document.
%% Problem 2

disp('Problem 2');
%a
%existing function works for multiple right-hand sides (RHS)
%b
%SimpleELim function
load testproblem.mat    
SimpleElimination= SimpleElim(A,b);

%Finding the Inverse of the matrix using Gausian Elimination
Inv=SimpleElim(A,eye(8));   
%Slicing the matrix into the A and A^(-1)
Inverse = Inv(:,9);
for ir=10:16
    Inverse=cat(2,Inverse,Inv(:,ir));
end

disp('Inverse/Gauss-Jordan Elimination:  ');
disp(Inverse);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(inv(A));
%% Problem 3
disp('Problem 3');
[Amod,Det]=Gauss_elim(A,b);
x2=backsub(Amod);

disp('Determinant/Gauss-Jordan Elimination:  ');
disp(Det);

disp('Matlab,GNU/Octave built-in solution:  ');
disp(det(A));
%% Functions
%Simple Forward Elimination
function [x] = SForwardElim(A,b,c)
%SForwardElim Simple Forward Elimination
% This function performs a simple forward elimination on a given square
% matrix and its RHS
%
% A= LHS matrix
% b= RHS matrix
% c= 0 if matrix needs to be catenated, and 1 if matrix is already catenated.
%The code used to write this function was written by Dr. Zettergen. Find it
%at https://github.com/Zettergren-Courses/EP501_matlab/blob/master/...
%linear_algebra/simple_elim_example.m

nref=length(b);             %size of system 
if c==0
    Awork=cat(2,A,b);           %Catenating A and B matrices
elseif c==1
    Awork=A;
else
    fprintf('Please provide 0 if the matrix needs to be catenated or 1 if catenation is not required')
end
for ir1=2:nref    %loop over rows from 2 to n performing elimination, this...
    %index marks what row we are starting the elimination from (i.e. using)...
    %for this particular column
    for ir2=ir1:nref   %this index marks the present position where...
        %elimination is being performed - i.e. where we are applying the...
        %elementary row operations
        fact=Awork(ir2,ir1-1); %multiplier of the variable we are...
        %attempting to eliminate, its ir-1 column of this row
        Awork(ir2,:)=Awork(ir2,:)-fact/Awork(ir1-1,ir1-1).*Awork(ir1-1,:);
        %subtract off previous row modified by a factor that eliminates the...
        %ir-1 column term in this row (so it has only super-diagonal elements)
    end %for
end %for
x=Awork;
end %function

%backsub function as provided by Dr. Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/...
%linear_algebra/backsub.m
function x=backsub(A)
% This function performs back substitution on an upper triangular matrix that has
% been modified by concatenating the RHS of the system.  
% Note that B is assumed to be upper triangular at this point.


n=size(A,1);                   %number of unknowns in the system
x=zeros(n,1);                  %space in which to store our solution vector
x(n)=A(n,n+1)/A(n,n);          %finalized solution for last variable,...
                                %resulting from upper triangular conversion

for ir1=n-1:-1:1
    x(ir1)=A(ir1,n+1);       %assume we're only dealing with a single...
                             %right-hand side here.
    fact=A(ir1,ir1);         %diagonal element to be divided through doing...
                             %subs for the ir2 row
    for ic=ir1+1:n
        x(ir1)=x(ir1)-A(ir1,ic)*x(ic);
    end %for
    x(ir1)=x(ir1)/fact;      %divide once at the end to minimize number of ops
end %for
end %function

%Forward Substituion for Lower Triangle
function  x = LTriFowardrSub(A,b)
% This function performs a forward elimination on a lower triangular
% matrix and its RHS
%
% A= LHS matrix
% b= RHS matrix
%The code used to write this function was written by Dr. Zettergen. Find it
%at https://github.com/Zettergren-Courses/EP501_matlab/blob/master/...
%linear_algebra/simple_elim_example.m
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

%Simple Backward Elimination
function [x] = SBackElim(A,b,c)
%SBackElim Simple Back Elimination
% This function performs a simple back elimination on a given square
% matrix and its RHS
%
% A= LHS matrix
% b= RHS matrix
% c= 0 if matrix needs to be catenated, and 1 if matrix is already catenated.
%The code used to write this function was modified from a code written by...
%Dr. Zettergen. Find it at
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/...
%linear_algebra/simple_elim_example.m

nref=length(b);                 %size of system
if c==0
    Awork=cat(2,A,b);           %Catenating A and B matrices
elseif c==1
    Awork=A;
else
    fprintf('Please provide 0 or 1 for third input')
end

for ir1=nref-1:-1:1         %loop over rows from 2 to n performing...
    %elimination, this index marks what row we are starting the elimination...
    %from (i.e. using) for this particular column
    for ir2=ir1:-1:1        %this index marks the present position where...
        %elimination is being performed - i.e. where we are applying the...
        %elementary row operations
        fact=Awork(ir2,ir1+1); %multiplier of the variable we are...
        %attempting to eliminate, its ir-1 column of this row
        Awork(ir2,:)=Awork(ir2,:)-fact/Awork(ir1+1,ir1+1).*Awork(ir1+1,:);    
        %subtract off previous row modified by a factor that eliminates the...
        %ir-1 column term in this row (so it has only super-diagonal elements)
    end %for
end %for
x=Awork;
end %function

%Simple Elimination/Gaussian Elimination
function [x] = SimpleElim(A,b)
%SimpleElim 
% This function performs both a simple back elimination and a simple 
%forward eliminatio on a given square matrix and its RHS
%
% A= LHS matrix
% b= RHS matrix
%The code used to write this function was modified from a code written by 
%Dr. Zettergen. Find it at
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/...
%linear_algebra/simple_elim_example.m

nref=length(b);             %Sizing the matrix
x = SForwardElim(A,b,0);    %Using Simple Forward Elimination function
x = SBackElim(x,b,1);       %Using Simple Back Elimination function
for ir=1:nref
    x(ir,:)=x(ir,:)/(x(ir,ir)); %Dividing through 
end
end %function

%Gauss elimination function as provided by Dr. Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/linear_algebra
function [Amod,Det]=Gauss_elim(A,b,verbose)

% [Amod,ord]=Gauss_elim(A,b,verbose)
%
% This function perform elimination with partial pivoting and scaling as
% described in Section 1.3.2 in the Hoffman textbook (viz. it does Gaussian
% elimination).  Note that the ordering which preserves upper triangularity
% is stored in the ord output variable, such that the upper triangular output
% is given by row-permuted matrix Amod(ord,:).  The verbose flag can be set to
% true or false (or omitted, default=false) in order to print out what the algirthm
% is doing for each elimination step.

%Parse the inputs, throw an error if something is obviously wrong with input data
narginchk(2,3);
if (nargin<3)
    verbose=false;
end %if

%Need to error check for square input.  

%Allocation of space and setup
Amod=cat(2,A,b);          %make a copy of A and modify with RHS of system
n=size(A,1);              %number of unknowns
ord=(1:n)';               %ord is a mapping from input row being operated...
            %upon to the actual row that represents in the matrix ordering
I=1;                      %Counter for number of pivots

%Elimination with scaled, partial pivoting for matrix Amod; note all row
%indices must be screen through ord mapping.
for ir1=1:n-1
    if (verbose)
        disp('Starting Gauss elimination from row:  ');
        disp(ir1);
        disp('Current state of matrix:  ');
        disp(Amod(ord,:));
    end %if
    
    %check scaled pivot elements to see if reordering should be done
    pivmax=0;
    ipivmax=ir1;      %max pivot element should never be higher than my...
                        %current position
    for ipiv=ir1:n    %look only below my current position in the matrix
        pivcurr=abs(Amod(ord(ipiv),ir1))/max(abs(Amod(ord(ipiv),:)));     
        %note that columns never get reordered
        if (pivcurr>pivmax)
            pivmax=pivcurr;
            ipivmax=ipiv;     %this stores the index into ord for row...
                                %having largest pivot element
            I=I+1;            %Counter for number of pivots
        end %if
    end %for
    
    %reorder if situation calls for it
    if (ipivmax ~= ir1)
        itmp=ord(ir1);
        ord(ir1)=ord(ipivmax);
        ord(ipivmax)=itmp;
        
        if (verbose)
            disp('Interchanging rows:  ');
            disp(itmp);
            disp(' and:  ');
            disp(ord(ir1));
            disp('Current matrix state after interchange:  ');
            disp(Amod(ord,:));
        end %if
    end %if
    
    %perform the elimination for this row, former references to ir1 are now
    %mapped through the ord array
    for ir2=ir1+1:n
        fact=Amod(ord(ir2),ir1);
        Amod(ord(ir2),ir1:n+1)=Amod(ord(ir2),ir1:n+1)-fact/...
            Amod(ord(ir1),ir1).*Amod(ord(ir1),ir1:n+1);    
        %only need columns ahead of where we are in matrix
    end %for
    
    if (verbose)
        disp('Following elimination for row:  ');
        disp(ir1);
        disp(' matrix state:  ');
        disp(Amod(ord,:));
    end %if
end %for
Amod=Amod(ord,:); %Reordering the Matrix
%Calculating the determinant
Det=Amod(1,1);
%Multiplying the diagonals
for ic=2:n
    Det=Det*Amod(ic,ic);
end
%Sign correction for pivoting
I=(-1)^I;
Det=Det*I;
end %function

