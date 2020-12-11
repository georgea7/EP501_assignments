function x=backsub(A)

% This function performs back substitution on an upper triangular matrix that has
% been modified by concatenating the RHS of the system.  
% Note that B is assumed to be upper triangular at this point.


n=size(A,1);                   %number of unknowns in the system
x=zeros(n,1);                  %space in which to store our solution vector
x(n)=A(n,n+1)/A(n,n);          %finalized solution for last variable, resulting from upper triangular conversion

for ir1=n-1:-1:1
    x(ir1)=A(ir1,n+1);       %assume we're only dealing with a single right-hand side here.
    fact=A(ir1,ir1);         %diagonal element to be divided through doing subs for the ir2 row
    for ic=ir1+1:n
        x(ir1)=x(ir1)-A(ir1,ic)*x(ic);
    end %for
    x(ir1)=x(ir1)/fact;      %divide once at the end to minimize number of ops
end %for

end %function

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

function fval=fRK(t,y,alpha)
   fval=-alpha*y;
end %function

function [Amod,ord]=Gauss_elim(A,b,verbose)

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
ord=(1:n)';               %ord is a mapping from input row being operated upon to the actual row that represents in the matrix ordering

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
    ipivmax=ir1;      %max pivot element should never be higher than my current position
    for ipiv=ir1:n    %look only below my current position in the matrix
        pivcurr=abs(Amod(ord(ipiv),ir1))/max(abs(Amod(ord(ipiv),:)));      %note that columns never get reordered...
        if (pivcurr>pivmax)
            pivmax=pivcurr;
            ipivmax=ipiv;     %this stores the index into ord for row having largest pivot element
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
        Amod(ord(ir2),ir1:n+1)=Amod(ord(ir2),ir1:n+1)-fact/Amod(ord(ir1),ir1).*Amod(ord(ir1),ir1:n+1);    %only need columns ahead of where we are in matrix
    end %for
    
    if (verbose)
        disp('Following elimination for row:  ');
        disp(ir1);
        disp(' matrix state:  ');
        disp(Amod(ord,:));
    end %if
end %for

end %function

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

function [root,it,success]=newton_exact(f,fprime,x0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
%
% finds a set of roots corresponding to the function f (input as a handle)
% given a function which computes the derivative

%% Error checking of input and setting of default values
narginchk(3,6);   %check for correct number of inputs to function
if (nargin<4)
    maxit=100;       %maximum number of iterations allowed
end %if
if (nargin<5)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<6)
    verbose=false;
end %if


%% Make sure we don't start at an inflection point with zero derivative
if (abs(fprime(x0))<tol)
    warning(' Attempting to start Newton iterations near an inflection point, you may wish to restart with a different guess...');
    x0=x0+1;   %bump the guess a ways off of initial value to see if we can get anything sensible
end %if


%% Newton iterations
it=1;
root=x0;
fval=f(root);
converged=false;
while(~converged && it<=maxit)
    derivative=fprime(root);
    if (abs(derivative)<100*tol)    %this (inflection point) will end up kicking the root really far away...
        converged=false;
        warning(' Derivative close to zero, terminating iterations with failed convergence... ');
        break;
    else
        root=root-fval./derivative;    % update root estimate
        fval=f(root);                  % see how far off we are from zero...
        if (verbose)
            fprintf(' iteration: %d; root:  %f + %f i; function value: %f, derivative:  %f \n',it,real(root),imag(root),fval,derivative);
        end %if
        it=it+1;
        converged=abs(fval)<tol;
    end %if
end %while
it=it-1;

if (~converged)
    warning('Used max number of iterations, or derivative near zero...')
    success=false;
else
    success=true;
end %if

end %function

function y=objfun1(x)

y=(32/3)*x.^4-(32/3)*x.^3+(8)*x.^2-4*x;

end %function

function y=objfun2(x)

y=(32/3)*x.^4-(32/3)*x.^3+(8)*x.^2-4*x+2;

end %function

function y=objgradfun1(x)

y=(128/3)*x.^3-(32)*x.^2+(32)*x-4;

end %function