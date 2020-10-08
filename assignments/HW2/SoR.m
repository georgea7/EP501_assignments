function [x,nit]=SoR(x0,A,b,tol,verbose,w)
%Successive over-Relaxation
%This code has been developed by modifying the Jacobi.m code provided by
%Dr.Zettergen in
%https://github.com/georgea7/EP501_matlab/blob/master/linear_algebra/Jacobi.m
%
%Inputs
%x0     : initial guess
%A      : LHS Matrix
%b      : RHS Matrix
%tol    : Tolerance
%w      : Relaxation parameter (Omega)
%verbose: 'true' to see steps, 'false' to hide steps
%
%Outputs
%x      : Solution
%nit    : Number of iterations
%% Check the inputs
narginchk(3,6);
if nargin<4
    tol=1e-6;
end %if
if nargin<5
    verbose=false;
end %if


%% Setup iterations
maxit=100;    %max number of iterations
n=size(A,1);  %system size
residual=10*ones(n,1);
difftot=1e3+tol;   %max sure we enter iterations
x=x0;


%% Perform iterations
it=1;
while(difftot>tol && it<=maxit)
    difftotprev=difftot;
    resprev=residual;
    xprev=x;
    for i=1:n
        residual(i)=b(i);
        for j=1:i-1
            residual(i)=residual(i)-A(i,j)*x(j);
        end %for
        for j=i:n
            residual(i)=residual(i)-A(i,j)*xprev(j);
        end %for
        x(i)=xprev(i)+w*residual(i)/A(i,i);
    end %for
    difftot=sum(abs(residual-resprev));
    
    if (verbose)
        fprintf('x= ');
        for i=1:n
            fprintf('%f   ',x(i));
        end %for
        fprintf('\n');
        fprintf('it=%d; difftot = %e\n',it,difftot);
    end %if
    
    if (difftot>difftotprev & it>2)
        error('Solution appears to be diverging, check diagonal dominance...')
    end %if
    it=it+1;
end %while

nit=it-1;
if (nit==maxit)
    warning('Solution may not have converged fully...')
end %if

end %function
