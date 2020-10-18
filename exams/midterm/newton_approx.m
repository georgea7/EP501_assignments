function [root,it,success]=newton_approx(f,x0,x0i1,maxit,tol,verbose)

% root=newton_approx(f)
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

%% Newton iterations
it=1;               %initiation
root=x0;            %Xi
rooti1=x0i1;        %Xi-1
fval=f(root);       %initiation
fvali1=f(rooti1);   %initiation
converged=false;
while(~converged && it<=maxit)
    derivative=(fval-fvali1)/(root-rooti1); %secant method formula
    if (abs(derivative)<100*tol)    %this (inflection point) will end up kicking the root really far away...
        converged=false;
        warning(' Derivative close to zero, terminating iterations with failed convergence... ');
        break;
    else
        rooti1=root;                   % update root_i-1 
        root=root-fval./derivative;    % update root estimate
        fval=f(root);                  % see how far off we are from zero...
        fvali1=f(rooti1);              % update fval_i-1
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