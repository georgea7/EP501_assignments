function [root,it,success]=newton_approx(f,x01,x0i1,maxit,tol,verbose)

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


% %% Make sure we don't start at an inflection point with zero derivative
% if (abs(fprime(x0))<tol)
%     warning(' Attempting to start Newton iterations near an inflection point, you may wish to restart with a different guess...');
%     x0=x0+1;   %bump the guess a ways off of initial value to see if we can get anything sensible
% end %if


%% Newton iterations
it=1;
root=x01;       %Xi
rooti1=x0i1;      %Xi-1
fval=f(root);
fvali1=f(rooti1);
converged=false;
while(~converged && it<=maxit)
    derivative=(fval-fvali1)/(root-rooti1); %secant method formula
    if (abs(derivative)<100*tol)    %this (inflection point) will end up kicking the root really far away...
        converged=false;
        warning(' Derivative close to zero, terminating iterations with failed convergence... ');
        break;
    else
        rooti1=root;
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