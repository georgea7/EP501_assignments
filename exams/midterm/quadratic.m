function x = quadratic(A)
%Quadratic equation solver
%x = quadratic(A)
%   ax^2+bx+c=0
%A - [a,b,c]' column vector

x(1)=(-A(2)+sqrt(A(2)^2-4*A(1)*A(3)))/(2*A(1));
x(2)=(-A(2)-sqrt(A(2)^2-4*A(1)*A(3)))/(2*A(1));
end