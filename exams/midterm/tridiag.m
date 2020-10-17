function x = tridiag(A,b)
%Tridiagonal Matrix Algorith or Thomas Algorithm
%x = tridiag(A,b)
%A - tridiagonal matrix
%b - b matrix

%Transforming A matrix to [3 x n] format
n=length(b);               %length of matrix
nullspce=0;                     %empty space to add to first column of Aprime
Aprime(:,1)=diag(A,-1);    %Aprime first coloumn
Aprime=cat(1,nullspce,Aprime);  %Correcting rows since first column has empty first row 
Aprime(:,2)=diag(A);       %Aprime second column
Aprime3=diag(A,1);         %Aprime thrid column done seperately
Aprime3(n,1)=0;                 %Correcting rows since third column has empty last row 
Aprime=cat(2,Aprime,Aprime3);   %Aprime third coloumn

%forward sweep 
for i=2:n
    proc=Aprime(i,1)/Aprime(i-1,2);
    Aprime(i,1)=proc;
    Aprime(i,2)=Aprime(i,2)-proc*Aprime(i-1,3);
    b(i)=b(i)-Aprime(i,1)*b(i-1);
end %for

%back substitution
x(n)=b(n)/Aprime(n,2);
for i=n-1:-1:1
    x(i)=(b(i)-Aprime(1,3)*x(i+1))/Aprime(i,2);
end
x=x';
end

    
