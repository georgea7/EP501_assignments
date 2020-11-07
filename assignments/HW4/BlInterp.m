function finterpmanual=BlInterp(xg,yg,f,xi,yi)
%f=BlInterp(xg,yg,f,x,y)
%Output
%   f: Interpolated values
%Input
%   xg: X grid
%   yg: Y grid
%    f: function
%    x: x values
%    y: y values

%Using example code provided by Dr. Zettergen

n=length(xi);               %Length of inputs
f2D=reshape(f,[96 96]);     %Reshaping f2D
finterpmanual=zeros(n,1);   %alotting space
for k=1:n                   %running for number of inputs
    x1=xi(k);               %kth x input
    y1=yi(k);               %kth y input
    [i,j]=InterpIndex2D(xg,yg,xi(k),yi(k));  %x and y indices 
    x=[xg(i),xg(i+1)];                       %xi, xi+1
    y=[yg(j),yg(j+1)];                       %yj, yj+1
    f=[f2D(i,j), f2D(i+1,j); f2D(i,j+1), f2D(i+1,j+1)]; %function values
    
    [X,Y]=meshgrid(x,y);
    fvec=f(:);
    xvec=X(:);
    yvec=Y(:);
    M=[ones(4,1),xvec(:),yvec(:),xvec(:).*yvec(:)];
    [Mmod,order]=Gauss_elim(M,fvec);
    avec=backsub(Mmod(order,:));
    finterpmanual(k)=avec(1)+avec(2)*x1+avec(3)*y1+avec(4)*x1*y1;
end
end