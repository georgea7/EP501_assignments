function [I] = InterpIndex1D(xg,x)
%[i] = InterpInd(xg,x)
%Output
%   i : Index
%   xg: grid
%   x : Input values
%This function identifies index i for x, where xgi<x<xgi+1

%Data Processing
n=length(x);
I=zeros(n,1);
for j=1:n
    %Initialisation
    fd=1;
    i=1;
    while fd>abs(xg(2)-xg(1))
        i=i+1;
        fd=abs(xg(i)-x(j));
    end
    I(j)=i;
end
end%function