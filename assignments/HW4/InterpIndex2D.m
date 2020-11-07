function [i,j] = InterpIndex2D(xg,yg,x,y)
%[i] = InterpInd(xg,x)
%Output
%   i : Index for x direction
%   j : Index for y direction
%   xg: grid in x direction
%   yg: grid in y direction
%   x : Input values
%This function identifies index i for x, where xgi<x<xgi+1 and index j for
%ygj<x<ygj+1

i=InterpIndex1D(xg,x);
j=InterpIndex1D(yg,y);

end%function