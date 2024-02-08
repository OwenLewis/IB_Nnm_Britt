%
% circle.m
%
% return equal spaced points on a circle of a given radius and center
%
%  input,  xc,yc -- center of the circle
%           r    -- radius of the circle
%          ds    -- desired arclength spacing in a circular rest state
%  output, X     -- point locations size (Nib,2)  
%          newds -- actual acrlength spacing of points
%
function [X, newds] = circle(xc,yc,r,ds);
    
  % place the points
  %
  Nib = round(2*pi*r/ds);
  theta = [linspace(0,2*pi,Nib+1)]';
  X = [ xc + r*cos(theta), yc + r*sin(theta)];
  X = X(1:Nib,:);
      
  % adjust ds to reflect rounding
  %
  newds = (2*pi*r)/Nib;

