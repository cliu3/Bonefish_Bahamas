% function probability 
%
% DESCRIPTION:
%    Calculate probability density function of a set of particles,
%    interpolated to grid points
%
% INPUT
%    xp, yp      : cordinates of data points
%    xm, ym      : mesh parameters from FVCOM grid
% OUTPUT:
%    PDF         : probability density values interpolated to grid
%
% Author(s):  
%    Chang Liu (University of Massachusetts Dartmouth)
%    
% Revision history
%   
%==============================================================================
function PDF = probability(xp, yp, xm, ym)

% compute PDF
xr=max(xp)-min(xp);
yr=max(yp)-min(yp);

% make the computation area square
a=0.25*xr;
b=0.75*xr-0.5*yr;
MIN_XY = [min(xp)-a,min(yp)-b];
MAX_XY = [max(xp)+a,max(yp)+b];

data = [xp;yp]';
[~,density,X,Y]=kde2d(data,128,MIN_XY,MAX_XY);
PDF = interp2(X,Y,density,xm,ym);



end
