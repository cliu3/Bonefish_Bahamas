function write_offline_lagfile_csv(file,lat,lon,z,trel)
% write position and release time data for particles to initialize offline Lagrangian code
%
% function write_offline_lagfile_csv(file,lat,lon,z,trel)
%
% DESCRIPTION:
%   write position and release time data for particles to initialize 
%   offline Lagrangian code
%
% INPUT
%   file:  filename
%   info:  string containing information about case (a global attribute)
%   x:     x-coordinate of initial particle position
%   y:     y-coordinate of initial particle position
%   v:     vertical-coordinate of initial particle position
%   vtype: vertical coordinate type (='z' for z-coordinate, ='s' for sigma coordinate)
%
% OUTPUT:
%    Initial particle position file
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

fprintf('creating initial lagrangian position file\n')

%------------------------------------------------------------------------------
% Output file format:
% the first column contains each particle?s longitudinal coordinate
% the second contains its latitudinal coordinate
% the third column contains the particle's depth (in meters from surface, e.g., -35.55)
% the fourth column contains the 'date of birth' which is delay in seconds from the 
%     beginning of the model run until the time that particle will start being tracked 
%     (i.e., start to move).
%------------------------------------------------------------------------------

% open boundary forcing
fid = fopen(file,'w');

for i=1:numel(lat)
  fprintf(fid,'%14.11f, \t %14.11f, \t %f, \t %f\n',lon(i),lat(i),z(i),trel(i));
end;
fclose(fid);

