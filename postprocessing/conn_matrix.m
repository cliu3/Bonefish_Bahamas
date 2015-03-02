% constructing PDFs and connectivity matricies for an offline Lagrangian tracking model (or IBM)
%
% DESCRIPTION:
%    Read output particle data from LTRANS and construct PDFs
%
% INPUT
%
% OUTPUT:
%    PDFv connectivity
%
% Author(s):
%    Chang Liu (University of Massachusetts Dartmouth)
%
%
% Revision history
%
%==============================================================================
clear all; close all;

fprintf('initializing...\n');


% open the mesh file
nc = netcdf('../input/roms_grd_rot_raw.nc','nowrite');
lon_rho = nc{'lon_rho'}(:);
lat_rho = nc{'lat_rho'}(:);
mask_rho = logical(nc{'mask_rho'}(:));
h = nc{'h'}(:);
pm = nc{'pm'}(:);
pn = nc{'pn'}(:);
area=1./(pm.*pn);




% open the particle data (LTRANS potput)
fname = '../output/output.nc';

nc = netcdf(fname,'nowrite');
lonp = nc{'lon'}(:);
latp = nc{'lat'}(:);
age = nc{'age'}(:);
time=nc{'model_time'}(:);

ntime=numel(time);


% load spawning zone array
load ../preprocessing/SpawningZone.mat

% PROJECTION
[xm, ym]=baham_project(lon_rho,lat_rho,'forward');
[xp, yp]=baham_project(lonp,latp,'forward');


%% compute PDF
for source=1:max(SpawningZone);
    pt_idx = (SpawningZone==source);
    
    %for i=1:ntime
    i = ntime;
        
        fprintf('calculating PDF for time frame %d\n',i);
        
        PDF = probability(xp(i,pt_idx), yp(i,pt_idx), xm, ym);
        PDF(isnan(PDF))=0;
        
        
        PDFv{source,i} = PDF;
        
    %end
end


save PDFv  PDFv %connectivity;

% read in settlement zone files
set_files = {'Abaco_Settlement','Eleu_Settlement','GBI_Settlement'};
nsink=numel(set_files);



% mark destination zones on mesh
c = {};
for target=1:nsink
    clear settlement
    load(set_files{target})
    fprintf('target %d\n',target);
    [xv, yv]=baham_project(settlement.LON,settlement.LAT,'forward');

    % determine which nodes are within the settlement zone polygon
    
    c{target} = double(inpolygons(xm,ym,xv,yv));
    
end

%
% TO DO: test inpolygons
%



%% compute  connectivity

for source=1:max(SpawningZone);
    
    for target=1:nsink
        fprintf('calculating connectivity from %d to %d\n',source,target);
        
        % Integral
        settpdf = c{target}.*PDFv{source,end};
        connectivity(source,target) = dot(settpdf(:),area(:));
    end;
end



save('data','PDFv','connectivity')






