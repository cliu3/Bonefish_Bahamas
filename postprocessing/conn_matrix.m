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
addpath('../preprocessing')
fprintf('initializing...\n');

%dropbox_dir = '/Users/gcowles/Dropbox/Bonefish_Bahamas_Project';

roms_mesh   = ['../input/roms_grd_rot_raw.nc'];
LTRANS_output = ['../output/output.nc'];
spawning_zones = ['../preprocessing/SpawningZone.mat'];

%set_files = {'Abaco_Settlement','Eleu_Settlement','GBI_Settlement'};

set_files = {'Abaco_settlement_final_wgs84', ...
             'Eleu_settlement_final_wgs84', ...
             'GBI_Settlement_final_wgs84', ...
             'Andros_settlement_final_wgs84'};

% open the mesh file
h  = ncread(roms_mesh,'h')';
pm = ncread(roms_mesh,'pm')';
pn = ncread(roms_mesh,'pn')';
area = 1./(pm.*pn);
lat_rho = ncread(roms_mesh,'lat_rho')';
lon_rho = ncread(roms_mesh,'lon_rho')';
mask_rho = ncread(roms_mesh,'mask_rho')';
mask_rho = logical(mask_rho);

% lon_rho = nc{'lon_rho'}(:);
% lat_rho = nc{'lat_rho'}(:);
% mask_rho = logical(nc{'mask_rho'}(:));
% h = nc{'h'}(:);
% pm = nc{'pm'}(:);
% pn = nc{'pn'}(:);
% area=1./(pm.*pn);




% open the particle data (LTRANS potput)

% nc = netcdf(LTRANS_output,'nowrite');
% lonp = nc{'lon'}(:);
% latp = nc{'lat'}(:);
% age = nc{'age'}(:);
% time=nc{'model_time'}(:);
%nc = ncread(LTRANS_output,'nowrite');
lonp = ncread(LTRANS_output,'lon')';  
latp = ncread(LTRANS_output,'lat')';  
age  = ncread(LTRANS_output,'age')';
time = ncread(LTRANS_output,'model_time')';
dob  = ncread(LTRANS_output,'dob')';
age  = ncread(LTRANS_output,'age')';
age  = age(:,1)/(24.*3600); %age in days

if(numel(unique(dob)) ~= 1)
  error('not all particles were released simultaneously');
end;

ntime=numel(time);
npart=size(lonp,2);

% setup the settlement probability distribution
sett_beg = 40;  %start at day 40
sett_end = 72;  %end at day 72
settprob = zeros(numel(age),1);
[mini,ibeg] = min(abs(age-sett_beg)); %output index marking settlement begin
[mini,iend] = min(abs(age-sett_end)); %output index marking settlement end
settprob(ibeg:iend) = 1./(iend-ibeg+1);
plot(age,settprob);
xlabel('time since release (days)');
ylabel('settlement probability (-/day)');
fprintf('sum of settlement probability should be 1, is %f\n',sum(settprob));




% load spawning zone array
load(spawning_zones);

% PROJECTION
[xm, ym]=baham_project(lon_rho,lat_rho,'forward');
[xp, yp]=baham_project(lonp,latp,'forward');

% specify the range of particle age based upon which the LPDF is computed
% e.g., to make the LPDF only replect probabilities from 40-72 days post-release,
% set age0 = 40 and age1 = 72
age0 = sett_beg;  % number in days
age1 = sett_end;  % number in days
time_idx = find(time >= (age0-1)*86400 & time <= (age1)*86400);


%% compute PDF
for source=1:max(SpawningZone);
    pt_idx = (SpawningZone(1:npart)==source);
    for i = time_idx;
        fprintf('calculating PDF for time frame %d\n',i); 
        PDF = probability(xp(i,pt_idx), yp(i,pt_idx), xm, ym);
        PDF(isnan(PDF))=0;
        PDF = PDF ./ dot(PDF(:),area(:));   % NORMALIZE
        PDFv{source,i} = PDF;
    end;  
end

save('PDFv',  'PDFv','-v7.3') %connectivity;

% read in settlement zone files
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





%% compute  connectivity

for source=1:max(SpawningZone);
    
    for target=1:nsink
        fprintf('calculating connectivity from %d to %d\n',source,target);
        settpdf = xm*0; 
        
        % Loop over output frames, accumulating PDF based on settlement
        % prob
        for i = time_idx;  
          if(settprob(i) > 0);
            fprintf('including PDF from frame %f with probability %f\n',i,settprob(i));
            settpdf = settpdf + c{target}.*PDFv{source,i}.*settprob(i);
          end;
        end;

        connectivity(source,target) = dot(settpdf(:),area(:));
    end;
end



save('data','PDFv','-v7.3','connectivity')

%% make some plots
%plot LPDF

spawning_names = {'Abaco','Eleu','GBI','Andros'};
for source=1:4;
    figure()
    PDF_plot = PDFv{source,end};
    PDF_plot(~mask_rho)=nan;
    pcolor(lon_rho,lat_rho,PDF_plot)
    shading flat
    colorbar;
    axis([-79.5 -75.5 24 27.5])
    title(['LPDF of individuals spawned in ',spawning_names{source}])
    saveas(gcf,['LPDF_',spawning_names{source},'.png'])
end

%plot connectivity
fig1=figure('Position',[100,100,400,300]);

imagesc(connectivity);
colorbar;
axis xy

tick=1:4;
ticklabel={'Abaco','Eleu','GBI','Andros'};
set(gca,'XTick',tick);
set(gca,'YTick',tick);
set(gca,'XTickLabel',ticklabel);
set(gca,'YTickLabel',ticklabel);

set(gcf,'Color','w')

saveas(gcf,'connectivity.png')


