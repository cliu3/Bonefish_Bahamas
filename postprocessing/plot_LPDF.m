% Plot LPDF
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
fname = '../output/output_16000.nc';

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


%% compute and plot LPDF
fprintf('Computing PDF...\n');
time_portion = [.01, .2, .4, .6,.8,1.0];  % portions of total time to calculate PDF

zoom_axis=[-80, -74,22,28.5];
fig1=figure('Position',[100,100,1500,300]);
fig2=figure('Position',[100,100,1800,300]);

load /research/data_utils/MyColormaps2.mat

for source=1:max(SpawningZone);
    pt_idx = (SpawningZone==source);
    
    for i=1:numel(time_portion)
        time_idx = ceil(time_portion(i).*ntime);
        %i = ntime;
        
        fprintf('calculating PDF for time frame %d/%d\n',time_idx,ntime);
        
        PDF = probability(xp(time_idx,pt_idx), yp(time_idx,pt_idx), xm, ym);
        PDF(isnan(PDF))=0;
        PDF = PDF ./ dot(PDF(:),area(:));   % NORMALIZE
        
        PDFv{source,i} = PDF;
        
        
        if source==2
            figure(fig1)
            subplot_tight(1,6,i);
            field=h;
            field(mask_rho <1)=nan;
            %pcolor(lon_rho,lat_rho,field);axis equal
            contour(lon_rho,lat_rho,mask_rho,'color',[.5 .5 .5]);axis equal
            hold on
            %         contour(lon_rho,lat_rho,mask_rho,[1],'k');
            %         hold on
            plot(lonp(time_idx,pt_idx), latp(time_idx,pt_idx),'r.','markersize',2)
            shading interp
            axis(zoom_axis);
            title(['day ',num2str(floor(time_idx/24)),'/60'])
            if i>1
                set(gca,'YTickLabel','')
            end
            
            figure(fig2)
            subplot_tight(1,7,i);
            field=PDF;
            field(mask_rho <1)=nan;
            pcolor(lon_rho,lat_rho,field*1e6);shading interp;axis equal
            hold on
            axis(zoom_axis)
            contour(lon_rho,lat_rho,mask_rho,'color',[.5 .5 .5])
            %set(gca,'Color',[220,156,26]/255)
            colormap(mycmap)
            caxis([0 6e-4]);
            title(['day ',num2str(floor(time_idx/24)),'/60'])
            if i>1
                set(gca,'YTickLabel','')
                
            end
        end
    end
end
figure(fig1)
set(gcf,'Color','w')
export_fig(['particles.png'])

figure(fig2)
subplot_tight(1,7,7);
colormap(mycmap)
caxis([0 6e-4]);axis equal;axis off;
cbar=colorbar('location','West');
initpos = get(cbar,'Position');
initpos(4)=initpos(4)*0.95;
set(cbar,'Position',initpos)
set(gcf,'Color','w')
export_fig(['LPDF.png'])





