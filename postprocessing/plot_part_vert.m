% plot verticle view of a single particle (time vs temp against background vertical temperature)

clear all; close all;

% paremeters
partID = 9926;  % which particle to plot

% open the mesh file
nc = netcdf('../input/roms_grd_rot_raw.nc','nowrite');
lon_rho = nc{'lon_rho'}(:);
lat_rho = nc{'lat_rho'}(:);
mask_rho = logical(nc{'mask_rho'}(:));
h = nc{'h'}(:);

% read vertical grid
nc = netcdf('../input/roms_his_rot_1.nc','nowrite');
sc_r = nc.sc_r(:);
Cs_r = nc.Cs_r(:);
hc = nc.hc(:);


% open the particle data (LTRANS output)
fname = '../output/output_16000.nc';

nc = netcdf(fname,'nowrite');
lonp = nc{'lon'}(:,partID);
latp = nc{'lat'}(:,partID);
zp = nc{'depth'}(:,partID);
age = nc{'age'}(:,partID);
time=nc{'model_time'}(:);

ntime=numel(time);

% finding nearest rho grid points along the track
for i=1:ntime
    tmp = sqrt((lon_rho(:)-lonp(i)).^2+(lat_rho(:)-latp(i)).^2);
    nearest(i)=find( tmp==min(tmp)  );
    
end
[iind,jind]=ind2sub(size(lon_rho),nearest);

bathp = h(nearest);
vertz = repmat((hc.*sc_r)',[1,numel(bathp)]) + ((bathp(:)-hc)*Cs_r)';
xx = repmat(1:ntime,[length(sc_r),1]);

% read forcing file
forcing_file='../input/roms_his_rot_1.nc';
nc = netcdf(forcing_file,'nowrite');
scrum_time = nc{'scrum_time'}(:);
time_comp=time+127440000;


temp_plot=[];
for i=1:ntime
    [mintimediff,timeidx] = min(abs(scrum_time-time_comp(i)));
    if (mintimediff>3600)
        forcing_file(end-3)=num2str(str2num(forcing_file(end-3))+1);
        nc = netcdf(forcing_file,'nowrite');
        scrum_time = nc{'scrum_time'}(:);
        [mintimediff,timeidx] = min(abs(scrum_time-time_comp(i)));
    end

    
    temp = nc{'temp'}(timeidx,:,iind(i),jind(i));
    temp_plot=[temp_plot,temp];
end

figure('Position',[100,100,900,400]);
pcolor(xx,vertz,temp_plot);shading interp;hold on
colorbar;
plot(zp,'k-','linewidth',2)
set(gcf,'PaperPositionMode','auto')
print(['Part_vert_',num2str(partID)],'-dpng','-r0')
export_fig(['Part_vert_',num2str(partID),'.png'])