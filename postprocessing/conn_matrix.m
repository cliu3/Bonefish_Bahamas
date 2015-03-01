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

% PROJECTION
[xm, ym]=baham_project(lon_rho,lat_rho,'forward');
[xp, yp]=baham_project(lonp,latp,'forward');


%% compute PDF

for i=1:ntime
    
    fprintf('calculating PDF for time frame %d\n',i);
    
    
    
    PDF = probability(xp(i,:), yp(i,:), xm, ym);
    PDF(isnan(PDF))=0;
    
    
    PDFv{i} = PDF;
    
end


save PDFv  PDFv %connectivity;


%% ====== below code to be adapted ========
  

%% compute  connectivity

for source=1:nsource;

    
    for target=1:nsink
        fprintf('calculating connectivity from %d to %d\n',source,target);
        
        % Integration
        
        mean_connectivity(source,target) = dot(c{target}.*wmPDF(source,:)',art1);
    end;
    
end


%% compute settlement success

deltaT = age(2);
S = zeros(nsource,nsink);
% compute hourly settlememnt coefficient

D = 0.1849*ones(size(age));
D(age>=14) = 0.1849.*erfc(0.4*(age(age>=14)-14));


xp=[];
yp=[];

%% hourly
for source=1:nsource;
    [~, nlag] = size(xp0(1,tid0==source));
    for target=1:nsink
        fprintf('calculating settlement success from %d to %d\n',source,target);
        % select settlement zone
        xv=settlement(target).X;
        yv=settlement(target).Y;
        r=age(age>=10);
        for t=r';
            fprintf('calculating age: %d\n',t);
            xp=xp0(age==t,tid0==source);
            yp=yp0(age==t,tid0==source);
            
            E = inpolygon(xp,yp,xv,yv);
            
            
            S(source,target) = S(source,target) + sum(D(age==t).*E.*deltaT);
            
        end
        S(source,target) = S(source,target)./nlag;
        
    end
end
S_hourly = S;
toc

%% daily
tic
for source=1:nsource;
    [~, nlag] = size(xp0(1,tid0==source));
    for target=1:nsink
        fprintf('calculating settlement success from %d to %d\n',source,target);
        % select settlement zone
        xv=settlement(target).X;
        yv=settlement(target).Y;
        
        for t=10:19;
            fprintf('calculating age: %d\n',t);
            xp=xp0(age==t,tid0==source);
            yp=yp0(age==t,tid0==source);
            
            E = inpolygon(xp,yp,xv,yv);
            
            
            S(source,target) = S(source,target) + sum(D(age==t).*E);
            
        end
        S(source,target) = S(source,target)./nlag;
        
    end
end
S_daily = S;
toc

save S S_hourly S_daily

save('data','wmPDF','mean_connectivity','-append')






