
% Generating initial position file for LTRANS
%
function initb()
%
% DESCRIPTION:
%    Dump initial particle positions and release (spawning) times to a csv file
%
% INPUT
%
% OUTPUT:
%    Initial particle position file
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Chang Liu
%
% Revision history
%
%==============================================================================
close all; clear all;

%=======================
% Define parameters
%=======================
nlag=12000;             %number of particles
% release date: many, but starting dec 08 2008 at 02:42 UTC-5
% 1st forcing file: roms_his_rot_dec_15.nc
% seconds from nov 30 2008 00:00 UTC to oct dec 08 2008 02:42 UTC-5: 718920
% release_beg=980339;        %time at which larval release begins (seconds after 1st forcing file)
release_beg=lt_start_calculator(08,2009,1,26,8,19,0);
release_end=release_beg;    %time at which larval release ends (seconds after 1st forcing file)
release_intv=0;    %time interval (in seconds) between two consecutive releases
fprintf('Release time (seconds after 1st forcing file, lt_start in LTRANS.data): %d\n',release_beg)

lag_pos_file = '../input/p12000_jan26.csv';

% Initialization of Variables
t_release = [];  %time of release in seconds after LTRANS initialization
lon_release = [];  %longitude of release location
lat_release = [];  %latitude of release location
z_release = [];  %depth of release location (in meters from surface, e.g., -35.55)

%------------------------------------------------------------------------------
% set a uniform distribution in a google-earth defined domain with
% release varying over
%------------------------------------------------------------------------------

% read in the kml file
%kml files = {'Abaco_Spawn.kml','Cape_Eleu_Spawn.kml','GBI_Spawn.kml'};
kml_files = {'ABACO_SPAWN_NEW.kml','ELEU_SPAWN_NEW.kml','GBI_SPAWN_NEW.kml','ANDROS_SPAWN_NEW.kml'};

latb=[];
lonb=[];
lat_release=[];
lon_release=[];

% divide total particles into (almost) equal groups
nspawnzone = numel(kml_files);
q=floor(nlag/nspawnzone);
r=mod(nlag,nspawnzone);
splits=[];
for i=1:nspawnzone
    splits(i)=q;
end
for i=1:r
    splits(i)=splits(i)+1;
end

for i=1:nspawnzone
    kml_file=kml_files(i);
%     kml_file{1}
    [latt,lont,dumz] = read_kml(kml_file{1});
    [lont, latt]=polycw(lont, latt);
    
    
    %     latb=[latb;latt;nan];
    %     lonb=[lonb;lont;nan];
    
    
    lat1 = min(latt);
    lat2 = max(latt);
    lon1 = min(lont);
    lon2 = max(lont);
    
    % create a uniform grid of points containing the kml perimeter
    % use twice the number of points because we want at least nlag inside the kml box
    % we will cull randomly later to get to nlag
    multifac = 2;
    ncull = -1;
    while(ncull < 0)
        lat_release0 = lat1 + (lat2-lat1).*rand(round(splits(i)*multifac),1);
        lon_release0 = lon1 + (lon2-lon1).*rand(round(splits(i)*multifac),1);
        
        % tag the points inside the box
        mark = inpolygons(lon_release0,lat_release0,lont,latt);
        pts = find(mark==1);
        
        
        % make sure we have enough
        ncull = length(pts) - splits(i);
        
        if (ncull < 0)
            multifac=multifac*1.1;
        end
    end;
    
    % report number
    fprintf('distributed %d points in the kml perimeter\n',length(pts))
    fprintf('culling %d points randomly\n',ncull);
    
    % randomly select points to remove
    pts = pts(randperm(end));
    % random_order = shuffle(1:length(pts));
    retain = pts(ncull+1:length(pts));
    
    lat_release0 = lat_release0(retain);
    lon_release0 = lon_release0(retain);
    
    lat_release=[lat_release;lat_release0];
    lon_release=[lon_release;lon_release0];
    
end





%plot release map
figure();hold on
plot(lonb,latb);
plot(lon_release,lat_release,'r.')
%------------------------------------------------------------------------
% set release time
%------------------------------------------------------------------------
% set release time using normal distribution
release_length=release_end-release_beg;
spawning_sigma=floor(release_length/6); %assuming 3-sigma 99%
if(release_beg ~= release_end);
    figure
    nrels = floor((release_end-release_beg)/release_intv);
    
    %params for normpdf()
    X=release_beg:release_intv:release_end;
    mu=(release_beg+release_end)*0.5;
    
    probs = normpdf(X,mu,spawning_sigma);
    probs = round(probs.*(nlag/sum(probs)));
    if(sum(probs)>nlag);
        if(probs(end)>(sum(probs)-nlag))
            probs(end)=probs(end)-(sum(probs)-nlag);
        else
            probs(end-1)=probs(end-1)-((sum(probs)-nlag)-probs(end));
            probs(end)=0;
        end
    end;
    plot(X,probs);
    xlabel('release time (seconds from start)');
    ylabel('number of particle released');
    
    % randomly select particles to be released on each day
    % issue is that randomly selecting integers from a small set of integers
    % can result in repeats, for example, randomly selecting five integers
    % from 1:10 can result in 1,4,2,4,5 being selected, thus we would
    % only be removing four integers from the set
    tvec = X(1:end-1);
    % shuffle a vector
    random_order = shuffle(1:nlag);
    i1 = 1;
    for ii=1:nrels
        if(probs(ii) > 0)
            i2 = i1 + probs(ii)-1;
            t_release(random_order(i1:i2)) = tvec(ii);
            i1 = i2 + 1;
        end
    end;
    % still some issue with shuffling, but minor, quick fix here
    pts = find(t_release<=0);
    if(~isempty(pts))
        t_release(pts) = mean(tvec);
    end;
    figure
    plot(t_release-min(t_release),1:nlag,'r+')
    %hist(t_release-min(t_release))
    %axis([release_beg,release_end,0,nlag]);
    xlabel('seconds after initial release')
    ylabel('particle # released on that day');
    
else %release start and end are the same, assign uniform t_release
    
    t_release = release_beg*ones(nlag,1);
    
end
%t_release = zeros(numel(lat_release));



%------------------------------------------------------------------------
% set release depth
% Release from 40-60 m or at the bottom if the depth is shallower
%------------------------------------------------------------------------

z_release = - (40+20.*rand(numel(lat_release),1));

%
% Save a zone number variable (i.e. which spawning zone a particle is from) as a mat file
% This will be used for computing connectivity matrix in postprocessing
% code
SpawningZone = zeros(nlag,1);
for i=1:numel(kml_files)
    [latt,lont,dumz] = read_kml(kml_files{i});
    in = inpolygon(lon_release,lat_release,lont,latt);
    SpawningZone(in) = i;

end
tabulate(SpawningZone)
fprintf('Each spawning zone should have (almost) equal number of particles.\n')
save SpawningZone SpawningZone

%------------------------------------------------------------------------------
% Output file format:
% the first column contains each particle?s latitudinal coordinate
% the second contains its longitudinal coordinate
% the third column contains the particle?s depth (in meters from surface, e.g., -35.55)
% the fourth column contains the 'date of birth' which is delay in seconds from the
%     beginning of the model run (1st forcing file) until the time that particle will start being tracked
%     (i.e., start to move).
%------------------------------------------------------------------------------
write_offline_lagfile_csv(lag_pos_file,lat_release, ...
    lon_release,z_release,t_release);
end

function [xout, yout]=polycw(xin, yin)
xin=xin(1:end-1);
yin=yin(1:end-1);

cx = mean(xin);
cy = mean(yin);

a = atan2(yin - cy, xin - cx);

[~, order] = sort(a);
order=flipud(order);
xout = xin(order);
yout = yin(order);

xout = [xout;xout(1)];
yout = [yout;yout(1)];

end

