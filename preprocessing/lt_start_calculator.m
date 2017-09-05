function lt_start=lt_start_calculator(filenum,yyyy,mm,dd,hh,mi,ss)

% Calculate larval release time/LTRANS start time
%
% OUTPUT:
%   lt_start: seconds after LTRANS initialization (as in LTRANS.data)
% INPUTS:
%   filenum:  forcing file number (as in LTRANS.data)
%   yyyy,mm,dd,hh,mi,ss: 
%             date and time of release
%
% Example Usage: (for 1/26/2009 8:19 UTC release)
%    lt_start_calculator(08,2009,1,26,8,19,0)
%

nc_fname = ['../input/roms_his_rot_',num2str(filenum,'%02d'),'.nc'];
first_time = ncread(nc_fname, 'scrum_time',1,1,1);

first_time_datenum = double(datenum(2005,1,1,0,0,0) + first_time/86400.);

lt_start = floor((datenum(yyyy,mm,dd,hh,mi,ss) - first_time_datenum) * 86400);


end