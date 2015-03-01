function [out_east,out_north] = baham_project(in_east,in_north,direction)

[out_east,out_north] = sp_proj('0901',direction,in_east,in_north,'m');

end