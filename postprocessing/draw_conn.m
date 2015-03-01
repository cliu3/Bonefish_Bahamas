clear all; close all;
addpath('/research/data_utils')


load('MyColormaps','mycmap')
%%
load data.mat
load S.mat

close all
imagesc(mean_connectivity);
hold on
axis xy;
axis square;
colorbar;
set(gcf,'Colormap',mycmap)
set(gca,'XTick',[1:27])
set(gca,'YTick',[1:25])
for i=[7.5, 14.5,23.5,25.5]
    plot(0.5:1:27.5,i*ones(1,28),'k--')
    plot(i*ones(1,26),0.5:1:25.5,'k--')
end
colorbar
caxis([0 0.25])
title('total connectivity for 2008 early case')

set(gcf,'color','white')
export_fig('total_connectivity.pdf')

figure
imagesc(S_hourly);
hold on
axis xy;
axis square;
colorbar;
set(gcf,'Colormap',mycmap)
set(gca,'XTick',[1:27])
set(gca,'YTick',[1:25])
for i=[7.5, 14.5,23.5,25.5]
    plot(0.5:1:27.5,i*ones(1,28),'k--')
    plot(i*ones(1,26),0.5:1:25.5,'k--')
end
colorbar
caxis([0 0.351])
title('transport success for 2008 early case')

set(gcf,'color','white')
export_fig('transport_success.pdf')