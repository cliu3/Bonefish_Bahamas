clear all; close all;


case_names = {'bahamas_test_apr10', 'bahamas_test_dec28', 'bahamas_test_dec8', 'bahamas_test_feb15', 'bahamas_test_feb8', 'bahamas_test_jan26', 'bahamas_test_mar11', 'bahamas_test_nov26', 'bahamas_test_nov4'};

%% Initialize empty cell array for total LPDF
PDFv_total=zeros(4,293,696);


for case_name=case_names
    
    load(['../',case_name{1},'/postprocessing/PDFv.mat']);
    for i=1:4
        PDFv_total(i,:,:) = squeeze(PDFv_total(i,:,:)) + PDFv{i,end};
    end

end




% for ii=1:4
%     for jj=926:1645
%         PDFv_total{ii,jj}=zeros(293,696);
%     end
% end
% 
% % Summation
% for i=1:4
%     load(['PDFv',num2str(i),'.mat']);
%     for ii=1:4
%         for jj=926:1645
%             PDFv_total{ii,jj}=PDFv_total{ii,jj}+PDFv{ii,jj};
%         end
%     end
% end

%save('PDFv_pooled', 'PDFv_total','-v7.3') 
%% 


%plot LPDF
%roms_mesh   = ['~/Documents/LTRANSv2b/run/Bahamas_test/input/roms_grd_rot_raw.nc'];
roms_mesh   = ['../bahamas_test_apr10/input/roms_grd_rot_raw.nc'];
mask_rho = ncread(roms_mesh,'mask_rho')';
mask_rho = logical(mask_rho);
lat_rho = ncread(roms_mesh,'lat_rho')';
lon_rho = ncread(roms_mesh,'lon_rho')';

spawning_names = {'Abaco','Eleu','GBI','Andros'};
for source=1:4
    figure()
    PDF_plot = squeeze(PDFv_total(source,:,:));
    PDF_plot(~mask_rho)=nan;
    pcolor(lon_rho,lat_rho,PDF_plot)
    caxis([0,0.0000000005])
    shading flat
    colorbar;
    axis([-79.5 -75.5 24 27.5])
    title(['LPDF of individuals spawned in ',spawning_names{source}])
    saveas(gcf,['pooled_LPDF_',spawning_names{source},'.png'])
end

%% pooled connectivity
pooled_connectivity = zeros(4,4);

for case_name=case_names
    
    load(['../',case_name{1},'/postprocessing/connectivity.mat']);
    
    pooled_connectivity=pooled_connectivity+connectivity;

end

% normalize
pooled_connectivity=pooled_connectivity./max(pooled_connectivity(:));

%plot connectivity
fig1=figure('Position',[100,100,400,300]);

imagesc(pooled_connectivity);
colorbar;
axis xy

tick=1:4;
ticklabel={'Abaco','Eleu','GBI','Andros'};
set(gca,'XTick',tick);
set(gca,'YTick',tick);
set(gca,'XTickLabel',ticklabel);
set(gca,'YTickLabel',ticklabel);

set(gcf,'Color','w')

saveas(gcf,'pooled_connectivity.png')
