% gf visum images 
% load
clear all
close all
clc
% load('/data/Technion_analysis/Amygdala/Vis_dimorphisim/adjpos.mat')
load('/data/Technion_analysis/goldfish/visium/gf_vis_alt/Sorted.mat','sample_sorted') % vdata
% load('/data/Technion_analysis/goldfish/visium/cropped/040122_dbscan_knn_n/Orignal.mat', 'geneid_all')%vgenes
% v_gen=geneid_all;
% v_data=data_sorted_all;
% load('/data/Technion_analysis/goldfish/g_matrix.mat','g_geneid','g_cluster_name','g_data','g_ca') % cgenes
% load('/data/Technion_analysis/goldfish/gf_mean_data1.mat','mean_data1') % cdata
% c_all_data=g_data;
% c_data=mean_data1;
% c_gen=g_geneid;
% clusteruni=g_cluster_name;
all_vis_id=natsort(unique(sample_sorted));
all_vis_id(5)=[];
set(0,'DefaultFigureVisible','on');

%% images
% namepdf='all_DAPI';
% cd('/data/runs/samples/Vis_collection')
f=figure('Name','all_dapis','NumberTitle','off','color','k');
ha = @(m,n,p) subtightplot (m, n, p,[.01 .01],[.01 .01],[.01 .01]);
% ax_link=[];
for vi=1:length(all_vis_id)
    vi
    imgv=imread(['/data/runs/samples/Vis_collection/',char(all_vis_id(vi)),'/',char(all_vis_id(vi)),'.tif']);
%     ax(vi)=ha(1,length(all_vis_id),vi);
ax(vi)=ha(3,8,vi);
    if vi>8
%     imgv=flip(img,2);
    imgx=imrotate(imadjust(imgv),90);
    imgx=flip(imgx,2);
    imshow(imgx,[])
    else
    imgx=imrotate(imadjust(imgv),90);
    imgx=flip(imgx,1);
    imshow(imgx,[])
    end
%      camroll(90) 
%      set(gca, 'YDir','reverse')
%     end
%     title(all_vis_id(vi))
%     pause
% cd('/data/Technion_analysis/goldfish/visium/dapi_images')
% imwrite(im2uint8(imgx),['/data/Technion_analysis/goldfish/visium/dapi_images/',cell2mat(all_vis_id(vi)),'_dapi.jpg'])

end
%  eval(['export_fig ','all_dapi.pdf -r 600']);    
% close all
