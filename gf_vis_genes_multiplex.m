%% run each visium genes
clear all
close all
clc
direct='/data/Technion_analysis/goldfish/visium/multiplex_genes';
cd(direct)


%% load
disp('gf_v')
load('/data/Technion_analysis/goldfish/visium/cropped/091922_three/Sorted.mat','data_orig_all_sorted','sample_sorted','bar_ar_sorted','T_cells_tmp') % vdata
load('/data/Technion_analysis/goldfish/visium/cropped/091922_three/Orignal.mat', 'geneid_all')%vgenes
v_gen=geneid_all;
% v_data=data_orig_sorted_all;
load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_new/non_normalize_all_data_15_8_2022.mat','g_geneid','g_cluster_name','g_data','g_ca') % cgenes
% load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_from_stav_june/gf_mean_data1.mat','mean_data1') % cdata
% c_all_data=g_data;
% c_data=mean_data1;
c_gen=g_geneid;
clusteruni=g_cluster_name;
all_vis_id=natsort(unique(sample_sorted));
all_vis_id(5,:)=[];
% v_gen=geneid_all;
v_data=data_orig_all_sorted;

%% plot genes
[p,~]=numSubplots(length(all_vis_id));
% adjpos = get(hf1, 'Position'); % get new fig position
% save('adjpos.mat','adjpos')
% load('/data/Technion_analysis/goldfish/visium/adjpos.mat')
set(0,'DefaultFigureVisible','on');
gsummary=readtable('/data/Technion_analysis/goldfish/visium/clusters_june/gsummary.csv');
% v_geni=string(table2array(gsummary(:,1)));
 v_geni=upper(["GAD2","SLC17A7";"SLC32A1","SLC17A6"]);% for specific gene list
% poolobj = parpool(8);
%% run on each gene combo
markergene=zeros(size(v_geni,2),length(sample_sorted));

for gnx=1:size(v_geni,2)
    %     search for specific gene
    gnxi=find(v_gen==v_geni(1,gnx)) % for specific gene
    gnxii=find(v_gen==v_geni(2,gnx)) % for specific gene
    %     gnxi=gnx;
    %   `  ha= tight_subplot(p(1),p(2),[.05 .05],[.01 .05],[.01 .01]);
    ha = @(m,n,p) subtightplot (m, n, p,[.01 .01],[.01 .01],[.01 .01]);
    markergene1=data_orig_all_sorted(gnxi,:);
    markergene2=data_orig_all_sorted(gnxii,:);
    markergene(gnx,:)=sum([markergene1;markergene2]);
end
    normmark=markergene./max(markergene);
    nans=logical(sum(isnan(normmark)));
    colmar=[normmark(1,:)', zeros(length(normmark),1) ,zeros(length(normmark),1)];% red
    colmag=[zeros(length(normmark),1),normmark(2,:)', zeros(length(normmark),1) ]; % green
    colmab=[zeros(length(normmark),1), zeros(length(normmark),1),normmark(2,:)']; % blue
    colmap=colmar+colmag;% combo of two 
    colmap(nans,:)=0.8.*ones(sum(nans),3);
%% run each vis alone
    hf1=figure('color','w','units','normalized','outerposition',[0 0 1 1]);

    ax_link=[];

    for vii=1:length(all_vis_id)
        %% find specifc visum
        % set(0,'DefaultFigureVisible','on'); % supress figure;
        % slice_path = fullfile(trial, ['corr_c2v_',char(all_vis_id(vi))]);
        vii
        %         cd(slice_path);
        curr_v=all_vis_id(vii);% example vis name
        v_id=find(string(sample_sorted)==curr_v);
        xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
  
            ax(vii)=ha(3,8,vii);


                        hold on
                        
%                         scatter(xyv(:,1),xyv(:,2),30,[0.8,0.8,0.8],'filled');
                        scatter(xyv(:,1),xyv(:,2), 30, colmap(v_id,:),'filled');

                        set(gca, 'XDir','reverse')
            axis equal;
            axis tight;
            axis off;
%             %         title(all_vis_id(vii))
%             %     linkaxes(ax_link,'xy');
if vii<9
            camroll(180)

  
end
            ax_link=[ax_link,ax(vii)];

%         end
    end % vis type
    % resacale
    sgtitle('INH-R | GLU-G','Interpreter','none')
    namepdf=char(v_gen(gnxi));
%     try
        allYLim = get(ax_link, {'YLim'});
        allYLim = cat(2, allYLim{:});
        set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
        allYLim = get(ax_link, {'XLim'});
        allYLim = cat(2, allYLim{:});
        set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
%     end
colormap
    %% save as png
%     gnx
% %     eval(['export_fig ',[direct,'/',namepdf],'.pdf -r 600']);
    %     set(gcf,'Renderer','OpenGL');
    %     set(gcf,'PaperPositionMode','auto');
    %     save2png([direct,'/',namepdf],gcf,700)

%     close(hf1)

 % gene
% delete(poolobj)