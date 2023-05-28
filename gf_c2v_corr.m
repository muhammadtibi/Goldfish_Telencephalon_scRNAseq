%% load for my own "rctd" (only rho corr.)
close all
clear all
clc
trial='/data/Technion_analysis/goldfish/visium';
set(0,'DefaultFigureVisible','off');% off / on
set(0,'DefaultFigureWindowStyle','normal')
direct='/data/Technion_analysis/goldfish/visium/cropped/091922_two';
cd(direct)
corrx=['/data/Technion_analysis/goldf' ...
    'ish/scRNAseq_gf/correlation/'];
% cd(corrx)
%% get col
load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_new/non_normalize_all_data_15_8_2022.mat', 'g_cluster_name')
nms=g_cluster_name;

ugf=unique(nms,'stable');
n_gi=contains(ugf,'g_');
n_g= regexprep(ugf(n_gi),'_','-');
n_g= regexprep(n_g,'g-','');

r1=n_g(contains(n_g,'GABA'));

% make colormap
gabacol=winter(length(r1));


r2=n_g(contains(n_g,'GLUT-'));
glutcol=spring(length(r2));

nncol=summer(length(nms)-length(r2)-length(r1));



colmap=[glutcol;gabacol;nncol];
%% visium
disp('load data')
load([direct,'/Sorted.mat'],'data_orig_all_sorted','sample_sorted','bar_ar_sorted', 'T_cells_tmp') % vdata
load([direct,'/Orignal.mat'], 'geneid_all')%vgenes
v_gen=geneid_all;
v_data=data_orig_all_sorted;
load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_new/non_normalize_all_data_15_8_2022.mat','g_geneid','g_cluster_name','non_normalize_data','g_ca') % cgenes
% load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_from_stav_june/gf_mean_data1.mat','mean_data1') % go to gf_violin %% calculate mean_data1
c_all_data=non_normalize_data;
c_gen=g_geneid;
clusteruni=g_cluster_name;
all_vis_id=natsort(unique(sample_sorted));
zlevel=str2double(extractAfter(all_vis_id,'_'));
% all_vis_id([5])=[];
% gfx=all_vis_id(1:8);
% all_vis_id(1:8)=[];
% all_vis_id=[all_vis_id;gfx];
%% caculate mean
disp('calc. mean')
mean_data1=zeros(length(c_gen),length(clusteruni));

for cc=1:length(clusteruni)
    cc
    gxx = find(strcmpi(clusteruni(cc),g_ca));
    mean_data1(:,cc)=mean(c_all_data(:,gxx),2);

end
c_data=mean_data1;

%% intersect genes and corr. scell to sspot
disp('intersect')
[inBoth,ic,iv]=intersect(c_gen,v_gen);% inBoth is the genes now
m_data_x=c_data(ic,:);
c_all_data_x=c_all_data(ic,:);
v_data_x=v_data(iv,:);

%% if : feature selection
% filter vis
% disp('feature selection')
% in_v = find(sum(v_data_x>0,2)>5 & sum(v_data_x>0,2)<length(v_data_x(1,:))*0.5);
% v_data_x=v_data_x(in_v,:);
% inBothv=inBoth(in_v);
% in_c = find(sum(c_all_data_x>0,2)>5 & sum(c_all_data_x>0,2)<length(c_all_data_x(1,:))*0.5);
% m_data_x=m_data_x(in_c,:);
% c_all_data_x=c_all_data_x(in_c,:);
% inBothc=inBoth(in_c);
% nix = contains(g_ca,'GABA','IgnoreCase',true);
% nex = contains(g_ca,'Glut','IgnoreCase',true);
% nnx = ~nix & ~nex;
% %
% % corr_filtni = cv_vs_m_selection(c_all_data_x(:,nix==1),inBothc,[],1,0);
% % corr_filtne = cv_vs_m_selection(c_all_data_x(:,nex==1),inBothc,[],1,0);
% % corr_filtnn = cv_vs_m_selection(c_all_data_x(:,nnx==1),inBothc,[],1,0);
% corr_filt = cv_vs_m_selection(c_all_data_x,inBothc,[2000],1,0);
% %
% % corr_filtx=[corr_filtni;corr_filtne;corr_filtnn];
% % corr_filt=unique(corr_filtx);
%
% %
%
% [inBothx,ic,iv]=intersect(inBothc(corr_filt),inBothv);% inBoth is the genes now
% m_data_y=m_data_x(corr_filt(ic),:);
% v_data_y=v_data_x(iv,:);
%%  or else : new-c orrelation based on markertable featured genes
% RHO_ALL=zeros(length(clusteruni),size(v_data_x,2));
% for typx=1:length(clusteruni)
%     typx
%     [RHO,~] = corr(cent_norm(m_data_x(:,typx)),cent_norm(v_data_x),'Type','Spearman');
%     RHO_ALL(typx,:)=RHO;
% end
% [RHO_ALL,~] = corr(cent_norm(m_data_y),cent_norm(v_data_y),'Type','Pearson');
disp('markertable- fs)')
RHO_ALL=zeros(length(clusteruni),size(v_data,2));
% cd('/data/Technion_analysis/goldfish/scRNAseq_gf/correlation')
[~,~,numn]=unique(g_ca,'stable');
[ind_gr_tmp_mark,vi0] = markertablefeatures_tmp(numn,c_all_data_x,10,1);% top_genes by me
%% then smooth : by or or or
%% or geomean
% for cc=1:length(clusteruni)
%     cc
%     marker_tmp=ind_gr_tmp_mark(:,cc);
%     expression=(v_data_x(marker_tmp,:));
%     w=vi0(:,cc);
%     w=w/sum(w);
%     ww=expression.*(repmat(w,1,size(expression,2)));
%     RHO_ALL(cc,:)=geomean(ww);%sum(v_data_x(marker_tmp,:).*repmat(w,1,size(v_data,2)));
% 
% end
% %% or impute vis
% prj=pca_wis(log2(v_data_y'+1),20);
% knbr=knn_mat(prj,10,1);
% v_data_imp=zeros(size(v_data));
% for k=1:size(v_data,2)
%     k
%     v_data_imp(:,k)=mean(v_data(:,logical(knbr(:,k))),2);
% end
% v_data_x=v_data_imp;
%% or sigmoid
sgmd=@(x) (1./(1+exp(-20*(x-0.6))));
for cc=1:length(clusteruni)
    cc
    marker_tmp=ind_gr_tmp_mark(:,cc);
    expression=(v_data_x(marker_tmp,:)>0);
    w=vi0(:,cc);
    w=w/sum(w);
    ww=expression.*(repmat(w,1,size(expression,2)));
    RHO_ALL(cc,:)=sgmd(sum(ww));%sum(v_data_x(marker_tmp,:).*repmat(w,1,size(v_data,2)));

end

% RHO_ALL=pdist2(m_data_y'>0,v_data_y'>0,'jaccard');
%     [RHO_ALL,~] = corr(m_data_x,v_data_x,'Type','Pearson');

writecell([cellstr(clusteruni),num2cell(RHO_ALL)],[direct,'/RHO_ALL.csv']);

%% create folders for all vises
% for vi= 1:length(all_vis_id)
% slice_path = fullfile(trial, ['corr_c2v_',char(all_vis_id(vi))]);
%     if ~exist(slice_path)
%         mkdir(slice_path)
%     end
% end
%% load rhos and normalize and create save dir
% RHO_norm=zeros(size(RHO_ALL));
% for typx=1:length(clusteruni)
%     RHO_N = RHO_ALL(typx,:)./max(RHO_ALL(typx,:));% per cell type
%     RHO_norm(typx,:)=RHO_N;
% end
% slice_path = fullfile(trial, 'corr_c2v_celltype');
% if ~exist(slice_path)
%     mkdir(slice_path)
% end
% RHO_norm=RHO_ALL./sum(RHO_ALL);
% RHO_norm=cent_norm(RHO_ALL);
%% plot corr.visiums per cell type
% [p,~]=numSubplots(length(all_vis_id));
% adjpos = get-(f, 'Position'); % get new fig position
% load('/data/Technion_analysis/Amygdala/Vis_dimorphisim/adjpos.mat')
cd(corrx);
set(0,'DefaultFigureVisible','on');% off / on

% load('/data/Technion_analysis/goldfish/scRNAseq_gf/sc_vis_pca/RHO_ALL','RHO_ALL');% saved above
% RHO_ALL(isnan(RHO_ALL))=0;
for typx=35%:length(clusteruni)% run on each scell type
    %%

    hf1=figure('Name',char(clusteruni(typx)),'NumberTitle','off','color','w','units','normalized','outerposition',[0 0 1 1]);
    ha = @(m,n,p) subtightplot (m, n, p,[.01 .01],[.01 .01],[.01 .01]);
    ax_link=[];

    for vii=1:length(all_vis_id)
        %% find specifc visum
        % set(0,'DefaultFigureVisible','on'); % supress figure;
        % slice_path = fullfile(trial, ['corr_c2v_',char(all_vis_id(vi))]);
        vii
        curr_v=all_vis_id(vii);% example vis name
        v_id=find(string(sample_sorted)==curr_v);
        xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
        % scatter(xyv(:,2),xyv(:,1));
        % numValues = length(xyv);
        %%
        markergene=RHO_ALL(typx,:);
        %         tmpthlow = prctile(markergene(markergene>0),50);
        %         tmpthhigh = prctile(markergene(markergene>0),98);
        %         markergene(markergene>tmpthhigh) = tmpthhigh;
        %         markergene(markergene<tmpthlow) = tmpthlow;
        %         tmpthhigh(isnan(tmpthhigh))=0;
        %         tmpthlow(isnan(tmpthlow))=0;
        c_rgb = [0,100,0]/255;
        %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
        %     ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
        try
            markergene_color = [interp1([min(markergene),max(markergene)],[0.8,c_rgb(1)],markergene'),...
                interp1([min(markergene),max(markergene)],[0.8,c_rgb(2)],markergene')...
                ,interp1([min(markergene),max(markergene)],[0.8,c_rgb(3)],markergene')];
        catch
            markergene_color=zeros(length(markergene),3);
        end
        %         ax(vii)=ha(1,length(all_vis_id),vii);
        ax(vii)=ha(2,8,vii);
        if vii<17
            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2), 40, markergene_color(v_id,:),'filled');
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2),45, markergene_color(v_id,:),'filled');
            camroll(180)

        else

            hold on
            scatter(xyv(:,1),-xyv(:,2), 30, markergene_color(v_id,:),'filled');
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2),30, markergene_color(v_id,:),'filled');

        end

        %         if vii>8
        %         end
        axis equal;
        axis tight;
        axis off;
        %         title(all_vis_id(vii))
        ax_link=[ax_link,ax(vii)];
        %     linkaxes(ax_link,'xy');
        %     camroll(90)


    end % vis type
    % resacale
    %     sgtitle(clusteruni(typx),'Interpreter','none')
    namepdf=char(clusteruni(typx));
    allYLim = get(ax_link, {'YLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
    allYLim = get(ax_link, {'XLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
    typx
    %     sgtitle(inBoth(ind_gr_tmp_mark(:,typx)))
    %%
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    save2png([corrx,namepdf],gcf,300)

    %     eval(['export_fig ',namepdf,'.pdf -r 300']);

    close(hf1)

end % cell type

%% plot corr.visiums per cell type
% [p,~]=numSubplots(length(all_vis_id));
% adjpos = get-(f, 'Position'); % get new fig position
% load('/data/Technion_analysis/Amygdala/Vis_dimorphisim/adjpos.mat')
cd(corrx);
set(0,'DefaultFigureVisible','on');% off / on

% load('/data/Technion_analysis/goldfish/scRNAseq_gf/sc_vis_pca/RHO_ALL','RHO_ALL');% saved above
% RHO_ALL(isnan(RHO_ALL))=0;
for typx=1:length(clusteruni)% run on each scell type
                 sec=0;

        for gf=1:2
                            hf1=figure;

        set(gcf,'color','w','position',[73           1        1109         192]);
        ha = @(m,n,p) subtightplot (m, n, p,[.001 .001],[.001 .001],[.001 .001]);
        ax_link=[];
        secend=sec+8;
        vc=0;
        for vii=(1+sec):secend
                        vc=vc+1;
        %% find specifc visum
        % set(0,'DefaultFigureVisible','on'); % supress figure;
        % slice_path = fullfile(trial, ['corr_c2v_',char(all_vis_id(vi))]);
        vii
        curr_v=all_vis_id(vii);% example vis name
        v_id=find(string(sample_sorted)==curr_v);
        xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
        % scatter(xyv(:,2),xyv(:,1));
        % numValues = length(xyv);
        %%
        markergene=RHO_ALL(typx,:);
        %         tmpthlow = prctile(markergene(markergene>0),50);
        %         tmpthhigh = prctile(markergene(markergene>0),98);
        %         markergene(markergene>tmpthhigh) = tmpthhigh;
        %         markergene(markergene<tmpthlow) = tmpthlow;
        %         tmpthhigh(isnan(tmpthhigh))=0;
        %         tmpthlow(isnan(tmpthlow))=0;
        c_rgb = [0 100 0]/255;
        %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
        %     ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
        try
            markergene_color = [interp1([min(markergene),max(markergene)],[0.8,c_rgb(1)],markergene'),...
                interp1([min(markergene),max(markergene)],[0.8,c_rgb(2)],markergene')...
                ,interp1([min(markergene),max(markergene)],[0.8,c_rgb(3)],markergene')];
        catch
            markergene_color=zeros(length(markergene),3);
        end
        %         ax(vii)=ha(1,length(all_vis_id),vii);
        ax(vc)=ha(1,8,vc);
        if vii<17
            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2), 20, markergene_color(v_id,:),'filled');
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2),20, markergene_color(v_id,:),'filled');
            camroll(180)

        else

            hold on
            scatter(xyv(:,1),-xyv(:,2), 30, markergene_color(v_id,:),'filled');
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2),30, markergene_color(v_id,:),'filled');

        end

        %         if vii>8
        %         end
        axis equal;
        axis tight;
        axis off;
        %         title(all_vis_id(vii))
        ax_link=[ax_link,ax(vc)];
        %     linkaxes(ax_link,'xy');
        %     camroll(90)


        end  % vis type
    % resacale
    %     sgtitle(clusteruni(typx),'Interpreter','none')
    namepdf=[num2str(gf),'_',char(clusteruni(typx))];
    allYLim = get(ax_link, {'YLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
    allYLim = get(ax_link, {'XLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
    typx
    %     sgtitle(inBoth(ind_gr_tmp_mark(:,typx)))
    %%
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    save2png([corrx,namepdf],gcf,300)

    %     eval(['export_fig ',namepdf,'.pdf -r 300']);

     close(hf1)
     sec=sec+8;
        end
end % cell type
%% BATCH correlations (for hannah "twist" region)

cd(corrx);
set(0,'DefaultFigureVisible','on');% off / on


groupx= [20, 21, 30, 31];
groupy= [5, 12, 16, 18, 19];
markergenex=mean(RHO_ALL(groupx,:));% mean groups
markergeney=mean(RHO_ALL(groupy,:));% mean groups

% for typx=1:2 % run on each BATCH scell type

    hf1=figure('color','w','units','normalized','outerposition',[0 0 1 1]);
    ha = @(m,n,p) subtightplot (m, n, p,[.01 .01],[.01 .01],[.01 .01]);
    ax_link=[];

    for vii=1:length(all_vis_id)
        %% find specifc visum
      
        vii
        curr_v=all_vis_id(vii);% example vis name
        v_id=find(string(sample_sorted)==curr_v);
        xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
        % scatter(xyv(:,2),xyv(:,1));
        % numValues = length(xyv);
        %%
%         markergene=mean(RHO_ALL(groupx,:));% mean groups

        c_rgbx = [0 200 0]/255;
        c_rgby = [200 0 0]/255;
        try
            markergene_colorx = [interp1([min(markergenex),max(markergenex)],[0.8,c_rgbx(1)],markergenex'),...
                interp1([min(markergenex),max(markergenex)],[0.8,c_rgbx(2)],markergenex')...
                ,interp1([min(markergenex),max(markergenex)],[0.8,c_rgbx(3)],markergenex')];
                  
            markergene_colory = [interp1([min(markergeney),max(markergeney)],[0.8,c_rgby(1)],markergeney'),...
                interp1([min(markergeney),max(markergeney)],[0.8,c_rgby(2)],markergeney')...
                ,interp1([min(markergeney),max(markergeney)],[0.8,c_rgby(3)],markergeney')];
            % combine and cauclate mean color 
            markergene_color= (markergene_colorx+ markergene_colory)/2;
        catch
            markergene_color=zeros(length(markergene),3);
        end
        %         ax(vii)=ha(1,length(all_vis_id),vii);
        ax(vii)=ha(2,8,vii);
            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2), 40, markergene_color(v_id,:),'filled');
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2),45, markergene_color(v_id,:),'filled');
            camroll(180)

      

        %         if vii>8
        %         end
        axis equal;
        axis tight;
        axis off;
        %         title(all_vis_id(vii))
        ax_link=[ax_link,ax(vii)];
        %     linkaxes(ax_link,'xy');
        %     camroll(90)


    end % vis type
    % resacale
    %     sgtitle(clusteruni(typx),'Interpreter','none')
%     namepdf=char(clusteruni(typx));
    allYLim = get(ax_link, {'YLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
    allYLim = get(ax_link, {'XLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
    typx
    %     sgtitle(inBoth(ind_gr_tmp_mark(:,typx)))
    %%
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    save2png([corrx,namepdf],gcf,300)

    %     eval(['export_fig ',namepdf,'.pdf -r 300']);

%     close(hf1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% region cluster 2 cell type
all_rois=[];
for typx=1: length(clusteruni)% run on each scell type
    G_all=[];
    for vz=1:length(all_vis_id)  % go over multidirectiores
        clear T
        vz
        curr_v=all_vis_id(vz);% example vis name
        v_id=find(string(sample_sorted)==curr_v);
        markergene=RHO_ALL(typx,v_id);
        xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
        T = array2table(xyv,'VariableNames',{'cent_x','cent_y'});
        %         vdir= ['/data/runs/samples/',char(all_vis_id(vz)),'/'];
        %         dot_roi=zeros(length(xyv),1);
        T.ROI=T_cells_tmp(v_id);
        T.rval=markergene';% get the avlues of RGO for each region
        G=groupsummary(T,"ROI","sum");
        G.vis = repmat(string(curr_v), 1, height(G))';
        G(:,[3 4]) = [];
        %         writetable(G,[char(curr_v),'_sum_mono','.csv']); % save raw excel table with centroids
        G_all=[G_all;G];
    end % region -all or spec

    F=G_all;
    F(:,[4]) = [];
    Fx=groupsummary(F,"ROI","sum");
    R=[Fx.ROI,Fx.sum_sum_rval./Fx.sum_GroupCount];% roi-id and ratio of sum rho value per region deiveid by sum of spots in this region(i.e. =mean)
    all_rois=[all_rois;Fx.ROI];
    %     all_rois=unique(all_rois);
    clear Fx F G_all
    typx
    save([char(clusteruni(typx)),'_ratio.mat'],'R')
end % cell type
%% then loop each R to conactante the matrices with unique roi
u_roi=unique(all_rois,'stable');
all_R= zeros(length(u_roi),length(clusteruni));
for z=1:length(u_roi)
    z
    for typx= 1:length(clusteruni)
        load([char(clusteruni(typx)),'_ratio.mat'],'R')
        rind=find(R(:,1)==u_roi(z));
        all_R(z,typx)=R(rind,2);
    end
end
% tabulate
cell_name=[{'ROI\C_Type'};clusteruni]';
roi_name=u_roi;
all_x=[roi_name,all_R];
all_x=[cell_name;all_x];

Tx=table(all_x);
writetable(Tx,'ratio_ct2roi.csv')
%% chose specific rows of matrix
% alltab=readtable('/data/Technion_analysis/Amygdala/Vis_dimorphisim/c2v_corr/ratio_ct2roi_all.csv'); % read allx/Tx
% all_fullnames= string(st.safe_name);
% all_acroid=zeros(length(sublist),1);
% for xi=1:length(sublist)
%     acroid=find(all_fullnames==sublist(xi));
%     all_acroid(xi)=acroid;
% end
% all_subid=[];
% for yi=1:length(sublist)
%     subid=find(string(table2array(alltab(:,1)))==st.acronym(all_acroid(yi)));
%     all_subid=[all_subid;subid];
% end
% subtab=alltab(all_subid,:);
% numsubtab=table2array(subtab(:,2:end));

%% heatmap
clusteruni=regexprep(clusteruni,'_','-');
clusteruni=regexprep(clusteruni,'g-','');
regstr=strings(1,size(all_R,1));
for rs=1:size(all_R,1)
    regstr(rs)=['region-',num2str(rs)];
end
numsubtab=all_R;
Z = linkage(numsubtab,'ward','correlation');
leaforder = optimalleaforder(Z,pdist(numsubtab));
xz=numsubtab(:);% makes matrix linear for caclculating perc%
figure;
set(gcf,'Color','w')
% subplot(1,2,1)
hden = dendrogram(Z,length(leaforder),'Reorder',leaforder);
set(gca, 'XTickLabel',regstr)
%     axis off

% subplot(1,2,2)
%%
figure;
set(gcf,'Color','w')
cmap=redblue(256);
colormap(cmap)
imagesc(numsubtab(leaforder,:),[prctile(xz,2),prctile(xz,95)])
numrois=1:numel(u_roi);
numtypes=1:numel(clusteruni);
set(gca, 'XTick',numtypes, 'XTickLabel',clusteruni)
set(gca, 'YTick',numrois, 'YTickLabel',regstr)
xtickangle(90)
namepdf='hm_c2roi';
eval(['export_fig ',namepdf,'.pdf -r 600']);
% save2png(['/data/Technion_analysis/goldfish/visium/spatialcluster2sc/',namepdf],gcf,700)
% svae reordered table
all_x=[string(u_roi((leaforder),1)),numsubtab(leaforder,:)];
all_x=[cell_name;all_x];
Ty=table(all_x);
writetable(Ty,'ratio_ct2roi.csv')
save('leaforder.mat','leaforder')



%     set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
% writetable(Ty,'ratio_ct2roi.csv')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot all_cor colormap
hf1=figure('Name','all_clusters','NumberTitle','off','color','w');

ha = @(m,n,p) subtightplot (m, n, p,[.01 .01],[.01 .01],[.01 .01]);
ax_link=[];

for vii=1:length(all_vis_id)
    vii

    ax(vii)=ha(3,8,vii);

    axis off;
    axis equal;
    hold on


    curr_v=all_vis_id(vii);% example vis name
    v_id=find(string(sample_sorted)==curr_v);
    %         v_binary=string(sample_sorted)==curr_v;
    xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
    rho_x=RHO_ALL(:,v_id);
    rho_x(rho_x<0)=0;% zeroing negative values
    %     rho_x=zeros(106,length(v_id));
    %     rho_x(:,1)=0.6*(rho_x(:,1)+1);
    %     [rv,ri]=max(rho_x);
    %     cmap=colmap(ri,:);
    clear cmap
    for x= 1:length(v_id)
        rho_y=rho_x(:,x)/sum(rho_x(:,x));% norm by sum
        rxc=rho_y .* colmap;% rho X colormap weigtheted color
        cmap(x,:)=sum(rxc,1);
    end
    scatter(xyv(:,1),-xyv(:,2),15,cmap,'filled');
    %                 drawnow
    if vii<17
        camroll(180)
    end

    ax_link=[ax_link,ax(vii)];

end % clusters

allYLim = get(ax_link, {'YLim'});
allYLim = cat(2, allYLim{:});
set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
allYLim = get(ax_link, {'XLim'});
allYLim = cat(2, allYLim{:});
set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
% eval(['export_fig ','all_clusters','.pdf -r 600']);
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D visium

figure('color','w');

% ax_link=[];
ilv=9:16;%interleave(1:8,9:16);
ilv_vis_id=all_vis_id(ilv);
ilv_z=zlevel(ilv);
yall=[];
xall=[];
zall=[];
cmapall=[];
mvisxall=[];
mvisyall=[];
for vii=length(ilv_vis_id):-1:1
    vii

    %     ax(vii)=ha(3,8,vii);

    axis off;
    axis equal;
    hold on


    curr_v=ilv_vis_id(vii);% example vis name
    v_id=find(string(sample_sorted)==curr_v);
    %         v_binary=string(sample_sorted)==curr_v;
    xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
    rho_x=RHO_ALL(:,v_id);
    rho_x(rho_x<0)=0;% zeroing negative values

    clear cmap
    for x= 1:length(v_id)
        rho_y=rho_x(:,x)/sum(rho_x(:,x));% norm by sum
        rxc=rho_y .* colmap;% rho X colormap weigtheted color
        cmap(x,:)=sum(rxc,1);
    end
    zcurr=ilv_z(vii).*ones(length(xyv),1);
    mvisx=mean(xyv(:,1));
    mvisy=mean(xyv(:,2));
    if vii==8
        refx=mvisx;
        refy=mvisy;
          xcurr=xyv(:,1);
        ycurr=xyv(:,2);
        scatter3(zcurr,xcurr,ycurr,100,cmap,'filled');
        scatter3(ilv_z(vii),refx,refy,100,'k','filled'); % centorid
      
    else
        deltax=mvisx-refx;
        deltay=mvisy-refy;
        xcurr=xyv(:,1)-deltax;
        ycurr=xyv(:,2)-deltay;
        scatter3(zcurr,xcurr,ycurr,100,cmap,'filled');
        scatter3(ilv_z(vii),mvisx-deltax,mvisy-deltay,100,'k','filled');         % centroid
    end


    xall=[xall;xcurr];
    yall=[yall;ycurr];
    cmapall=[cmapall;cmap];
    zall=[zall;zcurr];
    mvisxall=[mvisxall,mvisx];
    mvisyall=[mvisyall,mvisy];
end % clusters
view(30,45)
%%
shp = alphaShape(zall,xall,yall);
plot(shp,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.5);
for i=1:100
   camorbit(10,0,'camera')
   drawnow
end
%% dendro -v2c
Zpca = linkage(RHO_ALL,'ward','correlation');
Dpca = pdist(RHO_ALL,'correlation');
leaforder_pca = optimalleaforder(Zpca,Dpca);
disp("Dendogram tree")
figure('Name','dendogram','NumberTitle','off');
set(gcf,'color','w')
axes('position',[0.03,0.03,0.25,0.93])
% hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
[hden,T,outperm]  = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','left');
axis off
set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
axes('position',[0.35,0.03,0.63,0.93])
x=squareform(Dpca); imagesc(x(leaforder_pca,leaforder_pca));
colormap('summer')

set(gca,'ytick',[1:length(leaforder_pca)],'yticklabel',n_g(leaforder_pca),'xtick',[],'fontsize',5,'ydir','normal')
%%  part 1: connect dendros as in prj
% load('/data/Technion_analysis/goldfish/visium/outlines/GLUT_prj_sorted.mat','prj_sorted')
rgaba=RHO_ALL(contains(n_g,'GABA'),:);
tress=[[1 6];[7 10];[11 19];[20 23];[24 31];[32 40]];

brho_all=zeros(length(tress),size(RHO_ALL,2));
for i=1:length(tress);
    i
    branchx=tress(i,1):tress(i,2)
    brho_all(i,:)=mean(rgaba(branchx,:));

end
%% part 2: calculate all_zeroid

cd(corrx);
set(0,'DefaultFigureVisible','on');% off / on
idzero_all=zeros(size(brho_all));
colory=distinguishable_colors(length(tress));

for typx=1:length(tress)% run on each scell type
    %%

    hf1=figure('color','w','units','normalized','outerposition',[0 0 1 1]);
    ha = @(m,n,p) subtightplot (m, n, p,[.01 .01],[.01 .01],[.01 .01]);
    ax_link=[];

    for vii=1:length(all_vis_id)

        vii
        curr_v=all_vis_id(vii);% example vis name
        v_id=find(string(sample_sorted)==curr_v);
        xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));

        markergene=brho_all(typx,:);

        c_rgb =colory(typx,:);

        try
            markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
        catch
            markergene_color=zeros(length(markergene),3);
        end
        idzero_all(typx,:)=transpose(sum(markergene_color==0.7,2)==3);
        %         ax(vii)=ha(1,length(all_vis_id),vii);
        ax(vii)=ha(2,8,vii);
        if vii<17
            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            %             scatter(xyv(:,1),-xyv(:,2), 40, markergene_color(v_id,:),'filled');
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2),45, markergene_color(v_id,:),'filled');
            camroll(180)

        else

            hold on
            scatter(xyv(:,1),-xyv(:,2), 30, markergene_color(v_id,:),'filled');
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2),30, markergene_color(v_id,:),'filled');

        end

        %         if vii>8
        %         end
        axis equal;
        axis tight;
        axis off;
        %         title(all_vis_id(vii))
        ax_link=[ax_link,ax(vii)];
        %     linkaxes(ax_link,'xy');
        %     camroll(90)


    end % vis type
    % resacale
    %     sgtitle(clusteruni(typx),'Interpreter','none')
    %     namepdf=char(clusteruni(typx));
    allYLim = get(ax_link, {'YLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
    allYLim = get(ax_link, {'XLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
    typx
    %     sgtitle(inBoth(ind_gr_tmp_mark(:,typx)))
    %%
    %     set(gcf,'Renderer','OpenGL');
    %     set(gcf,'PaperPositionMode','auto');
    save2png([corrx,num2str(typx)],gcf,300)

    %     eval(['export_fig ',namepdf,'.pdf -r 300']);

    %     close(hf1)

end % cell type
%% part 3: plot the multiplexing
cd(corrx);
set(0,'DefaultFigureVisible','on');% off / on
% idzero_all=zeros(size(brho_all));


hf1=figure('color','w','units','normalized','outerposition',[0 0 1 1]);
ha = @(m,n,p) subtightplot (m, n, p,[.01 .01],[.01 .01],[.01 .01]);
ax_link=[];
colory=distinguishable_colors(length(tress));
idzero_all=brho_all>0.05;
% idzero_all(idzero_all==0,:)=0.7;
for vii=1:length(all_vis_id)
    vii
    curr_v=all_vis_id(vii);% example vis name
    v_id=find(string(sample_sorted)==curr_v);
    xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));

    for typx=1:length(tress)% run on each scell type




        c_rgb = idzero_all(typx,v_id)'*colory(typx,:);
        i_rgb=sum(c_rgb')>0;

        %         ax(vii)=ha(1,length(all_vis_id),vii);
        ax(vii)=ha(2,8,vii);
        if vii<17
            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            if typx==1
                scatter(xyv(:,1),-xyv(:,2), 40, [0.7 0.7 0.7],'filled');
            end
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(i_rgb,1),-xyv(i_rgb,2),45, c_rgb(i_rgb,:),'filled');

        else

            hold on
            %             scatter(xyv(:,1),-xyv(:,2), 30, markergene_color(v_id,:),'filled');
            %             set(gca, 'XDir','reverse')

            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            scatter(xyv(:,1),-xyv(:,2),30, markergene_color(v_id,:),'filled');

        end




    end % vis type


    camroll(180)

    axis equal;
    axis tight;
    axis off;
    ax_link=[ax_link,ax(vii)];



    % resacale
    %     sgtitle(clusteruni(typx),'Interpreter','none')
    %     namepdf=char(clusteruni(typx));
    allYLim = get(ax_link, {'YLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
    allYLim = get(ax_link, {'XLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
    typx
    %     sgtitle(inBoth(ind_gr_tmp_mark(:,typx)))
    %%
    %     set(gcf,'Renderer','OpenGL');
    %     set(gcf,'PaperPositionMode','auto');
    %         save2png([corrx,num2str(typx)],gcf,300)

    %     eval(['export_fig ',namepdf,'.pdf -r 300']);

    %     close(hf1)

end  % cell type





