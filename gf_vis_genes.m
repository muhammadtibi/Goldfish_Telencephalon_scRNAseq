%% run each visium genes
clear all
close all
clc
direct='/bigdata/web_pilot/trial';
cd(direct)


%% load
disp('gf_v')
visdirect='/data/Technion_analysis/goldfish/visium/cropped/091922_two';
load([visdirect,'/Sorted.mat'],'data_orig_all_sorted','sample_sorted','bar_ar_sorted','T_cells_tmp') % vdata
load([visdirect,'/Orignal.mat'], 'geneid_all')%vgenes
% v_gen=geneid_all;
% v_data=data_orig_sorted_all;
% load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_from_stav_june/g_matrix.mat','g_geneid','g_cluster_name','g_data','g_ca') % cgenes
% load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_from_stav_june/gf_mean_data1.mat','mean_data1') % cdata
% c_all_data=g_data;
% c_data=mean_data1;
% c_gen=g_geneid;
% clusteruni=g_cluster_name;
all_vis_id=natsort(unique(sample_sorted));
% all_vis_id(5)=[];
% all_vis_id([7 8])=all_vis_id([8 7]);
v_gen=geneid_all;
v_data=data_orig_all_sorted;
zlevel=str2double(extractAfter(all_vis_id,'_'));

%% plot genes
[p,~]=numSubplots(length(all_vis_id));
% adjpos = get(hf1, 'Position'); % get new fig position
% save('adjpos.mat','adjpos')
% load('/data/Technion_analysis/goldfish/visium/adjpos.mat')
set(0,'DefaultFigureVisible','off');
% gsummary=readtable('/data/Technion_analysis/goldfish/visium/clusters_june/gsummary.csv');
% v_geni=string(table2array(gsummary(:,1)));
v_geni=string(geneid_all);%upper(["GAD2","GFAP","SLC17A7","SLC17A6","GAD2","OTP","CCK","PYY","SST","OLIG2","SOX3","LHX8"]);% for specific gene list
% load('gen_list.mat')
% v_geni=upper(string(gen_list));
% poolobj = parpool(8);
for gnx=1:length(v_geni)% run on each scell type
    %%
    %     search for specific gene
    gnxi=find(v_gen==v_geni(gnx)) % for specific gene
    %     gnxi=gnx;
    hf1=figure;
    set(gcf,'color','w','position',[73,15,1100,400]);
    %     hf1=figure('Name',char(v_gen(gnxi)),'NumberTitle','off','color','w','units','normalized','outerposition',[0 0 1 1]);
    %   `  ha= tight_subplot(p(1),p(2),[.05 .05],[.01 .05],[.01 .01]);
    ha = @(m,n,p) subtightplot (m, n, p,[.001 .001],[.001 .001],[.001 .001]);
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
        % scatter(xyv(:,2),xyv(:,1));
        % numValues = length(xyv);
        %%
        markergene=log(data_orig_all_sorted(gnxi,:)+1);
        tmpthlow = prctile(markergene(markergene>0),50);
        tmpthhigh = prctile(markergene(markergene>0),90);
        markergene(markergene>tmpthhigh) = tmpthhigh;
        markergene(markergene<tmpthlow) = tmpthlow;
        c_rgb = [1,0,0];
        %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
        %     ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
        try
            markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
            ax(vii)=ha(2,8,vii);

            hold on
            scatter(xyv(:,1),-xyv(:,2),20, markergene_color(v_id,:),'filled');
            % %             {NO HEXAGON}

            
            axis equal;
            axis tight;
            axis off;
%             sgtitle(v_geni(gnx))
            %         title(all_vis_id(vii))
            ax_link=[ax_link,ax(vii)];
            %     linkaxes(ax_link,'xy');
            camroll(180)

        end
    end % vis type
    % resacale
    %     sgtitle(v_gen(gnxi),'Interpreter','none')
    namepdf=char(v_gen(gnxi));
    try
        allYLim = get(ax_link, {'YLim'});
        allYLim = cat(2, allYLim{:});
        set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
        allYLim = get(ax_link, {'XLim'});
        allYLim = cat(2, allYLim{:});
        set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
    end
    % save as png
    gnx
    %     eval(['export_fig ',[direct,'/',namepdf],'.pdf -r 600']);
    %     set(gcf,'Renderer','OpenGL');
    %     set(gcf,'PaperPositionMode','auto');
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    directx=[direct,'/VIS/'];
    % % %     save2png([directx,namepdf,'_visium'],gcf,300)
    %     close(hf1)

end % gene
% delete(poolobj)
%% Do it with each brain alone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot genes
[p,~]=numSubplots(length(all_vis_id));
% adjpos = get(hf1, 'Position'); % get new fig position
% save('adjpos.mat','adjpos')
% load('/data/Technion_analysis/goldfish/visium/adjpos.mat')
set(0,'DefaultFigureVisible','on');
% gsummary=readtable('/data/Technion_analysis/goldfish/visium/clusters_june/gsummary.csv');
% v_geni=string(table2array(gsummary(:,1)));
v_geni=upper(["CNR1","NR2F2","SST","ETV1","ELAVL4","MEIS2"]);% for specific gene list
% v_geni=string(geneid_all);% uncomment
% load('gen_list.mat')
% v_geni=upper(string(gen_list));
% poolobj = parpool(8);
%%
for gnx=1:length(v_geni)% run on each scell type
    %%
    %     search for specific gene
    gnx
    gnxi=find(v_gen==v_geni(gnx)); % for specific gene
    %     gnxi=gnx;

    % run each goldfish brain
    sec=0;
    for gf=2:3
        hf1=figure;

        set(gcf,'color','w','position',[73           1        1109         192]);
        %     hf1=figure('Name',char(v_gen(gnxi)),'NumberTitle','off','color','w','units','normalized','outerposition',[0 0 1 1]);
        %   `  ha= tight_subplot(p(1),p(2),[.05 .05],[.01 .05],[.01 .01]);
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
            %         cd(slice_path);
            curr_v=all_vis_id(vii);% example vis name
            v_id=find(string(sample_sorted)==curr_v);
            xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
            % scatter(xyv(:,2),xyv(:,1));
            % numValues = length(xyv);
            %%
            markergene=log(data_orig_all_sorted(gnxi,:)+1);
            tmpthlow = prctile(markergene(markergene>0),50);
            tmpthhigh = prctile(markergene(markergene>0),90);
            markergene(markergene>tmpthhigh) = tmpthhigh;
            markergene(markergene<tmpthlow) = tmpthlow;
            c_rgb = [1,0,0];
            %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
            %     ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
            try
                markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
                ax(vc)=ha(1,8,vc);


                %         if vii<3 % flip
                %             scatter(xyv(:,1),-xyv(:,2),20,'k');
                hold on
                %                         scatter(xyv(:,1),xyv(:,2), 20, markergene_color(v_id,:),'filled');
                %                         set(gca, 'XDir','reverse')
                %                     else
                %                         scatter(xyv(:,1),-xyv(:,2),1,'k');
                % %             {NOT HEXAGON}
                hold on
                if sec==0
                    scatter(-xyv(:,1),xyv(:,2),12, markergene_color(v_id,:),'filled');


                else
                    scatter(xyv(:,1),-xyv(:,2),20, markergene_color(v_id,:),'filled');
                end

           
                axis equal;
                axis tight;
                axis off;
                %         title(all_vis_id(vii))
                ax_link=[ax_link,ax(vc)];
                %     linkaxes(ax_link,'xy');
                camroll(180)

            end
        end % vis type
        % resacale
        %     sgtitle(v_gen(gnxi),'Interpreter','none')
        namepdf=char(v_gen(gnxi));
        try
            allYLim = get(ax_link, {'YLim'});
            allYLim = cat(2, allYLim{:});
            set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
            allYLim = get(ax_link, {'XLim'});
            allYLim = cat(2, allYLim{:});
            set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
        end
        % save as png
        gnx
        eval(['export_fig ',[direct,'/',namepdf],'.pdf -r 300']);
        %     set(gcf,'Renderer','OpenGL');
        %     set(gcf,'PaperPositionMode','auto');
        set(gcf,'Renderer','OpenGL');
        set(gcf,'PaperPositionMode','auto');
        directx=[direct,'/VIS',num2str(gf),'/'];
        %             save2png([directx,namepdf,'_visium_',num2str(gf)],gcf,300)
        close(hf1)
        sec=sec+8;
    end % gf
end % gene
% delete(poolobj)
%%
% Do it with hexagons instead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot genes
[p,~]=numSubplots(length(all_vis_id));
% adjpos = get(hf1, 'Position'); % get new fig position
% save('adjpos.mat','adjpos')
% load('/data/Technion_analysis/goldfish/visium/adjpos.mat')
set(0,'DefaultFigureVisible','on');
gsummary=readtable('/data/Technion_analysis/goldfish/visium/clusters_june/gsummary.csv');
% v_geni=string(table2array(gsummary(:,1)));
load('/data/Technion_analysis/goldfish/visium/spatial_fs/gen_list.mat', 'gen_list')
v_geni=upper(gen_list);% for specific gene list
% load('gen_list.mat')
% v_geni=upper(string(gen_list));
% poolobj = parpool(8);
for gnx=1:length(v_geni)% run on each scell type
    %%
    %     search for specific gene
    gnxi=find(v_gen==v_geni(gnx)) % for specific gene
    %     gnxi=gnx;
    hf1=figure;
    set(gcf,'color','w','position',[20,20,1020,120]);
    %     hf1=figure('Name',char(v_gen(gnxi)),'NumberTitle','off','color','w','units','normalized','outerposition',[0 0 1 1]);
    %   `  ha= tight_subplot(p(1),p(2),[.05 .05],[.01 .05],[.01 .01]);
    ha = @(m,n,p) subtightplot (m, n, p,[.01 .01],[.01 .01],[.01 .01]);
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
        % scatter(xyv(:,2),xyv(:,1));
        % numValues = length(xyv);
        %%
        markergene=data_orig_all_sorted(gnxi,:);
        tmpthlow = prctile(markergene(markergene>0),50);
        tmpthhigh = prctile(markergene(markergene>0),90);
        markergene(markergene>tmpthhigh) = tmpthhigh;
        markergene(markergene<tmpthlow) = tmpthlow;
        c_rgb = [1,0,0];
        %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
        %     ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
        try
            markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
                interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
                ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
            %             ax(vii)=ha(2,8,vii);


            %         if vii<3 % flip
            %             scatter(xyv(:,1),-xyv(:,2),20,'k');
            hold on
            %                         scatter(xyv(:,1),xyv(:,2), 20, markergene_color(v_id,:),'filled');
            %                         set(gca, 'XDir','reverse')
            %                     else
            %                         scatter(xyv(:,1),-xyv(:,2),1,'k');
            % %             {NOT HEXAGON}
            hold on
            pgon_all=[];
            for cor=1:length(xyv)
                pgon = nsidedpoly(6,'Center',xyv(cor,:),'Radius',30);
                pgon = rotate(pgon,30)
                pgon_all=[pgon_all,pgon];
            end
            plot(pgon_all)

            %                         scatter(xyv(:,1),-xyv(:,2),2, markergene_color(v_id,:),'filled');
            % %             {NO HEXAGON}

            % {HEXAGON}
            %             cData=mean(markergene_color(v_id),2);%ones(length(xyv),1);
            %             %         cData(v_id)=1;
            %
            %             if sum(cData)==0
            %                 map=[0.5, 0.5,0.5;0.5, 0.5,0.5];
            %
            %             else
            %                 map=[0.5, 0.5,0.5; 1, 0 ,0];
            %
            %             end
            %             mnlim=min(xyv);
            %             mxlim=max(xyv);
            %             rhex=40;
            %             hexScatter(xyv(:,1),xyv(:,2),cData,[mnlim(1)-2*rhex,mxlim(1)+2*rhex],[mnlim(2)-2*rhex,mxlim(2)+2*rhex],rhex,0) ;
            %             colormap(ax(vii),map)
            %             colorbar('off')
            %             set(gca, 'XDir','reverse')
            %             %         end
            %             %            if vii>8
            %             camroll(180)
            %         end
            axis off;
            axis equal;
            % {HEXAGON}
            %         end
            %         end
            axis equal;
            axis tight;
            axis off;
            %         title(all_vis_id(vii))
            ax_link=[ax_link,ax(vii)];
            %     linkaxes(ax_link,'xy');
            camroll(180)

        end
    end % vis type
    % resacale
    %     sgtitle(v_gen(gnxi),'Interpreter','none')
    namepdf=char(v_gen(gnxi));
    try
        allYLim = get(ax_link, {'YLim'});
        allYLim = cat(2, allYLim{:});
        set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
        allYLim = get(ax_link, {'XLim'});
        allYLim = cat(2, allYLim{:});
        set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
    end
    %% save as png
    gnx
    %     eval(['export_fig ',[direct,'/',namepdf],'.pdf -r 600']);
    %     set(gcf,'Renderer','OpenGL');
    %     set(gcf,'PaperPositionMode','auto');
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    directx=[direct,'/VIS/'];
    save2png([directx,namepdf,'_visium'],gcf,700)
    close(hf1)

end % gene
%% extra  visium 3d genes
ilv=9:16;%interleave(1:8,9:16);
ilv_vis_id=all_vis_id(ilv);
ilv_z=zlevel(ilv);
yall=[];
xall=[];
zall=[];
cmapall=[];
mvisxall=[];
mvisyall=[];
v_geni=["NEUROD6"];% for specific gene list
%     search for specific gene
gnxi=find(v_gen==v_geni) % for specific gene
markergene=data_orig_all_sorted(gnxi,:);
tmpthlow = prctile(markergene(markergene>0),50);
tmpthhigh = prctile(markergene(markergene>0),90);
markergene(markergene>tmpthhigh) = tmpthhigh;
markergene(markergene<tmpthlow) = tmpthlow;
c_rgb = [1,0,0];

markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
for vii=9:length(ilv_vis_id)+8
    vii
    curr_v=ilv_vis_id(vii-8);% example vis name
    v_id=find(string(sample_sorted)==curr_v);
    xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
    axis off;
    axis equal;
    hold on

    cmap=markergene_color(v_id,:);
    

    zcurr=ilv_z(vii-8).*ones(length(xyv),1);
    mvisx=mean(xyv(:,1));
    mvisy=mean(xyv(:,2));
    if vii==9
        refx=mvisx;
        refy=mvisy;
        xcurr=xyv(:,1);
        ycurr=xyv(:,2);
        scatter3(zcurr,xcurr,ycurr,100,cmap,'filled');
        %         scatter3(ilv_z(vii-8),refx,refy,100,'k','filled'); % centorid

    else
        deltax=mvisx-refx;
        deltay=mvisy-refy;
        xcurr=xyv(:,1)-deltax;
        ycurr=xyv(:,2)-deltay;
        scatter3(zcurr,xcurr,ycurr,100,cmap,'filled');
        %         scatter3(ilv_z(vii-8),mvisx-deltax,mvisy-deltay,100,'k','filled');         % centroid
    end


    xall=[xall;xcurr];
    yall=[yall;ycurr];
    cmapall=[cmapall;cmap];
    zall=[zall;zcurr];
    mvisxall=[mvisxall,mvisx];
    mvisyall=[mvisyall,mvisy];
    clear cmap

end % clusters
view(30,45)


%% add polyshape
% figure('color','w');

shp = alphaShape(zall,xall,yall);
plot(shp,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.5);
%% rotate
for i=1:100
    camorbit(10,0,'camera')
    drawnow
end

