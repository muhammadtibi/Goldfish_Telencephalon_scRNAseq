%% run each visium genes: find axial patterned genes
clear all
close all
clc
direct='/data/Technion_analysis/goldfish/visium/spatial_fs';
cd(direct)
%% update names and colors
load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_new/non_normalize_all_data_15_8_2022.mat','g_geneid','g_cluster_name','g_data','g_ca') % cgenes

geneid=g_geneid;
bx_f=ones(1,length(g_ca))';
data =g_data;%round(g_data./repmat(sum(g_data),length(g_data(:,1)),1)*nrmz);% normalize or no-nrmz = g_data
c=regexprep(g_ca,'_','-');
c=regexprep(c,'g-','');
% c=regexprep(c,'??','immune12');
% c=regexprep(c,'immune12immune12','immune12');
c=regexprep(c,'GLUT-','');
c=regexprep(c,'GABA-','');

cuni=regexprep(g_cluster_name,'_','-');
cuni=regexprep(cuni,'g-','');
% cuni=regexprep(cuni,'??','immune12');
% cuni=regexprep(cuni,'immune12immune12','immune12');
ri1=contains(cuni,'GABA-');
r1=cuni(ri1==1);
r1=regexprep(r1,'GABA-','');

% make colormap
gabacol=winter(length(r1));


ri2=contains(cuni,'GLUT-');
r2=cuni(ri2==1);
r2=regexprep(r2,'GLUT-','');
glut1col=spring(length(r2));


%% load
disp('gf_v')
load('/data/Technion_analysis/goldfish/visium/cropped/091922_two/Sorted.mat','data_orig_all_sorted','sample_sorted','bar_ar_sorted','T_cells_tmp') % vdata
load('/data/Technion_analysis/goldfish/visium/cropped/091922_two/Orignal.mat', 'geneid_all')%vgenes
v_gen=geneid_all;
% v_data=data_orig_sorted_all;
% load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_new/non_normalize_all_data_15_8_2022.mat','g_geneid','g_cluster_name','g_data','g_ca') % cgenes
% load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_new/gf_mean_data1.mat','mean_data1') % cdata
% c_all_data=g_data;
% c_data=mean_data1;
% c_gen=g_geneid;
% clusteruni=g_cluster_name;
all_vis_id=natsort(unique(sample_sorted));
% all_vis_id([1 2],:)=[];
% all_vis_id(1:9,:)=[];

v_gen=geneid_all;
v_data=data_orig_all_sorted;

%% HERE: find spatial genes(old)
%{
gen_list=[];
val_list=[];
samp_list=[];
bin_list=[];
rnk_list=[];
clear vis_info
all_gv=zeros(length(v_gen),length(all_vis_id));
for vii=1:length(all_vis_id)
    vii
    % find specifc visum
    curr_v=all_vis_id(vii);% example vis name
    v_id=find(string(sample_sorted)==curr_v);
    xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
    xv=xyv(:,1);
    yv=xyv(:,2);
    in = find(sum(v_data(:,v_id)>0,2)>10);% a gene that is expressed in more than 10 cells
    n_data=v_data(:,v_id);% in
    n_gen=v_gen(:);% in
    all_s=zeros(length(n_gen),1);
    all_p=zeros(length(n_gen),1);
    %%
    for gnx=1:length(n_gen)% run on each scell type
        if sum(gnx~=in)==length(in)
            %         gnxi=v_gen(gnx) % for specific gene
            markergene=n_data(gnx,:);
            sigx=std(xv(markergene>0));
            sigy=std(yv(markergene>0));
            sigs=sqrt(sigx^2 + sigy^2);
            spos=sum(markergene>0)/length(markergene);

            all_s(gnx)=sigs;% xy-std
            all_p(gnx)=spos;% percentage of red-spots
        end
    end

    %%
    %     P = polyfit(all_p,all_s,1);
    %     xi=0:0.001:1;
    %     yi=P(2)+ P(1).*all_p;
    %     figure;
    %     s=scatter(all_p,all_s,'.k');
    %     hold on
    %     %     s=scatter(all_p,all_s,'.r');
    %     row = dataTipTextRow("Gene",n_gen);
    %     s.DataTipTemplate.DataTipRows(end+1) = row;
    %     plot(all_p,yi,'r')
    %     ratio=all_s./yi;
    n_edges=10;
    [Y,E] = discretize(all_p,n_edges);
    sxy=[];
    for ii=1:n_edges
        %     sorted_p=sort(all_p(Y==ii));
        ii

        inn =find(Y==ii);
        s5=inn(all_s(inn)<prctile(all_s(inn),3));
        s95=inn(all_s(inn)>prctile(all_s(inn),97));
        sx=[s5;s95];
        sy=[find(all_s(inn)<prctile(all_s(inn),3));find(all_s(inn)>prctile(all_s(inn),97))];
        [~,ss] = sort(all_s(inn),'descend');
        r = 1:length(all_s(inn));
        r(ss) = r/length(r);% rank
        gen_list=[gen_list; n_gen(sx)];
        sxy=[sxy;sx];
        rnk_list=[rnk_list; r(sy)'];
        bin_list=[bin_list;ii*ones(length(n_gen(sx)),1)];
        %     val_list=[val_list; all_s(sx)];
        %           scatter(all_p(sx),all_s(sx),'or','filled');
        %       text(all_p(sx),all_s(sx),n_gen(sx));
    end
    samp=vii*ones(length(sxy),1);
    samp_list=[samp_list;samp];
    %     for ii=1:length(all_p)
    %         ii
    %     if ratio(ii)<mean([1.3,all_p(ii)]) && all_p(ii)>0.05 &&  all_p(ii)<0.85 && yi(ii)> all_s(ii)
    %       scatter(all_p(ii),all_s(ii),'or','filled');
    %       text(all_p(ii),all_s(ii),n_gen(ii));
    %       gen_list=[gen_list; n_gen(ii)];
    %     end
    %     end
    %     vis_info{vii,1}=[string(n_gen),all_p,all_s,ratio];
    %     xlabel('%')
    %     ylabel('STD')
end

%%
[gen_listx,ci,xi]=unique(gen_list,'stable');
rnk_listx=zeros(length(gen_listx),1);
bin_listx=zeros(length(gen_listx),1);
samp_listx=zeros(length(gen_listx),1);

for i=1:length(gen_listx)
    rid=find((xi==i));
    if median(rnk_list(rid))>0.97
        [mv,mi]=max(rnk_list(rid));
    else
        [mv,mi]=min(rnk_list(rid));

    end
    rnk_listx(i)=rnk_list(rid(mi));
    bin_listx(i)=bin_list(rid(mi));
    samp_listx(i)=samp_list(rid(mi));

end
t_list=table([string(gen_listx),rnk_listx,bin_listx,samp_listx]);
writetable(t_list,'gen_lst.csv');
%% get subid of binning
subidx=(bin_listx>=2 & bin_listx<=6);
subid=(bin_list>=2 & bin_list<=6);
gen_listx=gen_listx(subidx);
rnk_listx=rnk_listx(subidx);
bin_listx=bin_listx(subidx);
samp_listx=samp_listx(subidx);
gen_list=gen_list(subid);
rnk_list=rnk_list(subid);
bin_list=bin_list(subid);
samp_list=samp_list(subid);

%}

%% HERE: find spatial genes (new)
% gen_list=[];
% val_list=[];
% samp_list=[];
% bin_list=[];
% rnk_list=[];
% clear vis_info
all_gv=zeros(length(v_gen),length(all_vis_id));
all_gvx=zeros(length(v_gen),length(all_vis_id));
all_gvy=zeros(length(v_gen),length(all_vis_id));
all_mx=zeros(length(v_gen),length(all_vis_id));
all_my=zeros(length(v_gen),length(all_vis_id));

for vii=1:length(all_vis_id)
    vii
    % find specifc visum
    curr_v=all_vis_id(vii);% example vis name
    v_id=find(string(sample_sorted)==curr_v);
    xyv=cell2mat(bar_ar_sorted(v_id,[2 3]));
    xv=xyv(:,1);
    yv=xyv(:,2);
    %     MX= mean(xv);
    %     MY=mean(yv);
    %     MXY=mean([MX,MY]);
    in = find(sum(v_data(:,v_id)>0,2)>10);% a gene that is expressed in more than 10 cells
    n_data=v_data(:,v_id);% in
    n_gen=v_gen;% in
    all_s=zeros(length(n_gen),1);
    all_sx=zeros(length(n_gen),1);
    all_sy=zeros(length(n_gen),1);
    all_p=zeros(length(n_gen),1);
    %%
    for gnx=1:length(n_gen)% run on each scell type

        %         if sum(gnx~=in)==length(in)
        %         gnxi=v_gen(gnx) % for specific gene
        markergene=n_data(gnx,:);
        sigx=std(xv(markergene>0));
        sigy=std(yv(markergene>0));
        Mx=mean(xv(markergene>0));
        My=mean(yv(markergene>0));
        sigs=sqrt(sigx^2 + sigy^2);
        spos=sum(markergene>0)/length(markergene);

        all_s(gnx)=sigs;% xy-std
        mvisx=mean(xv);
        mvisy=mean(yv);
        all_mx(gnx,vii)=(Mx-mvisx)/(max(xv)-min(xv));
        all_my(gnx,vii)=(My-mvisy)/(max(yv)-min(yv));
        all_sx(gnx)=sigx;% x-std
        all_sy(gnx)=sigy;% xy-std
        all_p(gnx)=spos;% percentage of red-spots
        %         end
    end
   
    n_edges=10;
    [Y,E] = discretize(all_p,n_edges);
    for ii=2:6% take form, 0.2 :0.6
        %     sorted_p=sort(all_p(Y==ii));
        ii

        inn =find(Y==ii);
        % xy
        s5=inn(all_s(inn)<prctile(all_s(inn),2));
        s95=inn(all_s(inn)>prctile(all_s(inn),97));
        sxy=[s5;s95];
        %         all_gv(s95,vii)=-ii;
        all_gv(s5,vii)=ii;
        %         all_z(s5,vii)=all_s(s5)-Mx;%ii;
        % y
        s5=inn(all_sy(inn)<prctile(all_sy(inn),2));
        s95=inn(all_sy(inn)>prctile(all_sy(inn),97));
        sy=[s5;s95];
        %         all_gvy(s95,vii)=-ii;
        all_gvy(s5,vii)=ii;
        %         all_zy(s5,vii)=all_sy(s5);%ii;
        % x
        s5=inn(all_sx(inn)<prctile(all_sx(inn),2));
        s95=inn(all_sx(inn)>prctile(all_sx(inn),97));
        sx=[s5;s95];
        %         all_gvx(s95,vii)=-ii;
        all_gvx(s5,vii)=ii;
        %         all_zx(s5,vii)=all_sx(s5);%ii;
    end
    figure('color','w');
    scatter(all_p,all_s,10,'ok','filled');
    hold on;
    scatter(all_p(sx),all_s(sx),20,'or','filled');
    text(all_p(sx),all_s(sx),n_gen(sx));
    xlabel('%')
    ylabel('STD')
    
end
%% plot heatmap+bars of top20 down20
%% X(ML) +custom colorbar
[vx,ix]=sort(sum(all_gvx>0,2),'descend');% x
[vy,iy]=sort(sum(all_gvy>0,2),'descend'); % y
[vxy,ixy]=sort(sum(all_gv>0,2),'descend'); % xy
topg=100;
onoff=logical(all_gv(ixy(1:topg),:)>0);

hm=zeros(topg,16);
%%%%%%%%%%%
% x
genx=string(v_gen(ix(1:topg)));
genxy=string(v_gen(ixy(1:topg))); % choose genes of xy
hm(:,:)=all_mx(ixy(1:topg),:);%all_mx(gi,:);%all_gvx
hm=hm.*onoff;
nanid=isnan(hm);
hm(nanid)=0;
hmx=hm;
% for j=1:length(genxy)
%     gi=find(v_gen==genxy(j));
%     hm(j,:)=all_gvx(gi,:);%all_mx(gi,:);%all_gvx
% end
% sort genes
xx=logical(hm>0);
zmat=ones(size(hm)).*[1:8,1:8];
xmat=zeros(size(hm));
xmat(xx)=zmat(xx);
itrlv=interleave(1:8,9:16);
% [iv,id] = sort(sum(xmat,2) ./ sum(xmat~=0,2),'ascend');
zx = linkage([hmx,hmy],'ward','correlation');
Dx = pdist([hmx,hmy],'correlation');
leaforder = optimalleaforder(zx,Dx);
figure('color','w');
% subplot(1,2,1)
imagesc(hm(leaforder,itrlv));
% axis equal;axis tight;
% START of CUSTOM COLORBAR
L=10; %  number of data points
indexValue = 0;     % value for which to set a particular color
topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
indexColor = [0.9 0.9 0.9];       % color for indexed data value (white = [1 1 1])
bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
% Calculate where proportionally indexValue lies between minimum and
% maximum values
largest = max(max(hm));
smallest = min(min(hm));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
    linspace(bottomcolor(2),indexColor(2),100*index)',...
    linspace(bottomcolor(3),indexColor(3),100*index)'];
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
    linspace(indexColor(2),topColor(2),100*(L-index))',...
    linspace(indexColor(3),topColor(3),100*(L-index))'];
customCMap = [customCMap1;customCMap2];  % Combine colormaps
colormap(customCMap)
% psudo = pcolor(hm);
colorbar
set(gca, 'XTick',1:16, 'XTickLabel',interleave(1:8,1:8))
set(gca ,'YTick',1:length(genxy), 'YTickLabel',genxy(leaforder))
set(gca,'xtick',[])
% xtickangle(45)
title('ML')
xline([2.5:2:16],'k','LineWidth',5);
xline([1.5:1:16],'k','LineWidth',1);
yline([1.5:1:topg],'k','LineWidth',1);

% grid on
% cmap=redblue(256);
% colormap(redblue(256))
% colormap;
% colorbar;
%%%%%%%%%%%
% add sperating lines for puplation genes 
T = cluster(zx,'maxclust',15);% change number of clusters here 
dt=diff(T(leaforder))~=0;
vt=1:length(dt);
yline(vt(dt),'k','LineWidth',5)
hold on 
zhm=(hm(leaforder,itrlv)~=0);
idg=1:topg;
for d=size(hm,2):-1:1
scatter(d,idg(zhm(:,d)),'k','filled')
end
% xy

%% Y(DV) +custom colorbar
% y
hm=zeros(topg,16);
genxy=string(v_gen(ixy(1:topg))); % choose genes of xy
% geny=string(v_gen(iy(1:20))); choose genes of  y

hm(:,:)=all_my(ixy(1:topg),:);%all_mx(gi,:);%all_gvx
hm=hm.*onoff;
nanid=isnan(hm);
hm(nanid)=0;
hmy=hm;
% for j=1:length(genxy)
%     gi=find(v_gen==genxy(j));
%     hm(j,:)=all_my(gi,:);%all_gvy
% end
% sort genes
xx=logical(hm>0);
zmat=ones(size(hm)).*[1:8,1:8];
xmat=zeros(size(hm));
xmat(xx)=zmat(xx);
itrlv=interleave(1:8,9:16);
zx = linkage([hmx,hmy],'ward','correlation');
Dx = pdist([hmx,hmy],'correlation');
leaforder = optimalleaforder(zx,Dx);
% [iv,id] = sort((sum(xmat,2) ./ sum(xmat~=0,2)),'ascend');

figure('color','w');
% subplot(1,2,2)
imagesc(hm(leaforder,itrlv));% instead (id,itrlv)
% axis equal;axis tight;
% set(gca, 'XTick',1:16, 'XTickLabel',interleave(1:8,1:8))
% set(gca ,'YTick',1:length(genxy), 'YTickLabel',genxy)
% set(gca,'xtick',[])
% xtickangle(45)
title('DV')
xline([2.5:2:16],'k','LineWidth',5);
xline([1.5:1:16],'k','LineWidth',1);
yline([1.5:1:topg],'k','LineWidth',1);
% cmap=summer(256);
% START of CUSTOM COLORBAR
L=10; %  number of data points
indexValue = 0;     % value for which to set a particular color
topColor = [0 1 0];         % color for maximum data value (red = [1 0 0])
indexColor = [0.9 0.9 0.9];       % color for indexed data value (white = [1 1 1])
bottomcolor = [1 1 0];      % color for minimum data value (blue = [0 0 1])
% Calculate where proportionally indexValue lies between minimum and
% maximum values
largest = max(max(hm));
smallest = min(min(hm));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
    linspace(bottomcolor(2),indexColor(2),100*index)',...
    linspace(bottomcolor(3),indexColor(3),100*index)'];
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
    linspace(indexColor(2),topColor(2),100*(L-index))',...
    linspace(indexColor(3),topColor(3),100*(L-index))'];
customCMapx = [customCMap1;customCMap2];  % Combine colormaps
colormap(customCMapx)
% psudox = pcolor(hm);
colorbar

set(gca, 'XTick',1:16, 'XTickLabel',interleave(1:8,1:8))
set(gca ,'YTick',1:length(genxy), 'YTickLabel',genxy(leaforder))
set(gca,'xtick',[])
set(gca,'Position',[0.2042    0.1100    0.6063    0.8150])
title('DV')
% colorbar;

%%%%%%%%%%%
% add sperating lines for puplation genes 
T = cluster(zx,'maxclust',15);% change number of clusters here 
dt=diff(T(leaforder))~=0;
vt=1:length(dt);
yline(vt(dt),'k','LineWidth',5)
% add siginficance dots 
hold on 
zhm=(hm(leaforder,itrlv)~=0);
idg=1:topg;
for d=size(hm,2):-1:1
scatter(d,idg(zhm(:,d)),'k','filled')
end
% xy
%{ hm=zeros(20,16);

% genxy=string(v_gen(ixy(1:20)));
% for j=1:length(genxy)
%     gi=find(v_gen==genxy(j));
%     hm(j,:)=all_gv(gi,:);%all_gv
% end
% % sort genes
% xx=logical(hm>0);
% zmat=ones(size(hm)).*[1:8,1:8];
% xmat=zeros(size(hm));
% xmat(xx)=zmat(xx);
% itrlv=interleave(1:8,9:16);
% [iv,id] = sort(sum(xmat,2) ./ sum(xmat~=0,2),'ascend');
%
% subplot(1,4,3)
% imagesc(hm(id,itrlv));
% axis equal;axis tight;
% set(gca, 'XTick',1:16, 'XTickLabel',interleave(1:8,1:8))
% set(gca ,'YTick',1:length(genxy), 'YTickLabel',genxy(id))
% xtickangle(45)
% title('ML-DV')
% % cmap=redblue(256);
% colormap(summer)
% colormap;
% % colorbar;
%
% %%%%%%%%%%%
% % ap
% hm=zeros(20,16);
% % xy=logical(all_gv>0);
% % zmat=ones(size(all_gv)).*[1:8,1:8];
% % ymat=zeros(size(all_gv));
% ap=logical(all_gv>0);
% stend_all=zeros(length(v_gen),2);
% for i=1:length(v_gen)
% istr=strfind([0 ap(i,:)], [0 1])-1;  %gives indices of beginning of groups
% iend=strfind([ap(i,:) 0], [1 0]);    %gives indices of end of groups
% try
% stend_all(i,:)=[istr,iend];
% end
% end
% delta=stend_all(:,2)-stend_all(:,1);
% [vap,iap] = sort(delta,'descend');
% genxy=string(v_gen(iap(1:20)));
% for j=1:length(genxy)
%     gi=find(v_gen==genxy(j));
%     hm(j,:)=all_z(gi,:);
% end
% % sort genes
%
% xx=logical(hm>0);
% zmat=ones(size(hm)).*[1:8,1:8];
% xmat=zeros(size(hm));
% xmat(xx)=zmat(xx);
% itrlv=interleave(1:8,9:16);
% [iv,id] = sort(sum(xmat,2) ./ sum(xmat~=0,2),'ascend');
%
%
% subplot(1,4,4)
% imagesc(hm(id,itrlv));
% axis equal;axis tight;
% set(gca, 'XTick',1:16, 'XTickLabel',interleave(1:8,1:8))
% set(gca ,'YTick',1:length(genxy), 'YTickLabel',genxy(id))
% xtickangle(45)
% title('AP')
% cmap=[[1,1,1];summer(10)];
% colormap(cmap);
% colormap;
% colorbar;
%}%
% colormap(summer);
% sgtitle("stereo gene expression- 0.2:0.6")
%% create vis section  colormap
hf1=figure('color','w');
ax_link=[];
ML=1;
ha = @(m,n,p) subtightplot (m, n, p,[.001 .001],[.001 .001],[.001 .001]);
for vii=1:length(all_vis_id)
    %% find specifc visum
   
    vii
    curr_v=all_vis_id(vii);% example vis name
    v_id=find(string(sample_sorted)==curr_v);
    xv=cell2mat(bar_ar_sorted(v_id,2));
    yv=cell2mat(bar_ar_sorted(v_id,3));
    %%
if ML==1
    dots=round((max(xv)-min(xv))/60); %calculate number of dots in ecah X | Y 
    n_edges=dots;
    idots=round(linspace(1,length(customCMap),dots));
    clmp=customCMap(idots,:);   
    [Y,E] = discretize(xv,n_edges);
else % DV 
 dots=round((max(yv)-min(yv))/60); %calculate number of dots in ecah X | Y 
    n_edges=dots;
    idots=round(linspace(1,length(customCMapx),dots));
    clmp=customCMapx(idots,:);   
    [Y,E] = discretize(yv,n_edges);

end
%     markergene=log(data_orig_all_sorted(gnxi,:)+1);
   
        ax(vii)=ha(2,8,vii);


    for dis=1:dots
        scatter(xv(Y==dis),-yv(Y==dis),30,clmp(dis,:),'filled');
        hold on 
    end
        axis equal;
        axis tight;
        axis off;
        ax_link=[ax_link,ax(vii)];
        camroll(180)

end       % vis type

% resacale
%     sgtitle(v_gen(gnxi),'Interpreter','none')
% namepdf=char(v_gen(gnxi));
try
    allYLim = get(ax_link, {'YLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'YLim', [min(allYLim), max(allYLim)]);
    allYLim = get(ax_link, {'XLim'});
    allYLim = cat(2, allYLim{:});
    set(ax_link, 'XLim', [min(allYLim), max(allYLim)]);
end
% save as png

%     eval(['export_fig ',[direct,'/',namepdf],'.pdf -r 600']);
%     set(gcf,'Renderer','OpenGL');
%     set(gcf,'PaperPositionMode','auto');
set(gcf,'Renderer','OpenGL');
set(gcf,'PaperPositionMode','auto');
directx=[direct,'/VIS/'];
% % %     save2png([directx,namepdf,'_visium'],gcf,300)
%     close(hf1)

%% plot specific gene
close all
scoop=["GJA1"];
figure('color','w');
px=numSubplots(length(scoop));
for gg=1:length(scoop)
    gnxi=find(n_gen==scoop(gg)) % for specific gene
    markergene=n_data(gnxi,:);
    scatter(xv,yv,50,[0.8 0.8 0.8],'filled');
    hold on
    scatter(xv(markergene>0),yv(markergene>0),50,[1 0 0],'filled');
    title(scoop(gg))
    axis off; axis equal; axis tight;
end
% sgtitle("vis1")
%%