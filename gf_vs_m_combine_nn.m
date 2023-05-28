% combine gf & mm nn (newest)
%% gf_vs_m
clear all
close all
clc
nrmz=10000;
sp=1;

%% load mouse-new (my)
disp('mm-new')
load('/bigdata/all_mouse_dataset/loom_nn.mat','mnn_matrix','mnn_genes','mnn_region','mnn_marker','mnn_location','mnn_cl','mnn_cln','mnn_id','mnn_dev','mnn_ts','mnn_nt','mnn_t3');
m_nt=mnn_nt;
m_data=mnn_matrix';
m_name=mnn_cln;
m_genes=cellstr(mnn_genes);
m_sample=ones(length(mnn_cln),1);%mn_id;
m_region=mnn_region;
m_location=mnn_location;
m_marker=mnn_marker;
m_rank=string(mnn_t3);

%% random sampling mm
if sp==1
disp('sample mm')
n_m=unique(m_name,'stable');
y_all=[];
for rs=1:length(n_m)
    ig=find(m_name==n_m(rs));% for a cell cluster
    if  length(ig)>200
        rs
        y = randsample(ig,200); % here we take only XX cells of the real cell typoe cluster
        y_all=[y_all;y];
    else
        y_all=[y_all;ig];
    end

end
m_name=m_name(y_all);
m_region=m_region(y_all);
m_sample=m_sample(y_all);
m_rank=m_rank(y_all);
m_marker=m_marker(y_all);
m_location=m_location(y_all);
m_nt=m_nt(y_all);
m_data=m_data(:,y_all);
end

%% load gf- new
disp('gf-new')
load('/data/Technion_analysis/goldfish/scRNAseq_gf/gf_sc_new/non_normalize_all_data_15_8_2022.mat', 'g_geneid','g_data','g_ca');
% gdata =round(g_data./repmat(sum(g_data),length(g_data(:,1)),1)*nrmz);
nn=zeros(size(g_ca));
glut=contains(g_ca,"GLUT");
gaba=contains(g_ca,"GABA");
nn(~(glut+gaba))=1;
g_data=g_data(:,logical(nn));
g_name=g_ca(logical(nn));
g_genes=g_geneid;
g_sample=g_name;
%% random sampling gf
if sp==1
disp('sample gf')
y_all=[];
n_g=unique(g_name,'stable');
for rs=1:length(n_g)
    ig=find(g_name==n_g(rs));% for a cell cluster
    if  length(ig)>200
        rs
        y = randsample(ig,200); % here we take only XX cells of the real cell typoe cluster
        y_all=[y_all;y];
    else
        y_all=[y_all;ig];
    end

end
g_name=g_name(y_all);
g_ca=g_ca(y_all);
g_sample=g_sample(y_all);
g_data=g_data(:,y_all);
end

%% concatante MM & B mutual genes - for mouse bat
disp('conc m&g')

[inBoth,xi,xf] = intersect(upper(deblank(m_genes)), string(g_genes)) ;
m_data=m_data(xi,:);
m_data=double(round(m_data./repmat(sum(m_data),length(m_data(:,1)),1)*nrmz));
g_data=g_data(xf,:);
g_data=double(round(g_data./repmat(sum(g_data),length(g_data(:,1)),1)*nrmz));

geneid=inBoth;
%% neurotasmitor flag
disp('m g76 - nt flag')
gads=0;


ntf=ones(length(m_nt),1);% 0=> NN

flag_m=ntf;

%% create mouse m76 flag
disp('m g76')


nttp= m_rank=='Telencephalon projecting neurons';
ntti= m_rank=='Telencephalon interneurons';
nttx= m_rank=='Cholinergic, monoaminergic and peptidergic neurons';
ntty= m_rank=='Di- and mesencephalon neurons';
nttz= contains(m_region,"Pons") | contains(m_region,"Medulla") | contains(m_region,"Spinal cord")| contains(m_region,"Midbrain");

% nti=ones(length(m_rank),1);
% nti(ntgd | ntach | ntv)=0;% 0 = filter out

%
nti=ones(length(m_rank),1);
% nti((nttp | ntti | nttx | ntty) & ~nttz)=1;% 0 = filter in 
% choose INH /GLU or uncomment 
% ntf=contains(m_name,"INH");
% flag_m=ntf;

ntx=logical(nti);
%% create gf flag
disp('gf g76')

flag_gf=ones(length(g_name),1);% 0=> NN


%% concatante
disp('conc all')
all_name=[m_name(logical(ntx));g_name(logical(flag_gf))];
all_data=[m_data(:,logical(ntx)),g_data(:,logical(flag_gf))];
sampleid=[m_sample(logical(ntx));g_sample(logical(flag_gf))];
n_m = unique(m_name(logical(ntx)),'stable');
n_g = unique(g_name(logical(flag_gf)),'stable');
%% decoy ids & flags
disp('dec flg')
% % sampleid=cell(length(all_c),1);
cellid=cell(length(all_name),1);
flag_mgf=2.*ones(length(all_name),1);% gf=2 | b=1
flag_mgf(1:sum(ntx))=1;
flag_g76=[flag_m(ntx);flag_gf(logical(flag_gf))]; % 0 1 2
flag_mgfg76=[flag_m(ntx);flag_gf(logical(flag_gf))+3];% +3 to make indexes 0 1 2  | 3 4
n_mg=[n_m;n_g];
flag_rgn=[m_region(ntx);strings(sum(logical(flag_gf)),1)];
flag_loc=[m_location(ntx);strings(sum(logical(flag_gf)),1)];
flag_gen=[m_marker(ntx);strings(sum(logical(flag_gf)),1)];

%% save
disp('save')
cd('/data/Technion_analysis/goldfish/scRNAseq_gf/comparative/nn_gf_mouse/')
save("nn_mg_10x.mat","all_name","all_data","geneid","flag_mgf","flag_g76","flag_mgfg76","cellid","sampleid",'n_mg','n_m','n_g','flag_rgn','flag_loc','flag_gen','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at both datsets: co-scatter plot -log mean expression 
md=log2(mean(all_data(:,flag_mgf==1),2)+1);
gd=log2(mean(all_data(:,flag_mgf==2),2)+1);
figure;
scatter(md,gd)
%% co-scatter plot -percent expression 
% id=find(strcmp(geneid,'gad2'))
md=(mean(all_data(:,flag_mgf==1)>0,2));
gd=(mean(all_data(:,flag_mgf==2)>0,2));
%% mean expression 
md=log2(mean(all_data(:,flag_mgf==1),2)+1);
gd=log2(mean(all_data(:,flag_mgf==2),2)+1);
%%

smgd_all=zeros(length(n_g),length(n_m));
for i=1:length(n_g)
xx=all_name==n_g(i);
gd=log2(mean(all_data(:,flag_mgf==2 & xx),2)+1);
mgd=max(gd)/2;
i
for j=1:length(n_m)
yy=all_name==n_m(j);
md=log2(mean(all_data(:,flag_mgf==1 & yy),2)+1);
smgd=sum(md>mgd & gd>mgd);
smgd_all(i,j)=smgd;
% figure;
% scatter(md,gd)
end
end
%% heatmap of smgd_all
figure('Color','w');
cmap=redblue(256);
colormap(cmap)
% u_cx(55,:)=[]; % derbalek ana kemtehn I remove it because its problamatic
% u_cx(:,55)=[];
% u_cx=u_cx(1:length(n_mx),length(n_mx)+1:end);
% % Yup = prctile(u_cx(:),80);% thresholding imagesc 
% % Ydn = prctile(u_cx(:),50);
% Z1 = linkage(u_cx,'ward','euclidean');
% % leaforder1 = optimalleaforder(Z1,pdist(u_cx));
% leaforder1 = 1:length(n_m);
% Z2 = linkage(u_cx','ward','euclidean');
% leaforder2 = optimalleaforder(Z2,pdist(u_cx'));
imagesc(smgd_all,[0 20]);
ylabel('Goldfish')
xlabel('Mouse')


numm=1:numel(n_m);
numg=1:numel(n_g);
set(gca, 'YTick',numg, 'YTickLabel',n_g)
set(gca ,'XTick',numm, 'XTickLabel',n_m)
xtickangle(45)
%% save table
Tx=[n_g(leaforder2)';u_cx];
zz=([['Mouse\Goldfish';n_m(leaforder1)],['Marker';geni],['Location';loci],['Region';rgni]]);
Ty=[zz,Tx];
Tt=table(Ty);
writetable(Tt,'Tm2gf.csv')

colormap;
%%
[vmax,imax]=max(smgd_all,[],2);
z_all=zeros(size(smgd_all));
for i=1:length(n_g)
z_all(i,imax(i))=1;
end
figure('color','w');
imagesc(z_all,[0,1]);
ylabel('Goldfish')
xlabel('Mouse')
numm=1:numel(n_m);
numg=1:numel(n_g);
set(gca, 'YTick',numg, 'YTickLabel',n_g)
set(gca ,'XTick',numm, 'XTickLabel',n_m)
xtickangle(45)

cmap=redblue(256);
colormap(cmap)
colormap;
%
% n_mby=n_mbx;
% n_mby(55)=[];
% Tx=[n_g(leaforder2)';u_cx];
% zz=([['Mouse\Goldfish';n_m(leaforder1)],['Marker';geni],['Location';loci],['Region';rgni]]);
% Ty=[zz,Tx];
% Tt=table(Ty);
% writetable(Tt,'Tm2gf.csv')

%% co-scatter plot -percent expression 
% id=find(strcmp(geneid,'gad2'))
md=(mean(all_data(:,flag_mgf==1)>0,2));
gd=(mean(all_data(:,flag_mgf==2)>0,2));
figure;
scatter(md,gd)