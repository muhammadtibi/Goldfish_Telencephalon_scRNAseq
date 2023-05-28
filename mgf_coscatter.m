%% calculate and plot co-scatters of exp. and per. of genes of both species 
clear all
close all
clc
cd('/data/Technion_analysis/goldfish/scRNAseq_gf')
load('/bigdata/all_mouse_dataset/mn_markeragg.mat','mn_markeragg','mn_clagg')
load("n_mg_10x.mat","all_name","all_data","geneid","flag_mgf","flag_g76","flag_mgfg76","cellid","sampleid",'n_mg','n_m','n_g','flag_rgn','flag_loc','flag_gen')
%% look at both datsets: co-scatter plot -log mean expression 
md=log2(mean(all_data(:,flag_mgf==1),2)+1);
gd=log2(mean(all_data(:,flag_mgf==2),2)+1);
figure;
scatter(md,gd)
%% total: co-scatter plot -percent expression 
% id=find(strcmp(geneid,'gad2'))
md=mean(all_data(:,flag_mgf==1)>0,2);
gd=mean(all_data(:,flag_mgf==2)>0,2);
%% total: mean expression 
md=log2(mean(all_data(:,flag_mgf==1),2)+1);
gd=log2(mean(all_data(:,flag_mgf==2),2)+1);
%%
md_all=zeros(length(geneid),length(n_m));
for j=1:length(n_m)
yy=all_name==n_m(j);
md_all(:,j)=mean(all_data(:,flag_mgf==1 & yy)>0,2);
end
gd_all=zeros(length(geneid),length(n_g));
for j=1:length(n_g)
yy=all_name==n_g(j);
gd_all(:,j)=mean(all_data(:,flag_mgf==2 & yy)>0,2);
end 
%%
smgd_all=zeros(length(n_g),length(n_m));
rt=0.5;

for i=1:length(n_g)
gd=gd_all(:,i);
i

smgd=sum(md_all>rt & gd>rt,1);
smgd_all(i,:)=smgd;
% figure;
% scatter(md,gd)

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
imagesc(smgd_all,[99 100]);
ylabel('Goldfish')
xlabel('Mouse')


numm=1:numel(n_m);
numg=1:numel(n_g);
n_gx= regexprep(n_g,'_','-');
n_gx= regexprep(n_gx,'g-','');
set(gca, 'YTick',numg, 'YTickLabel',n_gx)
set(gca ,'XTick',numm, 'XTickLabel',n_m)
xtickangle(45)
%% plot scatter
figure('color','w')
j=1;
i=1;
rt=0.5;
scatter(md_all(:,j),gd_all(:,i),'.')
xline(0.5,'r','LineWidth',3)
yline(0.5,'r','LineWidth',3)
xlabel(n_m(j))
ylabel(n_gx(i))
title(num2str(sum(md_all(:,j)>rt & gd_all(:,i)>rt)))
%% save table
Tx=[n_gx';smgd_all'];
zz=(['Mouse\Goldfish';n_m,['Marker';geni],['Location';loci],['Region';rgni]]);
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
n_gx= regexprep(n_g,'_','-');
n_gx= regexprep(n_gx,'g-','');
set(gca, 'YTick',numg, 'YTickLabel',n_gx)
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