% % % 
% % % exp 1.
% % % 
clc; clear all; close all;

addpath(genpath('D:\MEPL\project\NWP\m_map'))

rname='D:\class\PhD\Objective_Data_Analysis\term_paper\ES_clim2.nc';
temp1=ncread(rname,'CLIM_TEMP');
u1=ncread(rname,'CLIM_U');
v1=ncread(rname,'CLIM_V');

rname='D:\class\PhD\Objective_Data_Analysis\term_paper\temp.nc';
temp2=ncread(rname,'ES_TEMP');
lon2=ncread(rname,'XT_OCEAN293_442');
lat2=ncread(rname,'YT_OCEAN281_470');
level2=ncread(rname,'ST_OCEAN1_33');
ilevel2 = -level2;


rname2='D:\class\PhD\Objective_Data_Analysis\term_paper\transport2.nc';
tr=ncread(rname2,'TR');
% y_tr=ncread(rname2,'Y_TR');

[Plg,Plt]=meshgrid(lon2,lat2);
for i=1:150
    for j=1:190
        for k=1:33
            for l=1:335
                if (temp2(i,j,k,l)< -1.0e+5)
                    temp2(i,j,k,l)=NaN;
                end
            end
        end
    end
end
% m_proj('UTM','long',[128 143],'lat',[33 52]);
% m_proj('mercator','lon',[128 143],'lat',[33 52]);
% m_pcolor(Plg,Plt,squeeze(temp2(1:150,1:190,1,335))');


%%% 120m whole temp
pcolor(Plg,Plt,squeeze(temp1(1:150,1:190,16,7))');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'Temp (^oC)';
xlabel('Longitude (^oE)','fontsize',27);
ylabel('Latitude (^oN)','fontsize',27);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mtemp.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');




%%%% 120m temp

pcolor(lon2(5:65),lat2(15:70),squeeze(temp1(5:65,15:70,16,7))');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'Temp (^oC)';
xlabel('Longitude (^oE)','fontsize',15);
ylabel('Latitude (^oN)','fontsize',15);
set(gca,'fontsize',10);
hold on;
[C hT]= contour(lon2(5:65),lat2(15:70),squeeze(temp2(5:65,15:70,16,335))', '-k','ShowText','on');
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mtemp_ksbcw.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');





% % EOF용 -> 129.5 : 131.5, 35.6
% %       129.7, 35 (point of warm water)
% % corr, 129.7, 35.6 (point of KSBCW)
% %       129.7, 37
% %       129, 38
% %       128.5, 39  (wonsan bay)
% %       129, 40
% %       130, 41
% %       132, 42 (vladivostok)
% % 
% % 울릉분지 130.5, 36
% %          130.5, 37.5
% % whole domain : 128 ~ 143, 33 ~ 52, 1~ 500, 1985 ~ 2011
% % size of temp2 : 150, 190, 33, 335

warm=squeeze(temp2(17,26,1,:));
ksbcw=squeeze(temp2(17,26,16,:));
ksbcw_2=squeeze(temp2(17,40,16,:));
ksbcw_3=squeeze(temp2(10,50,16,:));
wonsan_bay=squeeze(temp2(5,60,16,:));
wonsan_bay_2=squeeze(temp2(10,70,16,:));
wonsan_bay_3=squeeze(temp2(20,80,16,:));
vlad_surf=squeeze(temp2(40,90,1,:));
vlad_mid=squeeze(temp2(40,90,16,:));

for i=1:150
    for j=1:190
        corcoef1(i,j)=corr(ksbcw,squeeze(temp2(i,j,16,:)));
        corcoef1_2(i,j)=corr(ksbcw,squeeze(temp2(i,j,18,:)));
        corcoef1_3(i,j)=corr(ksbcw,squeeze(temp2(i,j,20,:)));
        corcoef1_4(i,j)=corr(ksbcw,squeeze(temp2(i,j,22,:)));
        corcoef1_5(i,j)=corr(ksbcw,squeeze(temp2(i,j,24,:)));
        corcoef2(i,j)=corr(ksbcw,squeeze(temp2(i,j,1,:)));
        corcoef3(i,j)=corr(tr',squeeze(temp2(i,j,16,:)));
    end
end

for j=1:120
    for k=1:33
        corcoef4(j,k)=corr(ksbcw,squeeze(temp2(17,j,k,:)));
        corcoef5(j,k)=corr(ksbcw,squeeze(temp2(22,j,k,:)));
        corcoef6(j,k)=corr(ksbcw,squeeze(temp2(27,j,k,:)));
        corcoef6_2(j,k)=corr(ksbcw,squeeze(temp2(32,j,k,:)));
        corcoef7(j,k)=corr(ksbcw,squeeze(temp2(12,j,k,:)));
        corcoef8(j,k)=corr(ksbcw,squeeze(temp2(7,j,k,:)));
    end
end


% % 120m corr
pcolor(Plg,Plt,corcoef1');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Longitude (^oE)','fontsize',10);
ylabel('Latitude (^oN)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mcorr.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');

% % 140m corr
pcolor(Plg,Plt,corcoef1_2');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Longitude (^oE)','fontsize',10);
ylabel('Latitude (^oN)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\140mcorr.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');


% % 160m corr
pcolor(Plg,Plt,corcoef1_3');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Longitude (^oE)','fontsize',10);
ylabel('Latitude (^oN)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\160mcorr.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');

% % 180m corr
pcolor(Plg,Plt,corcoef1_4');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Longitude (^oE)','fontsize',10);
ylabel('Latitude (^oN)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\180mcorr.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');

% % 200m corr
pcolor(Plg,Plt,corcoef1_5');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Longitude (^oE)','fontsize',10);
ylabel('Latitude (^oN)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\200mcorr.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');



% % cold water  & vertical section corr
pcolor(lat2(1:120)',ilevel2(1:33),corcoef4');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Latitude (^oN)','fontsize',10);
ylabel('depth(m)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mcorr_vert1.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
% axis equal;
axis tight;
saveas(gcf,filename,'tiff');

% % cold water  & vertical section corr2
pcolor(lat2(1:120)',ilevel2(1:33)',corcoef5');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Latitude (^oN)','fontsize',10);
ylabel('depth(m)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mcorr_vert2.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
% axis equal;
axis tight;
saveas(gcf,filename,'tiff');

% % cold water  & vertical section corr3
pcolor(lat2(1:120)',ilevel2(1:33)',corcoef6');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Latitude (^oN)','fontsize',10);
ylabel('depth(m)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mcorr_vert3.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
% axis equal;
axis tight;
saveas(gcf,filename,'tiff');

% % cold water  & vertical section corr3_2
pcolor(lat2(1:120)',ilevel2(1:33)',corcoef6_2');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Latitude (^oN)','fontsize',10);
ylabel('depth(m)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mcorr_vert3_2.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
% axis equal;
axis tight;
saveas(gcf,filename,'tiff');




% % cold water  & vertical section corr4
pcolor(lat2(1:120)',ilevel2(1:33)',corcoef7');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Latitude (^oN)','fontsize',10);
ylabel('depth(m)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mcorr_vert4.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
% axis equal;
axis tight;
saveas(gcf,filename,'tiff');

% % cold water  & vertical section corr5
pcolor(lat2(1:120)',ilevel2(1:33)',corcoef8');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
xlabel('Latitude (^oN)','fontsize',10);
ylabel('depth(m)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\120mcorr_vert5.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
% axis equal;
axis tight;
saveas(gcf,filename,'tiff');






% % surface corr
pcolor(Plg,Plt,corcoef2');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
xlabel('Longitude (^oE)','fontsize',10);
ylabel('Latitude (^oN)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\surf_corr.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');


% % transport corr
pcolor(Plg,Plt,corcoef3');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
xlabel('Longitude (^oE)','fontsize',10);
ylabel('Latitude (^oN)','fontsize',10);
set(gca,'fontsize',15);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\tr_corr.tif'
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');



startDate = datenum('01-15-1985');
endDate = datenum('11-15-2012');
xData = linspace(startDate,endDate,335);
grad=warm-ksbcw;
plot(xData,grad);
datetick('x','yyyy')
title_name = 'Difference of warm water and cold water';
title(title_name,'fontsize',30);
% xlabel('time','fontsize', 20);
ylabel('Temperature(^oC)','fontsize',20);
set(gca,'FontSize',10);
saveas(gcf,'D:\class\PhD\Objective_Data_Analysis\term_paper\dif_temp.tif','tiff')
close;

corr(ksbcw,ksbcw_2)
corr(ksbcw,ksbcw_2)
corr(ksbcw,ksbcw_3)
corr(ksbcw,wonsan_bay)
corr(ksbcw,wonsan_bay_2)
corr(ksbcw,wonsan_bay_3)
corr(ksbcw,vlad_surf)
corr(ksbcw,vlad_mid)
xcorr(ksbcw,vlad_mid,12,'coeff')
corr(ksbcw,grad)
plot(vlad_mid)
hold on
plot(ksbcw)


% % EOF

% sort data to one column
% % 129.5 ~ 131.5, 35.6, k=16
% % clear linetemp
% % for l=1:335
% %     for i=17:37
% %         linetemp((i-16)+(l-1)*21)=temp2(i,26,16,l);
% %     end
% % end

for i=17:37
    for l=1:335
        temp2anom(i-16,l)=temp2(i,26,16,l)-mean(temp2(i,26,16,:));
    end
end
% X=squeeze(temp2(17:37,26,16,:));
X=temp2anom;
[M N] = size(X);
 %solve eigen value problem
R = 1/N * X * X';
[V, D] = eig(R); %eigen vector , eigen value

%sort eigen values 
L = diag(D);L = sort(L,'descend'); %sorted eigen value
for i = 1:length(L);
    ind = find(L == L(i));
    E(:,i) = V(:,ind); %E : sorted eigen vector
end;

 % find PC 
Z = E' * X; %Z: PC

 %find variance / total variance
total_variance = sum(L) ;
mode_variance = abs(L)/total_variance*100; %각 mode별 설명정도

plot(mode_variance,'*','LineWidth',3);set(gca,'XTick',[1:1:M]);
title('Explained variance of Each Mode','fontsize',30);
xlabel('Mode number','fontsize', 20);
ylabel('variance/total variance (%)','fontsize',20);
set(gca,'FontSize',10);
filename='D:\class\PhD\Objective_Data_Analysis\term_paper\mode_var.tif'
hold off;
saveas(gcf,filename,'tiff');


%1) plot eigen vector(E)
%2) plot PC(Z)
K=2;
for i =1:K
    e = E(i,:);
    plot(e);
    title_name = strcat('Eigen Vector-','MODE',num2str(i));
    title(title_name,'fontsize',30);
    xlabel('station','fontsize', 20);ylabel('value(non-dimension)','fontsize',20);set(gca,'FontSize',10);
    saveas(gcf,strcat(title_name,'.tif'),'tiff')
    close;
end
 startDate = datenum('01-15-1985');
 endDate = datenum('11-15-2012');
 xData = linspace(startDate,endDate,335);
for i =1:K
    z = Z(i,:);
    plot(xData,z);
     datetick('x','yyyy')
    title_name = strcat('PC-','MODE',num2str(i));
    title(title_name,'fontsize',30);
    xlabel('time','fontsize', 20);ylabel('value(dimension)','fontsize',20);set(gca,'FontSize',10);
    saveas(gcf,strcat(title_name,'.tif'),'tiff')
    close;
end
