close all; clear all;  clc;
warning off;

testname = 'test06';
regionname = 'pollock_egg3';
earlyyear=[1983:1987];
lateyear=[1988:1992];
inputmonth = [1,2]; % % put month which you want to plot [month month ...]
allyear =[1983:1992];
refyear =[1983:1987];
caxisval=[0 150];
caxisval_diff=[-70 70];
addpath(genpath(['C:\Users\User\Dropbox\source\matlab\function\']))
[dropboxpath, erorr_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([dropboxpath, '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\MICT_pollack\paper\subroutine\']))
[byrmap3, error_status] = Func_0009_get_colormaps('byr3', dropboxpath);
[byrmap, error_status] = Func_0009_get_colormaps('byr2', dropboxpath);  
[refpolygon, lonlat, error_status] = Func_0007_get_polygon_data_from_regionname(regionname);

param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_', regionname, '.m'];
% figrawdir =strcat('Z:\내 드라이브\research\Ph_D_course\Particle tracking of the walleye pollock\figure\paper\'); % % where figure files will be saved
figrawdir =strcat('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\'); % % where figure files will be saved
filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
LTRANS_testname='Pollock6';

tifname = [figrawdir, 'Fig10_03.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0]; % right, up, width, height
correction_large_fig=[0,0,0,-0.28]; % right, up, width, height


temp_checktime=15;
% % %  83-87 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear egg_mask comb_egg_mask

for yearij = 1:length(earlyyear)
    tempyear = earlyyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));

        lastday_m=size(egg_mask,3);
        if (exist('comb_egg_mask')==0)
            comb_egg_mask=egg_mask;
            temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
            temp_surf(isnan(temp_surf))=50000;
            temp_surf(temp_surf<50000)=NaN;
            RCM_land=temp_surf;
            RCM_land(temp_surf==50000)=1;
            temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
            RCM_ocean=temp_surf;
            RCM_ocean(isfinite(temp_surf))=1;
            lon_RCM = ncread(ncname, 'lon_rho');
            lat_RCM = ncread(ncname, 'lat_rho');
        
            [indw, inde, inds, indn]=Func_0012_findind_Y(0.05, [127 130 38 38], lon_RCM, lat_RCM, 1);
            for j=1:size(RCM_land,2)
                coastal_grid_ind(j)=max(find(RCM_land(1:inde,j)==1))
                for i=1:size(RCM_land,1)
% % % %                     for raw lat
%                     coastal_distance(i,j) = m_lldist( ...
%                         [lon_RCM(coastal_grid_ind(j),j), lon_RCM(i,j)], ...
%                         [lat_RCM(coastal_grid_ind(j),j), lat_RCM(i,j)]);

% % %                     for fixed lat (38)
                    coastal_distance(i,j) = m_lldist( ...
                        [lon_RCM(coastal_grid_ind(j),j), lon_RCM(i,j)], ...
                        [lat_RCM(coastal_grid_ind(j),inds), lat_RCM(i,inds)]);
                end
            end
            coastal_distance=coastal_distance.*RCM_ocean;
            coastal_distance(1:indw,:)=NaN;
        else
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
        end
    end
end
% pcolor(coastal_distance'); shading flat; colorbar;
lon_rho = ncread(ncname, 'lon_rho');
lat_rho = ncread(ncname, 'lat_rho');
sum_data = sum(comb_egg_mask,3, 'omitnan');
sum_data(sum_data==0)=NaN;
testnameind=1;

uniq_coast_dist=unique(coastal_distance);
uniq_coast_dist=uniq_coast_dist(isfinite(uniq_coast_dist));
sum_data_with_dist(1:length(uniq_coast_dist))=0;
for k=1:length(uniq_coast_dist)
    refdist=uniq_coast_dist(k);
    for i=1:size(sum_data,1)
        for j=1:size(sum_data,2)
            if (coastal_distance(i,j)==refdist && isfinite(sum_data(i,j)))
                sum_data_with_dist(k)=sum_data_with_dist(k)+sum_data(i,j);
            end
        end
    end
end
% bar(uniq_coast_dist, sum_data_with_dist)
sum_data_with_dist_early_15d=sum_data_with_dist;

% % projection
testnameind=1;
sb1=subplot(4,8,[2 3 10 11]);  % Auto-fitted to the figure.
% sb1=subplot(2,8,[2 3]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

bar(uniq_coast_dist, sum_data_with_dist)
ylabel('# of individuals')
xlabel('Distance (km)')
xlim([0 500])
ylim([0 2000]); 
grid minor;
set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(300, 1500, '(A) 83-87', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(300, 1800, '15 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XAxisLocation', 'top')
set(ax{testnameind,1}, 'XTick',[0  200 400])
set(ax{testnameind,1}, 'XTicklabel', {'0', '200', '400'})


% % %  88-92 (late period, period II) 15d
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear egg_mask comb_egg_mask

for yearij = 1:length(lateyear)
    tempyear = lateyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));

        lastday_m=size(egg_mask,3);
        if (exist('comb_egg_mask')==0)
            comb_egg_mask=egg_mask;
            temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
            temp_surf(isnan(temp_surf))=50000;
            temp_surf(temp_surf<50000)=NaN;
            RCM_land=temp_surf;
            RCM_land(temp_surf==50000)=1;
            temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
            RCM_ocean=temp_surf;
            RCM_ocean(isfinite(temp_surf))=1;
            lon_RCM = ncread(ncname, 'lon_rho');
            lat_RCM = ncread(ncname, 'lat_rho');
        
            [indw, inde, inds, indn]=Func_0012_findind_Y(0.05, [127 130 38 38], lon_RCM, lat_RCM, 1);
            for j=1:size(RCM_land,2)
                coastal_grid_ind(j)=max(find(RCM_land(1:inde,j)==1))
                for i=1:size(RCM_land,1)
% % % %                     for raw lat
%                     coastal_distance(i,j) = m_lldist( ...
%                         [lon_RCM(coastal_grid_ind(j),j), lon_RCM(i,j)], ...
%                         [lat_RCM(coastal_grid_ind(j),j), lat_RCM(i,j)]);

% % %                     for fixed lat (38)
                    coastal_distance(i,j) = m_lldist( ...
                        [lon_RCM(coastal_grid_ind(j),j), lon_RCM(i,j)], ...
                        [lat_RCM(coastal_grid_ind(j),inds), lat_RCM(i,inds)]);
                end
            end
            coastal_distance=coastal_distance.*RCM_ocean;
            coastal_distance(1:indw,:)=NaN;
        else
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
        end
    end
end
% pcolor(coastal_distance'); shading flat; colorbar;
lon_rho = ncread(ncname, 'lon_rho');
lat_rho = ncread(ncname, 'lat_rho');
sum_data = sum(comb_egg_mask,3, 'omitnan');
sum_data(sum_data==0)=NaN;
testnameind=1;

uniq_coast_dist=unique(coastal_distance);
uniq_coast_dist=uniq_coast_dist(isfinite(uniq_coast_dist));
sum_data_with_dist(1:length(uniq_coast_dist))=0;
for k=1:length(uniq_coast_dist)
    refdist=uniq_coast_dist(k);
    for i=1:size(sum_data,1)
        for j=1:size(sum_data,2)
            if (coastal_distance(i,j)==refdist && isfinite(sum_data(i,j)))
                sum_data_with_dist(k)=sum_data_with_dist(k)+sum_data(i,j);
            end
        end
    end
end
% bar(uniq_coast_dist, sum_data_with_dist)
sum_data_with_dist_late_15d=sum_data_with_dist;

% % projection
testnameind=2;
sb1=subplot(4,8,[4 5 12 13]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

bar(uniq_coast_dist, sum_data_with_dist)
xlabel('Distance (km)')

xlim([0 500])
ylim([0 2000]); 
grid minor;
set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(300, 1500, '(B) 88-92', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(300, 1800, '15 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XAxisLocation', 'top')
set(ax{testnameind,1}, 'XTick',[0  200 400])
set(ax{testnameind,1}, 'XTicklabel', {'0', '200', '400'})
set(ax{testnameind,1}, 'YTicklabel', [])

% % %  get difference (late-early period, period II-I)
sum_data_with_dist_diff_15d=sum_data_with_dist_late_15d-sum_data_with_dist_early_15d;

% % projection
testnameind=3;
sb2=subplot(4,8,[6 7 14 15]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.

pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb2); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

bar(uniq_coast_dist, sum_data_with_dist_diff_15d, 'r')
xlabel('Distance (km)')
ylabel('# of individuals')

xlim([0 500])
ylim([-700 700]); 
grid minor;
set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(300, 350, '(C) B-A', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(300, 560, '15 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XAxisLocation', 'top')
set(ax{testnameind,1}, 'XTick',[0  200 400])
set(ax{testnameind,1}, 'XTicklabel', {'0', '200', '400'})
set(ax{testnameind,1}, 'YAxisLocation', 'right')


% % % % % % % % % % 
% % % % % % % % % %  row 2
% % % % % % % % % %

correction_upper_fig = [0 0.31 0 0];

temp_checktime=30;
% % %  83-87 (early period, period I) 30d

run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear egg_mask comb_egg_mask

for yearij = 1:length(earlyyear)
    tempyear = earlyyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));

        lastday_m=size(egg_mask,3);
        if (exist('comb_egg_mask')==0)
            comb_egg_mask=egg_mask;
            temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
            temp_surf(isnan(temp_surf))=50000;
            temp_surf(temp_surf<50000)=NaN;
            RCM_land=temp_surf;
            RCM_land(temp_surf==50000)=1;
            temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
            RCM_ocean=temp_surf;
            RCM_ocean(isfinite(temp_surf))=1;
            lon_RCM = ncread(ncname, 'lon_rho');
            lat_RCM = ncread(ncname, 'lat_rho');
        
            [indw, inde, inds, indn]=Func_0012_findind_Y(0.05, [127 130 38 38], lon_RCM, lat_RCM, 1);
            for j=1:size(RCM_land,2)
                coastal_grid_ind(j)=max(find(RCM_land(1:inde,j)==1))
                for i=1:size(RCM_land,1)
% % % %                     for raw lat
%                     coastal_distance(i,j) = m_lldist( ...
%                         [lon_RCM(coastal_grid_ind(j),j), lon_RCM(i,j)], ...
%                         [lat_RCM(coastal_grid_ind(j),j), lat_RCM(i,j)]);

% % %                     for fixed lat (38)
                    coastal_distance(i,j) = m_lldist( ...
                        [lon_RCM(coastal_grid_ind(j),j), lon_RCM(i,j)], ...
                        [lat_RCM(coastal_grid_ind(j),inds), lat_RCM(i,inds)]);
                end
            end
            coastal_distance=coastal_distance.*RCM_ocean;
            coastal_distance(1:indw,:)=NaN;
        else
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
        end
    end
end
% pcolor(coastal_distance'); shading flat; colorbar;
lon_rho = ncread(ncname, 'lon_rho');
lat_rho = ncread(ncname, 'lat_rho');
sum_data = sum(comb_egg_mask,3, 'omitnan');
sum_data(sum_data==0)=NaN;

uniq_coast_dist=unique(coastal_distance);
uniq_coast_dist=uniq_coast_dist(isfinite(uniq_coast_dist));
sum_data_with_dist(1:length(uniq_coast_dist))=0;
for k=1:length(uniq_coast_dist)
    refdist=uniq_coast_dist(k);
    for i=1:size(sum_data,1)
        for j=1:size(sum_data,2)
            if (coastal_distance(i,j)==refdist && isfinite(sum_data(i,j)))
                sum_data_with_dist(k)=sum_data_with_dist(k)+sum_data(i,j);
            end
        end
    end
end
% bar(uniq_coast_dist, sum_data_with_dist)
sum_data_with_dist_early_30d=sum_data_with_dist;

% % projection
testnameind=4;
sb1=subplot(4,8,[18 19 26 27]);  % Auto-fitted to the figure.
% sb1=subplot(2,8,[10 11]);  % Auto-fitted to the figure.

pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0]; % right, up, width, height
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

bar(uniq_coast_dist, sum_data_with_dist)
ylabel('# of individuals')
xlabel('Distance (km)')
xlim([0 500])
ylim([0 2000]); 
grid minor;
set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(300, 1500, '(D) 83-87', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(300, 1800, '30 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XTick',[0  200 400])
set(ax{testnameind,1}, 'XTicklabel', {'0', '200', '400'})


% % %  88-92 (lat period, period II) 30d
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear egg_mask comb_egg_mask

for yearij = 1:length(lateyear)
    tempyear = lateyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));

        lastday_m=size(egg_mask,3);
        if (exist('comb_egg_mask')==0)
            comb_egg_mask=egg_mask;
            temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
            temp_surf(isnan(temp_surf))=50000;
            temp_surf(temp_surf<50000)=NaN;
            RCM_land=temp_surf;
            RCM_land(temp_surf==50000)=1;
            temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
            RCM_ocean=temp_surf;
            RCM_ocean(isfinite(temp_surf))=1;
            lon_RCM = ncread(ncname, 'lon_rho');
            lat_RCM = ncread(ncname, 'lat_rho');
        
            [indw, inde, inds, indn]=Func_0012_findind_Y(0.05, [127 130 38 38], lon_RCM, lat_RCM, 1);
            for j=1:size(RCM_land,2)
                coastal_grid_ind(j)=max(find(RCM_land(1:inde,j)==1))
                for i=1:size(RCM_land,1)
% % % %                     for raw lat
%                     coastal_distance(i,j) = m_lldist( ...
%                         [lon_RCM(coastal_grid_ind(j),j), lon_RCM(i,j)], ...
%                         [lat_RCM(coastal_grid_ind(j),j), lat_RCM(i,j)]);

% % %                     for fixed lat (38)
                    coastal_distance(i,j) = m_lldist( ...
                        [lon_RCM(coastal_grid_ind(j),j), lon_RCM(i,j)], ...
                        [lat_RCM(coastal_grid_ind(j),inds), lat_RCM(i,inds)]);
                end
            end
            coastal_distance=coastal_distance.*RCM_ocean;
            coastal_distance(1:indw,:)=NaN;
        else
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
        end
    end
end
% pcolor(coastal_distance'); shading flat; colorbar;
lon_rho = ncread(ncname, 'lon_rho');
lat_rho = ncread(ncname, 'lat_rho');
sum_data = sum(comb_egg_mask,3, 'omitnan');
sum_data(sum_data==0)=NaN;

uniq_coast_dist=unique(coastal_distance);
uniq_coast_dist=uniq_coast_dist(isfinite(uniq_coast_dist));
sum_data_with_dist(1:length(uniq_coast_dist))=0;
for k=1:length(uniq_coast_dist)
    refdist=uniq_coast_dist(k);
    for i=1:size(sum_data,1)
        for j=1:size(sum_data,2)
            if (coastal_distance(i,j)==refdist && isfinite(sum_data(i,j)))
                sum_data_with_dist(k)=sum_data_with_dist(k)+sum_data(i,j);
            end
        end
    end
end
% bar(uniq_coast_dist, sum_data_with_dist)
sum_data_with_dist_late_30d=sum_data_with_dist;

% % projection
testnameind=5;
sb2=subplot(4,8,[20 21 28 29]);  % Auto-fitted to the figure.
% sb1=subplot(2,8,[10 11]);  % Auto-fitted to the figure.

pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0]; % right, up, width, height
delete(sb2); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

bar(uniq_coast_dist, sum_data_with_dist)
xlabel('Distance (km)')
xlim([0 500])
ylim([0 2000]); 
grid minor;
set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(300, 1500, '(E) 88-92', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(300, 1800, '30 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XTick',[0  200 400])
set(ax{testnameind,1}, 'XTicklabel', {'0', '200', '400'})
set(ax{testnameind,1}, 'YTicklabel', [])

% % %  get difference (late-early period, period II-I)
sum_data_with_dist_diff_30d=sum_data_with_dist_late_30d-sum_data_with_dist_early_30d;

% % projection
testnameind=6;
sb2=subplot(4,8,[22 23 30 31]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.

pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb2); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

bar(uniq_coast_dist, sum_data_with_dist_diff_30d, 'r')
ylabel('# of individuals')
xlabel('Distance (km)')
xlim([0 500])
ylim([-700 700]); 
grid minor;
set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(300, 350, '(F) B-A', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(300, 560, '30 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XTick',[0  200 400])
set(ax{testnameind,1}, 'XTicklabel', {'0', '200', '400'})
set(ax{testnameind,1}, 'YAxisLocation', 'right')


% % figure save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y/3]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
print('-dtiff','-r500',tifname); 
RemoveWhiteSpace([], 'file', tifname);
hold off;
close all;

% save ('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\data_dryad\Data06_numpar_lon.mat', ...
%             'lat_1d', ...
%             'lon_sum_early_data_15d', 'lon_sum_late_data_15d', 'lon_sum_diff_data_15d', ...
%             'lon_sum_early_data_30d', 'lon_sum_late_data_30d', 'lon_sum_diff_data_30d')

