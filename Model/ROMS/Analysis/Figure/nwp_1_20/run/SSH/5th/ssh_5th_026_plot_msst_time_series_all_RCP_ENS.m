close all; clear all;  clc;
warning off;

RCM_rcp26_testnames = {'test61', 'test62', 'test63', 'test64', 'ens09'};
RCM_rcp45_testnames = {'test57', 'test58', 'test59', 'test60', 'ens08'};
RCM_rcp85_testnames = {'test65', 'test66', 'test67', 'test68', 'ens10'};
% RCM_rcp26_testnames = { 'ens09'};
% RCM_rcp45_testnames = {'ens08'};
% RCM_rcp85_testnames = { 'ens10'};
% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

all_region2 ={'AKP4'};
% all_region2 ={'NWP'};


% scenname='rcp26';
% scenname='rcp45';
% scenname='rcp85';

close all;
% clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 testname testname2 sodatestname
% % % 
system_name=computer;

dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

regionind2 = 1;

shadlev = [0 35];
rms_shadlev = [0 4];
%     trendlev = [-3 3];  %% trend lev
trendlev = [-10 10];  %% trend lev
abstrendlev =[2 7];
reltrendlev =[-5 5];
conlev  = 0:5:35;
meanplotlev =[-0.2 0.2];
meanplotlev2 =[-0.1 0.9] .*100.0;

trendplotlev = [0 7];
trenddifflev = [-10 10];
sstlev =[-0.3 0.3];
sstdifflev = [0 20];

% for snu_desktopd
% testname=all_testname2{testnameind2}    % % need to change

pastyear = [2006:2020];
inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

varname ='temp';
variable='sst';
run('nwp_polygon_point.m');
regionname=all_region2{regionind2};


for folding=1:1
    switch(regionname)
        case('NWP') %% North western Pacific
            lonlat = [115, 164, 15, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('NWP2') %% North western Pacific
            lonlat = [115, 145, 25, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('ES') %% East Sea
            refpolygon=espolygon;
        case('NES') %% Northern East Sea
            refpolygon=nespolygon;
        case('SES') %% Southern East Sea
            refpolygon=sespolygon;
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        case('AKP') %% Around Korea Peninsula
            refpolygon=akppolygon;
        case('AKP2') %% Around Korea Peninsula
            refpolygon=akp2polygon;
        case('AKP3') %% Around Korea Peninsula
            refpolygon=akp3polygon;
        case('AKP4') %% Around Korea Peninsula
            refpolygon=akp4polygon;
        case('CA') %% Around Korea Peninsula
            refpolygon=capolygon;
        case('EKB') %% Around Korea Peninsula
            refpolygon=akp2polygon;
        otherwise
            ('?')
    end
    lonlat(1)=min(refpolygon(:,1));
    lonlat(2)=max(refpolygon(:,1));
    lonlat(3)=min(refpolygon(:,2));
    lonlat(4)=max(refpolygon(:,2));
end



for testind=1:length(RCM_rcp26_testnames)
    testname=RCM_rcp26_testnames{testind};
    drivename='D';
    filedir=strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\');
    RCM_rcp26_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp26_cmemsfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp26_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp26_movfilenames{testind} = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%     cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
%     cmems_filename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
end

for testind=1:length(RCM_rcp45_testnames)
    testname=RCM_rcp45_testnames{testind};
    drivename='D';
    filedir=strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\');
    RCM_rcp45_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp45_cmemsfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp45_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp45_movfilenames{testind} = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
end

for testind=1:length(RCM_rcp85_testnames)
    testname=RCM_rcp85_testnames{testind};
    drivename='D';
    filedir=strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\');
    RCM_rcp85_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp85_cmemsfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp85_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_rcp85_movfilenames{testind} = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
end



nyears=[1,2,5];
for testind=1:length(RCM_rcp26_testnames)
    RCM_rcp26_interped_sla(testind,:,:,:) = ncread(RCM_rcp26_interpedfilenames{testind}, 'interped_sst');
    xlen=size(RCM_rcp26_interped_sla,2); ylen = size(RCM_rcp26_interped_sla,3); tlen = size(RCM_rcp26_interped_sla,4);

    RCM_rcp26_interped_sla_2d(testind,:,:)=reshape(RCM_rcp26_interped_sla(testind,:,:,:), [xlen*ylen, tlen]);
    RCM_rcp26_interped_sla_3d(testind,:,:)=reshape(RCM_rcp26_interped_sla(testind,:,:),[xlen*ylen*12, tlen/12]);
    RCM_rcp26_interped_sla_yearly_mean(testind,:)=mean(RCM_rcp26_interped_sla_3d(testind,:,:),2,'omitnan') ; % (m -> cm)

    RCM_rcp26_interped_sla_mean(testind,:) = mean(RCM_rcp26_interped_sla_2d(testind,:,:),2,'omitnan');

    RCM_rcp26_interped_sla_4d(testind,:,:,:)=reshape(RCM_rcp26_interped_sla(testind,:,:,:), [xlen*ylen, 12, tlen/12]);
    RCM_rcp26_interped_sla_seasonal_mean=mean(mean(RCM_rcp26_interped_sla_4d,2,'omitnan'),4,'omitnan');
    
    for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    RCM_rcp26_interped_sla_seasonal_filtered(testind,t)=RCM_rcp26_interped_sla_mean(testind,t)-RCM_rcp26_interped_sla_seasonal_mean(testind,tt);
    end

%     for nyeari=1:length(nyears)
%         nyear=nyears(nyeari);
%         nc_varname=['interped_', num2str(nyear), 'y_movmean'];
%         eval([nc_varname, 's(testind,:,:,:)=squeeze(ncread(RCM_rcp26_movfilenames{testind},', '''', nc_varname, '''', ',[1,1,1], [inf inf inf]));']);
%     end
end

for testind=1:length(RCM_rcp45_testnames)
    RCM_rcp45_interped_sla(testind,:,:,:) = ncread(RCM_rcp45_interpedfilenames{testind}, 'interped_sst');
    xlen=size(RCM_rcp45_interped_sla,2); ylen = size(RCM_rcp45_interped_sla,3); tlen = size(RCM_rcp45_interped_sla,4);

    RCM_rcp45_interped_sla_2d(testind,:,:)=reshape(RCM_rcp45_interped_sla(testind,:,:,:), [xlen*ylen, tlen]);
    RCM_rcp45_interped_sla_3d(testind,:,:)=reshape(RCM_rcp45_interped_sla(testind,:,:),[xlen*ylen*12, tlen/12]);
    RCM_rcp45_interped_sla_yearly_mean(testind,:)=mean(RCM_rcp45_interped_sla_3d(testind,:,:),2,'omitnan') 

    RCM_rcp45_interped_sla_mean(testind,:) = mean(RCM_rcp45_interped_sla_2d(testind,:,:),2,'omitnan');

    RCM_rcp45_interped_sla_4d(testind,:,:,:)=reshape(RCM_rcp45_interped_sla(testind,:,:,:), [xlen*ylen, 12, tlen/12]);
    RCM_rcp45_interped_sla_seasonal_mean=mean(mean(RCM_rcp45_interped_sla_4d,2,'omitnan'),4,'omitnan');
    
    for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    RCM_rcp45_interped_sla_seasonal_filtered(testind,t)=RCM_rcp45_interped_sla_mean(testind,t)-RCM_rcp45_interped_sla_seasonal_mean(testind,tt);
    end

%     for nyeari=1:length(nyears)
%         nyear=nyears(nyeari);
%         nc_varname=['interped_', num2str(nyear), 'y_movmean'];
%         eval([nc_varname, 's(testind,:,:,:)=squeeze(ncread(RCM_rcp45_movfilenames{testind},', '''', nc_varname, '''', ',[1,1,1], [inf inf inf]));']);
%     end
end

for testind=1:length(RCM_rcp85_testnames)
    RCM_rcp85_interped_sla(testind,:,:,:) = ncread(RCM_rcp85_interpedfilenames{testind}, 'interped_sst');
    xlen=size(RCM_rcp85_interped_sla,2); ylen = size(RCM_rcp85_interped_sla,3); tlen = size(RCM_rcp85_interped_sla,4);

    RCM_rcp85_interped_sla_2d(testind,:,:)=reshape(RCM_rcp85_interped_sla(testind,:,:,:), [xlen*ylen, tlen]);
    RCM_rcp85_interped_sla_3d(testind,:,:)=reshape(RCM_rcp85_interped_sla(testind,:,:),[xlen*ylen*12, tlen/12]);
    RCM_rcp85_interped_sla_yearly_mean(testind,:)=mean(RCM_rcp85_interped_sla_3d(testind,:,:),2,'omitnan') 

    RCM_rcp85_interped_sla_mean(testind,:) = mean(RCM_rcp85_interped_sla_2d(testind,:,:),2,'omitnan');

    RCM_rcp85_interped_sla_4d(testind,:,:,:)=reshape(RCM_rcp85_interped_sla(testind,:,:,:), [xlen*ylen, 12, tlen/12]);
    RCM_rcp85_interped_sla_seasonal_mean=mean(mean(RCM_rcp85_interped_sla_4d,2,'omitnan'),4,'omitnan');
    
    for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    RCM_rcp85_interped_sla_seasonal_filtered(testind,t)=RCM_rcp85_interped_sla_mean(testind,t)-RCM_rcp85_interped_sla_seasonal_mean(testind,tt);
    end

%     for nyeari=1:length(nyears)
%         nyear=nyears(nyeari);
%         nc_varname=['interped_', num2str(nyear), 'y_movmean'];
%         eval([nc_varname, 's(testind,:,:,:)=squeeze(ncread(RCM_rcp85_movfilenames{testind},', '''', nc_varname, '''', ',[1,1,1], [inf inf inf]));']);
%     end
end


valnum=0;
run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
wrmap = bwrmap(51:100,:);

if (strcmp(system_name,'PCWIN64'))
    % % for windows
    figrawdir =strcat('D:\MEPL\project\SSH\figures\5th_year\figure\all\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
    param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
%     filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
elseif (strcmp(system_name,'GLNXA64'))
end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

run(param_script);

figdir=[figrawdir,'sst\'];
if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 
outfile = strcat(figdir,regionname);

% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
            xData2(i) = datenum([num2str(tempyear),'-',num2str(6,'%02i'),'-30',]);
        end
% end-------------------- make timedata for time series  

% start-------------------- msl time series (yearly mean) (RCM)
for folding=1:1
        jpgname=strcat(outfile, '_all_scen_ens_',regionname, '_RCM_msst_yearly_mean_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            
        
           %% Hue, Saturation, Value
            %% If S goes 0, color get to be transparent
            
            val_transparent = 0.3;
            
            cmap_rcp26 = [0,0,1]; % blue
            cmap_rcp26_b = rgb2hsv(cmap_rcp26);
            cmap_rcp26_b(:,2) =  val_transparent;
            cmap_rcp26_b = hsv2rgb(cmap_rcp26_b);
            
            cmap_rcp45 = [1, 165/255, 0];
            cmap_rcp45_b = rgb2hsv(cmap_rcp45);
            cmap_rcp45_b(:,2) =  val_transparent;
            cmap_rcp45_b = hsv2rgb(cmap_rcp45_b);
            
            cmap_rcp85 = [1, 0, 0];
            cmap_rcp85_b = rgb2hsv(cmap_rcp85);
            cmap_rcp85_b(:,2) =  val_transparent;
            cmap_rcp85_b = hsv2rgb(cmap_rcp85_b);
            
            hold on
            cmems_1y_mean=0;   
            pastlen=length(pastyear);
            alllen=length(inputyear);
            
            for testind=1:length(RCM_rcp26_testnames)
                RCM_rcp26_interped_1y_mean(testind)=mean(RCM_rcp26_interped_sla_yearly_mean(testind,1:5));
                if testind==length(RCM_rcp26_testnames)
                    mslplot_rcp26_all{testind}=plot(xData2(1:pastlen),RCM_rcp26_interped_sla_yearly_mean(testind,1:pastlen),'Color', cmap_rcp26)
                    mslplot_rcp26_all{testind+1}=plot(xData2(pastlen:end),RCM_rcp26_interped_sla_yearly_mean(testind, pastlen:end),'Color', cmap_rcp26)
                else
                    mslplot_rcp26_all{testind}=plot(xData2,RCM_rcp26_interped_sla_yearly_mean(testind,:),'Color', cmap_rcp26)                    
                end
                if testind == length(RCM_rcp26_testnames)
                    set(mslplot_rcp26_all{testind}, 'Linewidth', 9)
                    set(mslplot_rcp26_all{testind+1}, 'Linewidth', 9)
                end
            end
            
            for testind=1:length(RCM_rcp45_testnames)
                RCM_rcp45_interped_1y_mean(testind)=mean(RCM_rcp45_interped_sla_yearly_mean(testind,1:5));
                if testind==length(RCM_rcp45_testnames)
                    mslplot_rcp45_all{testind}=plot(xData2(1:pastlen),RCM_rcp45_interped_sla_yearly_mean(testind,1:pastlen),'Color', cmap_rcp45)
                    mslplot_rcp45_all{testind+1}=plot(xData2(pastlen:end),RCM_rcp45_interped_sla_yearly_mean(testind, pastlen:end),'Color', cmap_rcp45)
                else
                    mslplot_rcp45_all{testind}=plot(xData2,RCM_rcp45_interped_sla_yearly_mean(testind,:),'Color', cmap_rcp45)
                end
                if testind == length(RCM_rcp45_testnames)
                    set(mslplot_rcp45_all{testind}, 'Linewidth', 9)
                    set(mslplot_rcp45_all{testind+1}, 'Linewidth', 9)
                end
            end
            
            for testind=1:length(RCM_rcp85_testnames)
                RCM_rcp85_interped_1y_mean(testind)=mean(RCM_rcp85_interped_sla_yearly_mean(testind,1:5));
                if testind==length(RCM_rcp85_testnames)
                    mslplot_rcp85_all{testind}=plot(xData2(1:pastlen),RCM_rcp85_interped_sla_yearly_mean(testind,1:pastlen),'Color', cmap_rcp85)
                    mslplot_rcp85_all{testind+1}=plot(xData2(pastlen:end),RCM_rcp85_interped_sla_yearly_mean(testind, pastlen:end),'Color', cmap_rcp85)
                else
                    mslplot_rcp85_all{testind}=plot(xData2,RCM_rcp85_interped_sla_yearly_mean(testind,:),'Color', cmap_rcp85)                    
                end
                if testind == length(RCM_rcp85_testnames)
                    set(mslplot_rcp85_all{testind}, 'Linewidth', 9)
                    set(mslplot_rcp85_all{testind+1}, 'Linewidth', 9)
                end
            end
            
            for testind=1:length(RCM_rcp26_testnames)
                set(mslplot_rcp26_all{testind}, 'Color', cmap_rcp26_b)
                set(mslplot_rcp45_all{testind}, 'Color', cmap_rcp45_b)
                set(mslplot_rcp85_all{testind}, 'Color', cmap_rcp85_b)
                set(get(get(mslplot_rcp26_all{testind},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                set(get(get(mslplot_rcp45_all{testind},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                set(get(get(mslplot_rcp85_all{testind},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
                                    
            iiii=1;
            for xtli=2010:10:2100
               xtlabel{iiii}=num2str(xtli);
               iiii=iiii+1;
            end
            set(gca, 'xtick', xData2(5:10:end),'xticklabel',xtlabel);
            
%             set(mslplot_all{1},'Marker','*');
%             set(mslplot_all{2},'Marker','^');
%             set(mslplot_rcp26_all{5},'Marker','o');
%             set(mslplot_rcp45_all{5},'Marker','o');
%             set(mslplot_rcp85_all{5},'Marker','o');

                        
%             set(mslplot_all{4},'Marker','+');
            
            xlabel('Year')
            ylabel('Mean SST (°C)')
%             title([regionname, ', Mean sst(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
%             ylim(meanplotlev2)
            set(gca,'FontSize',15);
            grid on
            lgd=legend('RCP 2.6','RCP 4.5','RCP 8.5');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
%             for tind=1:size(RCM_interped_sla_yearly_mean,2)
%                 model_std(tind)=std(RCM_interped_sla_yearly_mean(:,tind));
%             end

            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver 22 17]) 
            hold off
            saveas(gcf,jpgname,'jpg');

            grid off
%         end
        close all;
        
        README='regional model SST, Time(xData2). xData2 should be used with datetick func';
        savefilename='D:\Data\Model\ROMS\nwp_1_20\all\RCP_SST_yearly_data.mat';
        save(savefilename, 'README', 'xData2', 'RCM_rcp26_interped_sla_yearly_mean', ...
            'RCM_rcp45_interped_sla_yearly_mean', 'RCM_rcp85_interped_sla_yearly_mean');
        
end
% end-------------------- msl time series (yearly mean) (RCM)

