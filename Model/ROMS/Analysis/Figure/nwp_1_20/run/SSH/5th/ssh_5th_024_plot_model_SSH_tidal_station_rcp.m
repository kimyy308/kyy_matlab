close all; clear all;  clc;
warning off;

all_testname2 = {'ens10'};
% all_testname2 = {'test53'};

% all_region2 ={'AKP4'}

% all_region2 ={'YS'};

% all_var2 = {'SST', 'SSS', 'SSH', 'BT'};
all_var2 = {'SSH'};

% all_region2 ={'NWP'}
for testnameind2=1:length(all_testname2)
        close all;
%         clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2
        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIAS
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\user\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end

%         shadlev = [0 35];
%         rms_shadlev = [0 4];
%     %     trendlev = [-3 3];  %% trend lev
%         trendlev = [-10 10];  %% trend lev
%         abstrendlev =[4 7];
%         reltrendlev =[-5 5];
%         conlev  = 0:5:35;
        meanplotlev =[-0.5 0.5];
        trendplotlev = [6 14];
%         sshlev =[-0.7 1.3];
%         sshdifflev = [40 70];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

%         varname ='zeta'
        run('nwp_polygon_point.m');
%         regionname=all_region2{regionind2};
        regionname='AKP4';

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

        % % % for EKB
        % regionname='EKB';
        % lonlat = [127, 129.5, 38, 40.5];

%         load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
%         filename = ['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\user\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        run(param_script);

        
        
        filename=[filedir, testname, '_', regionname, '_ssh_trend_', num2str(min(inputyear)), '_', num2str(max(inputyear)), '.nc'];
        trendmatname=[filedir, testname, '_', regionname, '_SSH_trend_sta_', num2str(min(inputyear)), '-', num2str(max(inputyear)), '.mat'];
        
        tidal_st= [ 126.5919, 37.4528; ...  % Incheon (Yellow Sea)
            126.1322, 36.6436; ...  % Anheung
            126.4858, 36.4065; ...  % Boryung
            126.5431, 35.9756; ...  % Gunsan
            126.3035, 35.6177; ...  % Wuido
            126.3756, 34.7297; ...  % Mokpo
            125.4356, 34.6842; ...  % Heuksando
            126.3003, 33.9619; ...  % Chujado (South Sea)
            126.7597, 34.2756; ...  % Wando
            127.8056, 34.7472; ...  % Yeosu
            128.4547, 34.8278; ...  % Tongyeong
            128.8008, 35.0242; ...  % Gadeokdo
            129.0453, 35.0564; ...  % Busan
            126.5431, 33.5675; ...  % Jeju
            126.5617, 33.2400; ...  % Seoguipo
            127.3089, 33.9983; ...  % Geomundo
            129.3872, 35.4719; ...  % Ulsan (East Sea)
            129.3839, 36.0372; ...  % Pohang
            129.1164, 37.6003; ...  % Mukho
            128.5942, 38.2072; ...  % Sokcho
            130.9136, 37.5614];  % Ulleungdo
        
        tidal_name = { 'West-01-Incheon', 'West-02-Anheung', 'West-03-Boryung', 'West-04-Gunsan', 'West-05-Wuido', 'West-06-Mokpo', 'West-07-Heuksando', ...
                       'South-08-Chujado', 'South-09-Wando', 'South-10-Yeosu', 'South-11-Tongyeong', 'South-12-Gadeokdo',  ...
                       'South-13-Busan', 'South-14-Jeju', 'South-15-Seoguipo', 'South-16-Geomundo' ...
                       'East-17-Ulsan', 'East-18-Pohang', 'East-19-Mukho', 'East-20-Sokcho', 'East-21-Ulleungdo'};
        info=ncinfo(filename);
        lon_rho=ncread(filename, 'lon_rho');
        lat_rho=ncread(filename, 'lat_rho');
        mask_rho=ncread(filename, 'trend_filtered');
        mask_rho(~isnan(mask_rho))=1;
        lon_rho=lon_rho.*mask_rho;
        lat_rho=lat_rho.*mask_rho;
%         raw_ssh=ncread(filename, 'raw_ssh');
        
%         for ind_sta=1:size(tidal_st,1)
%             for i=1:size(lon_rho,1)
%                 for j=1:size(lat_rho,2)
%                     dist(i,j)=m_lldist([lon_rho(i,j), tidal_st(ind_sta,1)], [lat_rho(i,j), tidal_st(ind_sta,2)]);
%                 end
%             end
%             ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
%             ind_sta_model_lon(ind_sta)=mod(ind_sta_model(ind_sta),size(lon_rho,1));
%             ind_sta_model_lat(ind_sta)=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
% %             model_msl(:,ind_sta)=raw_ssh(ind_sta_model_lon(ind_sta),ind_sta_model_lat(ind_sta),:);
%             model_msl(:,ind_sta)=ncread(filename, 'raw_ssh', [ind_sta_model_lon(ind_sta), ind_sta_model_lat(ind_sta), 1], [1, 1, inf]);
%         end
% %         lon_rho(ind_sta_model(1))
% %         lat_rho(ind_sta_model(1))
% %         lon_rho(ind_sta_model_lon(1), ind_sta_model_lat(1))
% %         lat_rho(ind_sta_model_lon(1), ind_sta_model_lat(1))
%         clim_model_msl=reshape(model_msl,[12, size(model_msl,1)/12, size(tidal_st,1)]);
%         clim_mean_msl=squeeze(mean(clim_model_msl,2));
%         
%         trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
%         trend(1:ind_sta)=NaN;
%         for ind_sta=1:size(tidal_st,1)
%             for ind_time=1:size(model_msl,1)
% %                 mod(ind_time-1,12)+1
%                 clim_filtered_msl(ind_time,ind_sta)=model_msl(ind_time,ind_sta)-clim_mean_msl(mod(ind_time-1,12)+1,ind_sta);
%             end
%             p=polyfit(trendtime,squeeze(clim_filtered_msl(:,ind_sta))',1);
%             trend(ind_sta)=p(1) * 1000.0;
%         end
%         plot(squeeze(clim_filtered_msl(:,1))*100)
        
        figdir=[figrawdir,'Tidal\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
   
    % start-------------------- plot trend along station
    load(trendmatname);
        outfile = strcat(figdir,'Trends');
        jpgname=strcat(outfile, '_', testname, '_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
            trendplot=plot(trend,'-o','MarkerFaceColor','k')
            xticks(1:21)
            xticklabels(tidal_name)
            xtickangle(45)
            xlabel('Station')
            ylabel('Trend (mm/year)')
            title([testname, ' trends (',num2str(min(inputyear),'%04i'), ...
                    '-',num2str(max(inputyear),'%04i'),'), ','M = ',num2str(round(mean(trend),2)), ' mm/y'])
%             axis tight;
            ylim([7 9])
            set(gcf,'PaperPosition', [0 0 36 12]) 
            set(trendplot,'LineWidth',2);
            set(gca,'FontSize',20);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all
%         end

end