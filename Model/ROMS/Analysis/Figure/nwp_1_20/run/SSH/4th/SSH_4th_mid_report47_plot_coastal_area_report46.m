close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'test53','test54','test55','test56','ens03'};
% all_testname2 = {'test57', 'test58','test59','test60'};
all_testname2 = {'ens03'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP', 'AKP2'}

all_region2 ={'AKP2'};
% all_dist2 = [0.05, 0.1, 0.5];
all_dist2 = [0.1];

% all_region2 ={'NWP'}
for distind2= 1:length(all_dist2)
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_dist2 distind2 mask_coast mask_coast_all
        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIAS
        dist=all_dist2(distind2);
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\KYY\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end

        shadlev = [0 35];
        rms_shadlev = [0 4];
    %     trendlev = [-3 3];  %% trend lev
        trendlev = [-10 10];  %% trend lev
        abstrendlev =[4 7];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.3 0.3];
        trendplotlev = [4 6.5];
        sshlev =[-0.2 1.0];
        sshdifflev = [30 70];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1976:2005]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

        varname ='zeta';
        var='SSH';
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
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

%         load(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        filename = ['E:\Data\Model\ROMS\nwp_1_20\',testname,'\',testname,'_',regionname, ...
                    '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);

    % %     valid cell number
    %      for vi=1:size(comb_spatial_meanressh,1)
    %          for vj=1:size(comb_spatial_meanressh,2)
    %              if (isnan(comb_spatial_meanressh(vi,vj,1))~=1)
    %                 valnum=valnum+1;
    %              end
    %          end
    %      end

    %     plot(squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan')))
    %     isize = size(comb_data_filtered,1)
    %     jsize = size(comb_data_filtered,2)
    %     lsize = size(comb_data_filtered,3)
    %     comb_yearly_data_filtered=reshape(comb_data_filtered,[isize, jsize, 12, lsize/12]);
    %     mean_yearly_data_filtered=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
    %     trendtime=14:29;
    %     p=polyfit(trendtime,mean_yearly_data_filtered(14:29)',1);
    %     yearly_interped_trend=p(1);
    %     yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y


        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\', testname, '\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))
        
        run(param_script);
        
        figdir=[figrawdir,'Trend\SSH\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        
        
        lon_rho=ncread(filename,'lon_rho');
        lat_rho=ncread(filename,'lat_rho');
        trend_filtered = ncread(filename,'trend_filtered');
        mean_trend_filtered = ncread(filename,'mean_trend_filtered');
    

%         comb_data=ncread(filename, 'raw_ssh');   
%         comb_spatial_data=reshape(comb_data, [size(lon_rho,1) size(lat_rho,2) 12 length(inputyear)]);
%         climtrendtime=inputyear(1):inputyear(end);
%         for i=1:size(lon_rho,1)
%             for j=1:size(lat_rho,2)
%                 for k=1:12  % month
%                     comb_spatial_data(i,j,k,:)=comb_spatial_data(i,j,k,:)-mean(comb_spatial_data(i,j,k,:),'omitnan');
%                     p=polyfit(climtrendtime,squeeze(comb_spatial_data(i,j,k,:))',1);
%                     clim_trend(i,j,k)=p(1);
%                 end
%             end
%         end

        clim_trend = ncread(filename,'clim_ssh_trend');
        mean_clim_trend = squeeze(mean(mean(clim_trend, 1, 'omitnan'), 2, 'omitnan'));
        
        colval_west=1.7;
        colval_south = 2.5;
        colval_east = 4.0;
        if exist('mask_coast') ==0
            load('C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Roms_tools\coastline_i.mat')
            numx=size(lon_rho,1);
            numy=size(lat_rho,2);
            for i=1:numx
                for j=1:numy
                  mask_rho(i,j)=double(isfinite(trend_filtered(i,j)));
                end
            end
            mask_coast=mask_rho;
%             dist=0.1;
            % for ncsti=1:size(ncst,1)
            %     mask_coast(find(lon_rho>= ncst(ncsti,1)-dist & lon_rho<=ncst(ncsti,1)+dist ...
            %         & lat_rho<=ncst(ncsti,2)+dist & lat_rho>=ncst(ncsti,2)-dist))=2;
            % end
            % mask_coast(mask_rho==0)=0;
            for i=1:numx
                for j=1:numy
                    if mask_rho(i,j)==0
                        mask_coast(find(lon_rho>= lon_rho(i,j)-dist & lon_rho<=lon_rho(i,j)+dist ...
                    & lat_rho<=lat_rho(i,j)+dist & lat_rho>=lat_rho(i,j)-dist))=2;
                    end
                end
            end
            mask_coast(mask_rho==0)=0;
            
            temp_mask_coast_west = inpolygon(lon_rho,lat_rho,yspolygon2(:,1), yspolygon2(:,2));
            mask_coast_west=mask_coast .* temp_mask_coast_west;
    %         pcolor(mask_coast_west')
    %         shading flat;

            temp_mask_coast_south = inpolygon(lon_rho,lat_rho,sspolygon2(:,1), sspolygon2(:,2));
            mask_coast_south=mask_coast .* temp_mask_coast_south;
    %         pcolor(mask_coast_south')
    %         shading flat;

            temp_mask_coast_east = inpolygon(lon_rho,lat_rho,espolygon2(:,1), espolygon2(:,2));
            mask_coast_east=mask_coast .* temp_mask_coast_east;
            pcolor(mask_coast_east')
            shading flat;

            mask_coast_all=mask_rho;

            mask_coast_all(mask_coast_west==2)=colval_west;
            mask_coast_all(mask_coast_south==2)=colval_south;
            mask_coast_all(mask_coast_east==2)=colval_east;
    %         mask_coast_all(mask_coast_all==1)=NaN;
            pcolor(lon_rho(160:330, 40:200)', lat_rho(160:330, 40:200)',mask_coast_all(160:330, 40:200)')
            shading flat;
            grid on
            colormap(jet)
%             colorbar
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_mask_coast_all_dist_',num2str(dist*100), 'km.tif'); %% ~_year_month.jpg
            saveas(gcf,jpgname,'jpg');
            close all;
        end
        
        % save('mask_coast_all.mat', 'mask_rho', 'mask_coast_all')

        trend_filtered_west=trend_filtered(mask_coast_all==colval_west);
        m_west=mean(trend_filtered_west(:));
        disp([testname, ', Yellow Sea trend : ', num2str(round(m_west,2)), ' mm/year, dist : ', num2str(dist), 'deg'])
        trend_filtered_south=trend_filtered(mask_coast_all==colval_south);
        m_south=mean(trend_filtered_south(:));
        disp([testname, ', South Sea trend : ', num2str(round(m_south,2)), ' mm/year, dist : ', num2str(dist), 'deg'])
        trend_filtered_east=trend_filtered(mask_coast_all==colval_east);
        m_east=mean(trend_filtered_east(:));
        disp([testname, ', East Sea trend : ', num2str(round(m_east,2)), ' mm/year, dist : ', num2str(dist), 'deg'])
        trend_filtered_coast= [trend_filtered_west(:)', trend_filtered_south(:)', trend_filtered_east(:)'];
        m_coast=mean(trend_filtered_coast(:));
        disp([testname, ', Korean Coastal Area trend : ', num2str(round(m_coast,2)), ' mm/year, dist : ', num2str(dist), 'deg'])

% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
% end-------------------- make timedata for time series  

% % start-------------------- msl time series (seasonal filtered)
%         jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_', ...
%             num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             for varind=1:length(inputyear)*12
%                 msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%             end
%     %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
%             msl_filt=msl_filt-mean(msl_filt);    
%             p=polyfit(xData,msl_filt,1);
%             msl2=xData*p(1)+p(2);
%             mslplot=plot(xData,msl_filt,'k')
%             hold on
%             mslplot2=plot(xData,msl2,'Color','r')
%             xlabel('year')
%             ylabel('Mean SSH (m)')
%             title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend_filtered,2)), ' mm/y'])
%             datetick('x','yyyy','keepticks')
%             axis tight;
%             ylim(meanplotlev)
%             set(mslplot,'LineWidth',2);
%             set(gca,'FontSize',20);
%             grid on
%             hold off
%             saveas(gcf,jpgname,'jpg');
%             grid off
%         end
%         close all;
% % end-------------------- msl time series (seasonal filtered)
% 
% % start-------------------- msl time series
%         jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             for varind=1:length(inputyear)*12
%                 msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%             end
%             msl=msl-mean(msl);    
%             p=polyfit(xData,msl,1);
%             msl2=xData*p(1)+p(2);
%             mslplot=plot(xData,msl,'k')
%             hold on
%             mslplot2=plot(xData,msl2,'Color','r')
%             xlabel('year')
%             ylabel('Mean SSH (m)')
%             mean_trend=ncread(filename,'mean_trend');
%             title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend,2)), ' mm/y'])
%             datetick('x','yyyy','keepticks')
%             axis tight;
%             ylim(meanplotlev)
%             set(mslplot,'LineWidth',2);
%             set(gca,'FontSize',20);
%             grid on
%             hold off
%             saveas(gcf,jpgname,'jpg');
%             grid off
%             close all;
%         end
% % end-------------------- msl time series
% 
% % start-------------------- climatological msl time series
%         climdir = [figdir,'\CLIM\'];
%         if (exist(strcat(climdir) , 'dir') ~= 7)
%             mkdir(strcat(climdir));
%         end 
%         climoutfile = strcat(climdir,regionname);
%         for monthij=1:12
%             jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_', ...
%                 num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
%             if (exist(jpgname , 'file') ~= 2)
%                 if (exist('msl')==0)
%                     for varind=1:length(inputyear)*12
%                         msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%                     end
%                     climmsl=reshape(msl,[12,length(inputyear)]);
%                 else
%                     climmsl=reshape(msl,[12,length(inputyear)]);
%                 end
%                 tempmsl=squeeze(climmsl(monthij,:));
%                 tempmsl=tempmsl-mean(tempmsl);
%                 p=polyfit(inputyear,tempmsl,1);
%                 msl2=inputyear*p(1)+p(2);
%                 mslplot=plot(inputyear,tempmsl,'k')
%                 hold on
%                 mslplot2=plot(inputyear,msl2,'Color','r')
%                 xlabel('year')
%                 ylabel('Mean SSH (m)')
%                 title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
%                     ', ', calendarname{monthij}(1:3), '), ',num2str(round(mean_clim_trend(monthij),2)), ' mm/y'])
% %                 datetick('x','yyyy','keepticks')
%                 axis tight;
%                 ylim(meanplotlev)
%                 set(mslplot,'LineWidth',2);
%                 set(gca,'FontSize',20);
%                 grid on
%                 hold off
%                 saveas(gcf,jpgname,'jpg');
%                 grid off
%                 close all;
%             end
%         end
% % end-------------------- climatological msl time series

% % start-------------------- climatological msl trend
%         climdir = [figdir,'\CLIM\'];
%         if (exist(strcat(climdir) , 'dir') ~= 7)
%             mkdir(strcat(climdir));
%         end 
%         climoutfile = strcat(climdir,regionname);
%         jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_trend_', ...
%             num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
%         if abs(mean_clim_trend(1))<0.01
%             mean_clim_trend=mean_clim_trend*1000.0  %% m -> mm
%         end
%         if (exist(jpgname , 'file') ~= 2)
%             mslplot=plot(1:12,mean_clim_trend,'k')
%             hold on
%             xlabel('month')
%             ylabel('trend (mm/yr)')
%             title([regionname, ', climatological trend(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
%             axis tight;
%             ylim(trendplotlev)
%             set(mslplot,'LineWidth',2);
%             set(gca,'FontSize',20);
%             grid on
%             hold off
%             saveas(gcf,jpgname,'jpg');
%             grid off
%             close all;
%         end
% % end-------------------- climatological msl trend    
        

    end
end
    clear mask_coast
end