close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
all_testname2 = {'test11', 'test12'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP4', 'YS', 'YSECS', 'ECS2', 'NES', 'SES', 'ES'}

% all_region2 ={'BOH'};

all_region2 ={'AKP4'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2
        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIAS
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
        abstrendlev =[2 7];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.15 0.15];
        trendplotlev = [0 7];
        trenddifflev = [-10 10];
        sshlev =[-0.3 0.3];
        sshdifflev = [0 20];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1994:2014]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]

        varname ='zeta'
        variable='SSH'
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
             case('ECS2') %% East China Sea
                refpolygon=ecs2polygon;
            case('YSECS') %% East China Sea
                refpolygon=ysecspolygon;
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
            case('BOH') %% Around Korea Peninsula
                refpolygon=bohpolygon;
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

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);
        
        cmems_trend=ncread(filename, 'cmems_trend');
        cmems_mask=ones(size(cmems_trend));
        cmems_mask(isnan(cmems_trend))=NaN;
            
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
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_10\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_', regionname, '.m']
            filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
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
        cmems_trend_filtered = ncread(filename,'cmems_trend_filtered');
        interped_trend_filtered = ncread(filename,'interped_trend_filtered');

        
% start-------------------- corr_interped plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_interped_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            corr_interped=ncread(filename,'corr_interped');
            corr_interped=corr_interped.*cmems_mask;
            lon_cmems=ncread(filename,'lon_cmems');
            lat_cmems=ncread(filename,'lat_cmems');
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(lon_cmems',lat_cmems',squeeze(corr_interped(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('corr, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped(:), 'omitnan'),2)));  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(wrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([0.2 0.8]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- corr_interped plot plot

% start-------------------- corr_interped_filtered plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_interped_filtered_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            corr_interped_filtered=ncread(filename,'corr_interped_filtered');
            lon_cmems=ncread(filename,'lon_cmems');
            lat_cmems=ncread(filename,'lat_cmems');
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(lon_cmems',lat_cmems',squeeze(corr_interped_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('corr clim filtered, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped_filtered(:), 'omitnan'),2)));  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(wrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([0.2 0.8]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- corr_interped_filtered plot plot


% start-------------------- corr_spatial mean plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_spatial_mean_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            corr_spatial_mean=ncread(filename,'corr_spatial_mean');
            lon_cmems=ncread(filename,'lon_cmems');
            lat_cmems=ncread(filename,'lat_cmems');
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(lon_cmems',lat_cmems',squeeze(corr_spatial_mean(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('corr climatology, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_spatial_mean(:), 'omitnan'),2)));  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(wrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([0.2 0.8]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- corr_spatial mean plot



% start-------------------- corr_interped_detrended plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_interped_detrended_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            corr_interped_detrended=ncread(filename,'corr_interped_detrended');
            lon_cmems=ncread(filename,'lon_cmems');
            lat_cmems=ncread(filename,'lat_cmems');
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(lon_cmems',lat_cmems',squeeze(corr_interped_detrended(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('corr filt detrend, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped_detrended(:), 'omitnan'),2)));  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(wrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([0.2 0.8]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- corr_interped_detrended plot plot


% start-------------------- corr_ lowpassed plot
        nyears=[2,3,5];
        nc_varname_prefixes={'interped_', 'interped_detrended_'};
        nc_titlename_prefixes={'lowpass', 'lowpass-det'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
                    eval(['corr_',nc_varname, '=ncread(filename,','''','corr_', nc_varname,'''',');']);
                    eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);

                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                    tit_prefix=nc_titlename_prefixes{nc_varnameij};
                    titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                        num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr,2)));  

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(wrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    caxis([0.2 0.8]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['corr_', nc_varname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
            end
        end
% end-------------------- corr_ lowpassed plot

% start-------------------- corr_ lowpassed plot (corrected)
        nyears=[2,5];
        nc_varname_prefixes={'corrected_interped_sla_'};
        nc_titlename_prefixes={'lowpass-c'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
                    eval(['corr_',nc_varname, '=ncread(filename,','''','corr_', nc_varname,'''',');']);
                    eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);

                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                    tit_prefix=nc_titlename_prefixes{nc_varnameij};
                    titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                        num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr,2)));  

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(wrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    caxis([0.2 0.8]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['corr_', nc_varname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
            end
        end
% end-------------------- corr_ lowpassed plot (corrected)

% start-------------------- corr_ lowpassed plot (corrected_diff)
        nyears=[2,5];
        nc_varname_prefixes={'corrected_interped_sla_'};
        nc_varname_prefixes2={'interped_'};
        nc_titlename_prefixes={'lowpass-c'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', 'diff','_',regionname, '_corr_diff_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
                    diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                    diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                    eval(['corr_',nc_varname, '=ncread(diffname,','''','corr_', nc_varname2,'''',');']);
                    eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);
                    eval(['corr2_',nc_varname, '=ncread(diffname2,','''','corr_', nc_varname,'''',');']);
                    eval(['corr2_',nc_varname, '=corr2_', nc_varname, '.*cmems_mask;']);
                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''','-','corr2_', nc_varname, '(:,:)', '''',  '));']);
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                    eval(['mcorr2=mean(corr2_',nc_varname,'(:),', '''', 'omitnan', '''',');']);

                    tit_prefix=nc_titlename_prefixes{nc_varnameij};
                    titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                        num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr-mcorr2,2)));  

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    caxis([-0.5 0.5]);
% corr_corrected_interped_sla_2y_lowpass-corr2_corrected_interped_sla_2y_lowpass
                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['corr_', nc_varname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
            end
        end
% end-------------------- corr_ lowpassed plot (corrected_diff)

% start-------------------- corr_ lowpassed plot (corrected_diff, normalized)
        nyears=[2,5];
        nc_varname_prefixes={'corrected_interped_sla_'};
        nc_varname_prefixes2={'interped_'};
        nc_titlename_prefixes={'lowpass-c'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', 'diff','_',regionname, '_corr_diff_norm_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%                 if (exist(jpgname , 'file') ~= 2)
                    diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                    diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                    eval(['corr_',nc_varname, '=ncread(diffname,','''','corr_', nc_varname2,'''',');']);
                    eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);
                    eval(['corr2_',nc_varname, '=ncread(diffname2,','''','corr_', nc_varname,'''',');']);
                    eval(['corr2_',nc_varname, '=corr2_', nc_varname, '.*cmems_mask;']);
                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    eval(['mean_corr_diff=', 'squeeze(corr_', nc_varname, '(:,:)', '-', 'corr2_', nc_varname, '(:,:));'])
                    std_mean_corr_diff=std(mean_corr_diff(:), 'omitnan');
                    n_mean_corr_diff=mean_corr_diff/std_mean_corr_diff;
                    eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',n_mean_corr_diff','''',  ');']);

%                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''','-','corr2_', nc_varname, '(:,:)', '''',  '));']);
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                    eval(['mcorr2=mean(corr2_',nc_varname,'(:),', '''', 'omitnan', '''',');']);

                    tit_prefix=nc_titlename_prefixes{nc_varnameij};
                    titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                        num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(n_mean_corr_diff(:), 'omitnan'),2)));  

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    caxis([-6 6]);
% corr_corrected_interped_sla_2y_lowpass-corr2_corrected_interped_sla_2y_lowpass
                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['corr_', nc_varname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
%                 end
            end
        end
% end-------------------- corr_ lowpassed plot (corrected_diff, normalized)


% start-------------------- diff mask between lowpassed corr and rms 
        nyears=[2,5];
        nc_varname_prefixes={'corrected_interped_sla_'};
        nc_varname_prefixes2={'interped_'};
        nc_titlename_prefixes={'lowpass-c'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', 'diff','_',regionname, '_rms_corr_diff_diff_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)


                    matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_diff_',nc_varname, '_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    load(matname)
                    rms_mask=NaN(size(mean_rms_diff));
                    rms_mask(mean_rms_diff<0)=-1;
                    rms_mask(mean_rms_diff>0)=1;
                    diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                    diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                    eval(['corr_',nc_varname, '=ncread(diffname,','''','corr_', nc_varname2,'''',');']);
                    eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);
                    eval(['corr2_',nc_varname, '=ncread(diffname2,','''','corr_', nc_varname,'''',');']);
                    eval(['corr2_',nc_varname, '=corr2_', nc_varname, '.*cmems_mask;']);
                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    eval(['mean_corr_diff=', 'squeeze(corr_', nc_varname, '(:,:)', '-', 'corr2_', nc_varname, '(:,:));'])
                    
                    corr_mask=NaN(size(mean_corr_diff));
                    corr_mask(mean_corr_diff<0)=-1;
                    corr_mask(mean_corr_diff>0)=1;
                    
                    diff_mask = rms_mask+corr_mask;
%                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',mean_corr_diff','''',  ');']);
                    m_pcolor(lon_cmems', lat_cmems', diff_mask');
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                    eval(['mcorr2=mean(corr2_',nc_varname,'(:),', '''', 'omitnan', '''',');']);

                    tit_prefix=nc_titlename_prefixes{nc_varnameij};
                    titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                        num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr-mcorr2,2)));  

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    caxis([-2 2]);
% corr_corrected_interped_sla_2y_lowpass-corr2_corrected_interped_sla_2y_lowpass
                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['corr_', nc_varname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
            end
        end
% end-------------------- diff mask between lowpassed corr and rms 

% start-------------------- diff between lowpassed corr and rms 
        nyears=[2,5];
        nc_varname_prefixes={'corrected_interped_sla_'};
        nc_varname_prefixes2={'interped_'};
        nc_titlename_prefixes={'diff-val-norm'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', 'diff','_',regionname, '_rms_corr_diff_norm_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)


                    matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_diff_',nc_varname, '_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    load(matname)

                    std_mean_rms_diff=std(mean_rms_diff(:),'omitnan');
                    n_mean_rms_diff=mean_rms_diff/std_mean_rms_diff;
                    
                    diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                    diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                    eval(['corr_',nc_varname, '=ncread(diffname,','''','corr_', nc_varname2,'''',');']);
                    eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);
                    eval(['corr2_',nc_varname, '=ncread(diffname2,','''','corr_', nc_varname,'''',');']);
                    eval(['corr2_',nc_varname, '=corr2_', nc_varname, '.*cmems_mask;']);
                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    eval(['mean_corr_diff=', 'squeeze(corr_', nc_varname, '(:,:)', '-', 'corr2_', nc_varname, '(:,:));'])
                    
                    std_mean_corr_diff=std(mean_corr_diff(:),'omitnan');
                    n_mean_corr_diff=mean_corr_diff/std_mean_corr_diff;
                    
%                     diff_mask = rms_mask+corr_mask;
                    diff_value = n_mean_rms_diff + n_mean_corr_diff;
%                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',mean_corr_diff','''',  ');']);
                    m_pcolor(lon_cmems', lat_cmems', diff_value');
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                    eval(['mcorr2=mean(corr2_',nc_varname,'(:),', '''', 'omitnan', '''',');']);

                    tit_prefix=nc_titlename_prefixes{nc_varnameij};
                    titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                        num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(diff_value(:),'omitnan'),2)));  

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    caxis([-4 4]);
% corr_corrected_interped_sla_2y_lowpass-corr2_corrected_interped_sla_2y_lowpass
                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['corr_', nc_varname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
            end
        end
% end-------------------- diff between lowpassed corr and rms 


% start-------------------- corr_clim plot
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        for monthij = 1:12
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_corr_clim_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
                corr_clim=ncread(filename,'corr_clim');
                lon_cmems=ncread(filename,'lon_cmems');
                lat_cmems=ncread(filename,'lat_cmems');

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',squeeze(corr_clim(:,:,monthij))');
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                
                tempcorr= corr_clim(:,:,monthij);
                titlename = strcat('corr clim, ', num2str(monthij, '%02i'), ', ', testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(tempcorr(:), 'omitnan'),2)));  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([0.2 0.8]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
        end
% end-------------------- corr_clim plot plot

        clim_trend = ncread(filename,'clim_ssh_trend');
        mean_clim_trend = squeeze(mean(mean(clim_trend, 1, 'omitnan'), 2, 'omitnan'));


% start-------------------- absolute trend plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(wrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis(abstrendlev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'absolute_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- absolute trend plot
 
% start-------------------- absolute trend plot (bwrmap)
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_bwrmap_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(bwrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([-4 4]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'absolute_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- absolute trend plot (bwrmap)


% start-------------------- relative trend plot (bwrmap)
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_relative_ssh_trend_bwrmap_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'-mean_trend_filtered));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(bwrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([-4 4]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'absolute_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- relative trend plot (bwrmap)



%  start-------------------- absolute trend difference plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_diff_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            lon_cmems=ncread(filename,'lon_cmems');
            lat_cmems=ncread(filename,'lat_cmems');
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            trend_diff = interped_trend_filtered-cmems_trend_filtered;
            m_pcolor(double(lon_cmems)',lat_cmems',trend_diff(:,:)');
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(trend_diff(:), 'omitnan'),2)), ' mm/y');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(bwrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis(trenddifflev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'absolute_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- absolute trend difference plot


% start-------------------- absolute trend_corrected plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_corrected_ssh_trend_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            correction_trend=mean(cmems_trend_filtered(:), 'omitnan')-mean_trend_filtered;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)')+correction_trend);
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered+correction_trend,2)), ' mm/y');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(wrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([0 7]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'absolute_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- absolute trend_corrected plot


% start-------------------- climatological absolute trend plot
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        for monthij = 1:12
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(lon_rho)',lat_rho',squeeze(clim_trend(:,:,monthij)'.*1000));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                    ', ', calendarname{monthij}(1:3), '), ','M=',num2str(round(mean_clim_trend(monthij)*1000,2)),'mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
%                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis(abstrendlev);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', num2str(monthij), '_', 'relative_ssh_trend', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
        end
        close all;
% end-------------------- climatological absolute trend plot
% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
% end-------------------- make timedata for time series  

% start-------------------- msl time series (seasonal filtered)
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            for varind=1:length(inputyear)*12
                msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
            end
    %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
            msl_filt=msl_filt-mean(msl_filt);    
            p2=polyfit(1:length(msl_filt),msl_filt,1);
            p2=p2(1)*1000.0*12.0;
            
            p=polyfit(xData,msl_filt,1);
            msl2=xData*p(1)+p(2);
            mslplot=plot(xData,msl_filt,'k')
            hold on
            mslplot2=plot(xData,msl2,'Color','r')
            xlabel('year')
            ylabel('Mean SSH (m)')
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
        end
        close all;
% end-------------------- msl time series (seasonal filtered)




% start-------------------- msl time series
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            for varind=1:length(inputyear)*12
                msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
            end
            msl=msl-mean(msl);     
            p2=polyfit(1:length(msl),msl,1);
            p2=p2(1)*1000.0*12.0;
            
            p=polyfit(xData,msl,1);
            msl2=xData*p(1)+p(2);
            mslplot=plot(xData,msl,'k')
            hold on
            mslplot2=plot(xData,msl2,'Color','r')
            xlabel('year')
            ylabel('Mean SSH (m)')
            mean_trend=ncread(filename,'mean_trend');
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
% end-------------------- msl time series

% start-------------------- msl time series (seasonal filtered) (corrected, SROCC)
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_corrected_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            for varind=1:length(inputyear)*12
                msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
            end
    %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
            
            for tind=1:length(msl_filt)
                msl_filt(tind)=msl_filt(tind)+(0.56 + 0.46 + 0.29 + 0.09)/1000.0/12.0*(tind-1);
            end
            msl_filt=msl_filt-mean(msl_filt);    
            p2=polyfit(1:length(msl_filt),msl_filt,1);
            p2=p2(1)*1000.0*12.0;
            
            
            p=polyfit(xData,msl_filt,1);
            msl2=xData*p(1)+p(2);
            mslplot=plot(xData,msl_filt,'k')
            hold on
            mslplot2=plot(xData,msl2,'Color','r')
            xlabel('year')
            ylabel('Mean SSH (m)')
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
        end
        close all;
% end-------------------- msl time series (seasonal filtered) (corrected, SROCC)

% start-------------------- msl time series (corrected, SROCC)
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_corrected_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            for varind=1:length(inputyear)*12
                msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
            end
            for tind=1:length(msl)
                msl(tind)=msl(tind)+(0.56 + 0.46 + 0.29 + 0.09)/1000.0/12.0*(tind-1);
            end
            msl=msl-mean(msl);     
            p2=polyfit(1:length(msl),msl,1);
            p2=p2(1)*1000.0*12.0;
            
            p=polyfit(xData,msl,1);
            msl2=xData*p(1)+p(2);
            mslplot=plot(xData,msl,'k')
            hold on
            mslplot2=plot(xData,msl2,'Color','r')
            xlabel('year')
            ylabel('Mean SSH (m)')
            mean_trend=ncread(filename,'mean_trend');
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
% end-------------------- msl time series (corrected, SROCC)

% start-------------------- lowpassed msl time series
        nyears=[2,3,5];
        nc_varname_prefixes={'cmems_', 'interped_'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', testname, '_',regionname, '_',nc_varname,'_msl_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
                    for varind=1:length(inputyear)*12
                        eval(['tempmsl=squeeze(ncread(filename,', '''', nc_varname, '''', ',[1,1,varind], [inf inf 1]));']);
                        msl(varind)=mean(tempmsl(:),'omitnan');
                        if (strcmp(nc_varname_prefix, 'interped_')==1)
                            temp_raw_msl=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
                        elseif (strcmp(nc_varname_prefix, 'cmems_')==1)
                            temp_raw_msl=ncread(filename,'cmems_sla',[1 1 varind], [inf inf 1]);
                        end
                        raw_msl(varind)=squeeze(mean(temp_raw_msl(:),'omitnan'));
                    end
                    msl=msl-mean(msl, 'omitnan');    
                    raw_msl=raw_msl-mean(raw_msl(:));

                    p=polyfit(xData,msl,1);
                    msl2=xData*p(1)+p(2);
                    mslplot=plot(xData,msl,'k')
                    hold on
                    mslplot2=plot(xData,msl2,'Color','r')
                    mslplot3=plot(xData,raw_msl,'Color', [0.8 0.8 0.8])
                    xlabel('year')
                    ylabel('Mean SSH (m)')
%                     mean_trend=ncread(filename,'mean_trend');
                    title([regionname, ', ', num2str(nyear), 'y low-passed MSL(', ...
                        num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
                    datetick('x','yyyy','keepticks')
                    axis tight;
                    ylim(meanplotlev)
                    set(mslplot,'LineWidth',2);
                    set(gca,'FontSize',15);
                    grid on
                    hold off
                    saveas(gcf,jpgname,'jpg');
                    grid off
                    close all;
                end
            end
        end
% end-------------------- lowpassed msl time series

% start-------------------- detrended lowpassed msl time series
        nyears=[2,3,5];
        nc_varname_prefixes={'cmems_detrended_', 'interped_detrended_'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', testname, '_',regionname, '_',nc_varname,'_msl_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
                    for varind=1:length(inputyear)*12
                        eval(['tempmsl=squeeze(ncread(filename,', '''', nc_varname, '''', ',[1,1,varind], [inf inf 1]));']);
                        msl(varind)=mean(tempmsl(:),'omitnan');
                        if (strcmp(nc_varname_prefix, 'interped_detrended_')==1)
                            temp_raw_msl=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
                        elseif (strcmp(nc_varname_prefix, 'cmems_detrended_')==1)
                            temp_raw_msl=ncread(filename,'cmems_sla',[1 1 varind], [inf inf 1]);
                        end
                        raw_msl(varind)=squeeze(mean(temp_raw_msl(:),'omitnan'));
                    end
                    msl=msl-mean(msl, 'omitnan');    
                    raw_msl=raw_msl-mean(raw_msl(:));

                    p=polyfit(xData,raw_msl,1);
                    msl2=xData*p(1)+p(2);
                    mslplot=plot(xData,msl,'k')
                    hold on
%                     mslplot2=plot(xData,msl2,'Color','r')
                    mslplot2=plot(xData,raw_msl-msl2,'Color', [0.8 0.8 0.8])
                    xlabel('year')
                    ylabel('Mean SSH (m)')
%                     mean_trend=ncread(filename,'mean_trend');
                    title([regionname, ', ', num2str(nyear), 'y low-passed det MSL(', ...
                        num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
                    datetick('x','yyyy','keepticks')
                    axis tight;
                    ylim(meanplotlev)
                    set(mslplot,'LineWidth',2);
                    set(gca,'FontSize',15);
                    grid on
                    hold off
                    saveas(gcf,jpgname,'jpg');
                    grid off
                    close all;
                end
            end
        end
% end-------------------- detrended lowpassed msl time series



% start-------------------- climatological msl time series
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        for monthij=1:12
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
                if (exist('msl')==0)
                    for varind=1:length(inputyear)*12
                        msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
                    end
                    climmsl=reshape(msl,[12,length(inputyear)]);
                else
                    climmsl=reshape(msl,[12,length(inputyear)]);
                end
                tempmsl=squeeze(climmsl(monthij,:));
                tempmsl=tempmsl-mean(tempmsl);
                p=polyfit(inputyear,tempmsl,1);
                msl2=inputyear*p(1)+p(2);
                mslplot=plot(inputyear,tempmsl,'k')
                hold on
                mslplot2=plot(inputyear,msl2,'Color','r')
                xlabel('year')
                ylabel('Mean SSH (m)')
                title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                    ', ', calendarname{monthij}(1:3), '), ',num2str(round(mean_clim_trend(monthij)*1000,2)), ' mm/y'])
%                 datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
        end
% end-------------------- climatological msl time series

% start-------------------- climatological msl trend
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_trend_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
        if abs(mean_clim_trend(1))<0.01
            mean_clim_trend=mean_clim_trend*1000.0  %% m -> mm
        end
        if (exist(jpgname , 'file') ~= 2)
            mslplot=plot(1:12,mean_clim_trend,'k')
            hold on
            xlabel('month')
            ylabel('trend (mm/yr)')
            title([regionname, ', climatological trend(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
            axis tight;
            ylim(trendplotlev)
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
% end-------------------- climatological msl trend

% start-------------------- seasonal msl
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_msl_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
        
%         clim_ssh=ncread(filename,'clim_ssh');
%         for ijij=1:12
%             temp_clim_ssh=clim_ssh(:,:,ijij);
%             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
%         end
%         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
        
        if (exist(jpgname , 'file') ~= 2)
            comb_interped_data=ncread(filename,'interped_ssh');
            lon_cmems=ncread(filename, 'lon_cmems');
            lat_cmems=ncread(filename, 'lat_cmems');
            len_lon=length(lon_cmems);
            len_lat=length(lat_cmems);
            comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
            clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
            cmems_trend=ncread(filename, 'cmems_trend');
            cmems_mask=ones(size(cmems_trend));
            cmems_mask(isnan(cmems_trend))=NaN;
            clim_ssh = clim_ssh .* cmems_mask;
            for ijij=1:12
                temp_clim_ssh=clim_ssh(:,:,ijij);
                mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
            end
            mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
            
            mslplot=plot(1:12,mean_clim_ssh,'k')
            hold on
            xlabel('month')
            ylabel('trend (mm/yr)')
            title([regionname, ', seasonal msl(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',20);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
% end-------------------- seasonal msl


% start-------------------- cmems seasonal msl
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_cmems_msl_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
        
%         clim_ssh=ncread(filename,'clim_ssh');
%         for ijij=1:12
%             temp_clim_ssh=clim_ssh(:,:,ijij);
%             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
%         end
%         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
        
        if (exist(jpgname , 'file') ~= 2)
            comb_interped_data=ncread(filename,'interped_ssh');
            lon_cmems=ncread(filename, 'lon_cmems');
            lat_cmems=ncread(filename, 'lat_cmems');
            len_lon=length(lon_cmems);
            len_lat=length(lat_cmems);
            clim_ssh=ncread(filename, 'clim_cmems_ssh');
            clim_ssh(clim_ssh>1000000)=NaN;
            cmems_trend=ncread(filename, 'cmems_trend');
            cmems_mask=ones(size(cmems_trend));
            cmems_mask(isnan(cmems_trend))=NaN;
            clim_ssh = clim_ssh .* cmems_mask;
            for ijij=1:12
                temp_clim_ssh=clim_ssh(:,:,ijij);
                mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
            end
            mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
            
            mslplot=plot(1:12,mean_clim_ssh,'k')
            hold on
            xlabel('month')
            ylabel('trend (mm/yr)')
            title([regionname, ', seasonal msl(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',20);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
% end-------------------- cmems seasonal msl

% start-------------------- cmems seasonal msl amplitude
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_amp_cmems_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
        
%         clim_ssh=ncread(filename,'clim_ssh');
%         for ijij=1:12
%             temp_clim_ssh=clim_ssh(:,:,ijij);
%             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
%         end
%         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
        
        if (exist(jpgname , 'file') ~= 2)
            comb_interped_data=ncread(filename,'interped_ssh');
            lon_cmems=ncread(filename, 'lon_cmems');
            lat_cmems=ncread(filename, 'lat_cmems');
            len_lon=length(lon_cmems);
            len_lat=length(lat_cmems);
            clim_ssh=ncread(filename, 'clim_cmems_ssh');
            clim_ssh(clim_ssh>1000000)=NaN;
            cmems_trend=ncread(filename, 'cmems_trend');
            cmems_mask=ones(size(cmems_trend));
            cmems_mask(isnan(cmems_trend))=NaN;
            clim_ssh = clim_ssh .* cmems_mask;
            for ijij=1:12
                temp_clim_ssh=clim_ssh(:,:,ijij);
                mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
            end
            mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
            for loni=1:size(clim_ssh,1)
                for lati=1:size(clim_ssh,2)
                    amp_clim_ssh(loni,lati)=(max(clim_ssh(loni,lati,:))-min(clim_ssh(loni,lati,:)))/2.0;
                end
            end
%             pcolor(amp_clim_ssh'*100)
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(lon_cmems',lat_cmems',squeeze(amp_clim_ssh(:,:)'*100.0));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('seasonal amp, ','cmems',',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(amp_clim_ssh(:)*100.0, 'omitnan'),2)));  

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([2 16]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'cmems seasonal msl amplitude', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- cmems seasonal msl amplitude


% start-------------------- seasonal msl amplitude
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_amp_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
        
%         clim_ssh=ncread(filename,'clim_ssh');
%         for ijij=1:12
%             temp_clim_ssh=clim_ssh(:,:,ijij);
%             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
%         end
%         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
        
        if (exist(jpgname , 'file') ~= 2)
            
            lon_cmems=ncread(filename, 'lon_cmems');
            lat_cmems=ncread(filename, 'lat_cmems');
            len_lon=length(lon_cmems);
            len_lat=length(lat_cmems);
            comb_interped_data=ncread(filename,'interped_ssh');
            comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
            clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
            clim_ssh(clim_ssh>1000000)=NaN;
            cmems_trend=ncread(filename, 'cmems_trend');
            cmems_mask=ones(size(cmems_trend));
            cmems_mask(isnan(cmems_trend))=NaN;
            clim_ssh = clim_ssh .* cmems_mask;
            for ijij=1:12
                temp_clim_ssh=clim_ssh(:,:,ijij);
                mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
            end
            mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
            for loni=1:size(clim_ssh,1)
                for lati=1:size(clim_ssh,2)
                    amp_clim_ssh(loni,lati)=(max(clim_ssh(loni,lati,:))-min(clim_ssh(loni,lati,:)))/2.0;
                end
            end
%             pcolor(amp_clim_ssh'*100)
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(lon_cmems',lat_cmems',squeeze(amp_clim_ssh(:,:)'*100.0));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('seasonal amp, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(amp_clim_ssh(:)*100.0, 'omitnan'),2)));  

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis([2 16]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'cmems seasonal msl amplitude', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- seasonal msl amplitude


% start-------------------- seasonal msl amplitude diff
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_amp_diff', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
        
%         clim_ssh=ncread(filename,'clim_ssh');
%         for ijij=1:12
%             temp_clim_ssh=clim_ssh(:,:,ijij);
%             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
%         end
%         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
        
        if (exist(jpgname , 'file') ~= 2)
            
            lon_cmems=ncread(filename, 'lon_cmems');
            lat_cmems=ncread(filename, 'lat_cmems');
            len_lon=length(lon_cmems);
            len_lat=length(lat_cmems);
            clim_cmems_ssh=ncread(filename, 'clim_cmems_ssh');
            clim_cmems_ssh(clim_cmems_ssh>1000000)=NaN;
            comb_interped_data=ncread(filename,'interped_ssh');
            comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
            clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
            clim_ssh(clim_ssh>1000000)=NaN;
            cmems_trend=ncread(filename, 'cmems_trend');
            cmems_mask=ones(size(cmems_trend));
            cmems_mask(isnan(cmems_trend))=NaN;
            clim_ssh = clim_ssh .* cmems_mask;
            clim_cmems_ssh=clim_cmems_ssh .* cmems_mask;

%             for ijij=1:12
%                 temp_clim_ssh=clim_ssh(:,:,ijij);
%                 mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
%             end
%             mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
            for loni=1:size(clim_ssh,1)
                for lati=1:size(clim_ssh,2)
                    amp_clim_ssh(loni,lati)=(max(clim_ssh(loni,lati,:))-min(clim_ssh(loni,lati,:)))/2.0;
                    amp_clim_cmems_ssh(loni,lati)=(max(clim_cmems_ssh(loni,lati,:))-min(clim_cmems_ssh(loni,lati,:)))/2.0;
                end
            end
%             pcolor(amp_clim_ssh'*100)
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(lon_cmems',lat_cmems',amp_clim_cmems_ssh(:,:)'*100-amp_clim_ssh(:,:)'*100);
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('seasonal amp diff, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(amp_clim_cmems_ssh(:)*100.0, 'omitnan')-mean(amp_clim_ssh(:)*100.0, 'omitnan'),2)));  

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis([2 16]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'cmems seasonal msl amplitude', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- seasonal msl amplitude diff


% start-------------------- get distance weighted msl
%         for i= 1:length(recon_lon)
%             for j=1:length(recon_lat)
%                 for l=1:size(comb_data_filtered,3)
%                     comb_data_filtered_weight(i,j,l)=comb_data_filtered(i,j,l).*cos(recon_lat(j)/180.0*pi);
%                 end
%             end
%         end
%         
%         temp_varvar=comb_data_filtered_weight(:,:,1);
%         for i= 1:length(recon_lon)
%             for j=1:length(recon_lat)
%                 temp_varvar(i,j)=cos(recon_lat(j)/180.0*pi).*temp_varvar(i,j)./temp_varvar(i,j);
%             end
%         end
%         for i= 1:length(recon_lon)
%             for j=1:length(recon_lat)
%                 temp_varvar(isnan(temp_varvar))=NaN;
%             end
%         end
%         m_cidfw=squeeze(sum(sum(comb_data_filtered_weight,1,'omitnan'),2,'omitnan'));
%         m_cidfw=m_cidfw./sum(sum(temp_varvar,1,'omitnan'),2,'omitnan');
% end-------------------- get distance weighted msl        
        



% start-------------------- climatological msl time series

        for i =1:length(inputyear) 
            tempyear=inputyear(i);
    %         for month=1:12
                xData2(i) = datenum([num2str(6,'%02i'),'-30-',num2str(tempyear)]);
    %         end
        end
        %     plot(squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan')))
%         isize = size(comb_data_filtered,1)
%         jsize = size(comb_data_filtered,2)
%         lsize = size(comb_data_filtered,3)
%         comb_yearly_data_filtered=reshape(comb_data_filtered,[isize, jsize, 12, lsize/12]);
%         mean_yearly_data_filtered=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
        
%         trendtime=1:length(xData2);
%         p=polyfit(trendtime,mean_yearly_data_filtered(1:length(xData2))',1);
%         yearly_interped_trend=p(1);
%         yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y
% 
%         yearly_msl=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         yearly_msl=yearly_msl-mean(yearly_msl);    
%         p=polyfit(xData2,yearly_msl',1);
%         yearly_msl2=xData2*p(1)+p(2);
%         yearly_mslplot=plot(xData2,yearly_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
%         hold on
%         yearly_mslplot2=plot(xData2,yearly_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('Year')
%         ylabel('Mean SSH (m)')
%         title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(yearly_interped_trend,2)), ' mm/y'])
%         ylim(meanplotlev)
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(yearly_mslplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
%         close all;



% start-------------------- cmems MSL plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
                    end
                end
                cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data= mean_data .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' mean, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
                caxis(colorbar_lev);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
% end-------------------- cmems MSL plot


% start-------------------- RMS plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));s
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_rms_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                            mean_model_data=mean_data;
                        end
                        mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
                        mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
                    end
                end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
                mean_data= mean_data .* mask_model;
                mean_model_data = mean_model_data .* mask_model;
                cmems_msl=mean(mean_data(:), 'omitnan');
                model_msl=mean(mean_model_data(:), 'omitnan');
%                 mean_data=mean_data-mean(mean_data(:),'omitnan');
%                 mean_model_data=mean_model_data-mean(mean_model_data(:),'omitnan');
%                 mean_rms = sqrt((mean_model_data - mean_data).^2) * 100.0;

                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('sq_diff' , 'var') ~= 1)
                            sq_diff=zeros(size(data));
                        end
                            sq_diff=sq_diff + ((model_data-model_msl)*100-(data-cmems_msl)*100).^2;
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_rms=sqrt(sq_diff./(length(inputyear) * length(inputmonth)));
                mean_rms=mean_rms.*mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([5 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data sq_diff
            end
% end-------------------- RMS plot


% start-------------------- RMS_corrected plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));s
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_rms_corrected_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'corrected_interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                            mean_model_data=mean_data;
                        end
                        mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
                        mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
                    end
                end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
                mean_data= mean_data .* mask_model;
                mean_model_data = mean_model_data .* mask_model;
                cmems_msl=mean(mean_data(:), 'omitnan');
                model_msl=mean(mean_model_data(:), 'omitnan');
%                 mean_data=mean_data-mean(mean_data(:),'omitnan');
%                 mean_model_data=mean_model_data-mean(mean_model_data(:),'omitnan');
%                 mean_rms = sqrt((mean_model_data - mean_data).^2) * 100.0;

                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'corrected_interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('sq_diff' , 'var') ~= 1)
                            sq_diff=zeros(size(data));
                        end
                            sq_diff=sq_diff + ((model_data-model_msl)*100-(data-cmems_msl)*100).^2;
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_rms=sqrt(sq_diff./(length(inputyear) * length(inputmonth)));
                mean_rms=mean_rms.*mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([5 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data sq_diff
            end
% end-------------------- RMS_corrected plot

% start-------------------- BIAS plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_bias_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                            mean_model_data=mean_data;
                        end
                        mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
                        mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
                    end
                end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
                mean_data= mean_data .* mask_model;
                mean_model_data = mean_model_data .* mask_model;
                cmems_msl=mean(mean_data(:), 'omitnan');
                model_msl=mean(mean_model_data(:), 'omitnan');
%                 mean_data=mean_data-mean(mean_data(:),'omitnan');
%                 mean_model_data=mean_model_data-mean(mean_model_data(:),'omitnan');
%                 mean_rms = sqrt((mean_model_data - mean_data).^2) * 100.0;

                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_bias' , 'var') ~= 1)
                            mean_bias=zeros(size(data));
                        end
                            mean_bias=mean_bias + ((model_data-model_msl)*100-(data-cmems_msl)*100)/(length(inputyear) * length(inputmonth));
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_bias=mean_bias.*mask_model;
               
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_bias');
%                 m_pcolor(cut_lon_rho(1:30,1:20)',cut_lat_rho(1:30,1:20)',mean_bias(1:30,1:20)');


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' bias',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_bias(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([-15 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data mean_bias
            end
% end-------------------- BIAS plot

% start-------------------- cmems variance plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
%                 for yearij=1:length(inputyear)
%                     tempyear=inputyear(yearij);
%                     yearstr=num2str(tempyear, '%04i');
%                     for monthij=1:length(inputmonth)
%                         tempmonth=inputmonth(monthij);
%                         monthstr=num2str(tempmonth, '%02i');
%                         varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
%                         data(:,:)=ncread(filename,'cmems_sla',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         if (exist('mean_data' , 'var') ~= 1)
%                             mean_data=zeros(size(data));
%                         end
%                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
%                     end
%                 end
                cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                for loni=1:size(cmems_sla,1)
                    for lati=1:size(cmems_sla,2)
                        cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:).*100);
                    end
                end
%                 cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
%                 mean_data= mean_data .* mask_model;
                cmems_sla_var=cmems_sla_var .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',cmems_sla_var');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' VAR, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis([0 2]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
% end-------------------- cmems variance plot


% start-------------------- cmems std plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'STD','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
%                 for yearij=1:length(inputyear)
%                     tempyear=inputyear(yearij);
%                     yearstr=num2str(tempyear, '%04i');
%                     for monthij=1:length(inputmonth)
%                         tempmonth=inputmonth(monthij);
%                         monthstr=num2str(tempmonth, '%02i');
%                         varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
%                         data(:,:)=ncread(filename,'cmems_sla',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         if (exist('mean_data' , 'var') ~= 1)
%                             mean_data=zeros(size(data));
%                         end
%                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
%                     end
%                 end
                cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                for loni=1:size(cmems_sla,1)
                    for lati=1:size(cmems_sla,2)
                        cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:).*100);
                    end
                end
%                 cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
%                 mean_data= mean_data .* mask_model;
                cmems_sla_var=cmems_sla_var .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',sqrt(cmems_sla_var)');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' STD, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis([0 2]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
% end-------------------- cmems std plot



% start-------------------- cmems variance_seasonal plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '_', num2str(tempmonth,'%02i'),'.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2)        
                    run(param_script);
    %                 for yearij=1:length(inputyear)
    %                     tempyear=inputyear(yearij);
    %                     yearstr=num2str(tempyear, '%04i');
    %                     for monthij=1:length(inputmonth)
    %                         tempmonth=inputmonth(monthij);
    %                         monthstr=num2str(tempmonth, '%02i');
    %                         varind=((yearij-1)*12)+monthij
                            if (exist('lon_min' , 'var') ~= 1)
                                lon_cmems=ncread(filename, 'lon_cmems');
                                lat_cmems=ncread(filename, 'lat_cmems');
                                [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
    %                         data_info = ncinfo(filename, varname); 
    %                         data(:,:)=ncread(filename,'cmems_sla',[lon_min(1) lat_min(1) varind], ...
    %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
    %                         if (exist('mean_data' , 'var') ~= 1)
    %                             mean_data=zeros(size(data));
    %                         end
    %                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
    %                     end
    %                 end
                    for yearij=1:length(inputyear)
                        cmems_sla(:,:,yearij)=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) (yearij-1)*12+tempmonth], ...
                                    [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                    end

                    for loni=1:size(cmems_sla,1)
                        for lati=1:size(cmems_sla,2)
                            cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:).*100);
                        end
                    end
    %                 cmems_msl=mean(mean_data(:), 'omitnan');
                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
    %                 mean_data= mean_data .* mask_model;
                    cmems_sla_var=cmems_sla_var .* mask_model;
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                    m_pcolor(cut_lon_rho',cut_lat_rho',cmems_sla_var');

                    shading(gca,m_pcolor_shading_method);   

                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat(variable, ' VAR, ','cmems,', calendarname{tempmonth}(1:3), ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
                    caxis([0 1]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    close all;
                    clear lon_rho mean_data
                end
            end
% end-------------------- cmems variance_seasonal plot


% start-------------------- model variance plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(filename, 'lon_cmems');
                    lat_cmems=ncread(filename, 'lat_cmems');
                    [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                    cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end

                model_sla=ncread(filename,'interped_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_var(loni,lati)=var(model_sla(loni,lati,:).*100);
                    end
                end
%                 cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
%                 mean_data= mean_data .* mask_model;
                model_sla_var=model_sla_var .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_var');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' VAR, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 2]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
% end-------------------- model variance plot

% start-------------------- model variance_seasonal plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                pngname=strcat(outfile, '_', testname,'_',regionname, '_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'),'_', num2str(tempmonth,'%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2)        
                    run(param_script);
                    if (exist('lon_min' , 'var') ~= 1)
                        lon_cmems=ncread(filename, 'lon_cmems');
                        lat_cmems=ncread(filename, 'lat_cmems');
                        [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end

    %                 model_sla=ncread(filename,'interped_sla', [lon_min(1) lat_min(1) 1], ...
    %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
    %                         
                    for yearij=1:length(inputyear)
                        model_sla(:,:,yearij)=ncread(filename,'interped_sla', [lon_min(1) lat_min(1) (yearij-1)*12+tempmonth], ...
                                    [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                    end        

                    for loni=1:size(model_sla,1)
                        for lati=1:size(model_sla,2)
                            model_sla_var(loni,lati)=var(model_sla(loni,lati,:).*100);
                        end
                    end
    %                 cmems_msl=mean(mean_data(:), 'omitnan');
                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
    %                 mean_data= mean_data .* mask_model;
                    model_sla_var=model_sla_var .* mask_model;
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                    m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_var');

                    shading(gca,m_pcolor_shading_method);   

                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat(variable, ' VAR, ',testname,',', calendarname{tempmonth}(1:3), ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
                    caxis([0 1]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    close all;
                    clear lon_rho mean_data
                end
            end
% end-------------------- model variance_seasonal plot


% start-------------------- cmems variance_lowpass plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            nyears=[1:5];
            nc_varname_prefixes={'cmems_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                
                    pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', nc_varname,'_','VAR_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2)        
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['cmems_sla', '=ncread(filename,', '''', nc_varname,'''',');']);
%                         cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
%                                     [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                        cmems_sla=cmems_sla(:,:,nyear*6+1:end-(nyear*6)-1);
                        for loni=1:size(cmems_sla,1)
                            for lati=1:size(cmems_sla,2)
                                cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:));
                            end
                        end
        %                 cmems_msl=mean(mean_data(:), 'omitnan');
                        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                        hold on;

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
        %                 mean_data= mean_data .* mask_model;
                        cmems_sla_var=cmems_sla_var .* mask_model;
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                        m_pcolor(cut_lon_rho',cut_lat_rho',cmems_sla_var'.*100);

                        shading(gca,m_pcolor_shading_method);   

                        m_gshhs_i('color',m_gshhs_line_color)  
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                        titlename = strcat(nc_varname, ' VAR, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(jet);
                        set(h,'fontsize',colorbar_fontsize);
                        title(h,'(cm)','fontsize',colorbar_title_fontsize);
                        caxis([0 1]);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                        saveas(gcf,pngname,'tif');
                        close all;
                        clear lon_rho mean_data
                    end
                end
            end
% end-------------------- cmems variance_lowpass plot

% start-------------------- model variance_lowpass plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            nyears=[1:5];
            nc_varname_prefixes={'interped_sla_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                
                    pngname=strcat(outfile, '_', testname,'_',regionname, '_', nc_varname,'_','VAR_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2)        
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['model_sla', '=ncread(filename,', '''', nc_varname,'''',');']);
%                         cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
%                                     [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                        model_sla=model_sla(:,:,nyear*6+1:end-(nyear*6)-1);
                        for loni=1:size(model_sla,1)
                            for lati=1:size(model_sla,2)
                                model_sla_var(loni,lati)=var(model_sla(loni,lati,:));
                            end
                        end
        %                 cmems_msl=mean(mean_data(:), 'omitnan');
                        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                        hold on;

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
        %                 mean_data= mean_data .* mask_model;
                        model_sla_var=model_sla_var .* mask_model;
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                        m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_var'.*100);

                        shading(gca,m_pcolor_shading_method);   

                        m_gshhs_i('color',m_gshhs_line_color)  
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                        titlename = strcat(nc_varname, ' VAR, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(jet);
                        set(h,'fontsize',colorbar_fontsize);
                        title(h,'(cm)','fontsize',colorbar_title_fontsize);
                        caxis([0 1]);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                        saveas(gcf,pngname,'tif');
                        close all;
                        clear lon_rho mean_data
                    end
                end
            end
% end-------------------- model variance_lowpass plot
    end
    
% start-------------------- model pattern correlation

cmems_adt=ncread(filename,'cmems_adt');
interped_ssh=ncread(filename,'interped_ssh');
cmems_sla=ncread(filename,'cmems_sla');
interped_sla=ncread(filename,'interped_sla');
m_cmems_adt=mean(cmems_adt,3,'omitnan');
m_interped_ssh=mean(interped_ssh,3,'omitnan');
m_cmems_sla=mean(cmems_sla,3,'omitnan');
m_interped_sla=mean(interped_sla,3,'omitnan');
m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));

corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))))
corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))))
corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
    cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)

% end-------------------- model pattern correlation


% start-------------------- model pattern correlation (yearly)
    pngname=strcat(outfile, '_', testname,'_',regionname, '_', nc_varname,'_','yearly_pattern_corr_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    if (exist(pngname , 'file') ~= 2)        
        for yearij=1:length(inputyear)
            tempyear=inputyear(yearij);
            
            yearly_cmems_adt=ncread(filename,'cmems_adt',[1 1 (yearij-1)*12+1], [inf inf 12]);
            yearly_interped_ssh=ncread(filename,'interped_ssh',[1 1 (yearij-1)*12+1], [inf inf 12]);
            yearly_cmems_sla=ncread(filename,'cmems_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);
            yearly_interped_sla=ncread(filename,'interped_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);
            
            m_cmems_adt=mean(yearly_cmems_adt,3,'omitnan');
            m_interped_ssh=mean(yearly_interped_ssh,3,'omitnan');
            m_cmems_sla=mean(yearly_cmems_sla,3,'omitnan');
            m_interped_sla=mean(yearly_interped_sla,3,'omitnan');
%             m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
%             m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));

            ssh_corr(yearij) = corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))));
            sla_corr(yearij) = corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))));
%             corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
%                 cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)
        end
        plot(ssh_corr);
        hold on
        plot(sla_corr);
        hold off
        close all
    end
% end-------------------- model pattern correlation (yearly)


end