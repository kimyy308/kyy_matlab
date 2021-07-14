close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
all_testname2 = {'test11', 'test12'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

all_region2 ={'AKP4'};
% all_region2 ={'YSECS', 'ECS2', 'ES', 'YS', 'NES', 'SES'}

% all_region2 ={'AKP4'}
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
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

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
            case('YSECS') %% Yellow sea and East China Sea
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
        
        cmems_filename = ['E:\Data\Observation\CMEMS\',regionname, ...
        'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc']
        comb_cmems_err = ncread(cmems_filename, 'cmems_err');
        comb_cmems_data = ncread(cmems_filename, 'cmems_sla');
        m_err=mean(comb_cmems_err,3);
        mean_1y=comb_cmems_data(:,:,1:12);
        mean_1y=mean(mean_1y(:), 'omitnan');
        mean_ly=comb_cmems_data(:,:,end-12:end);
        mean_ly=mean(mean_ly(:), 'omitnan');

        for y_ind=1:length(inputyear)
            yearly_sealevel(:,:,y_ind)=mean(comb_cmems_data(:,:,(y_ind-1)*12+1:y_ind*12),3,'omitnan');
        end
        clear var
        for xind=1:size(yearly_sealevel,1)
            for yind=1:size(yearly_sealevel,2)                
                var_yearly_sealevel(xind,yind)=var(yearly_sealevel(xind,yind,:));
            end
        end
        diff_all_y=mean(sqrt(var_yearly_sealevel(:)).*100, 'omitnan');
        m_err(m_err.*100<diff_all_y.*2)=NaN;
        err_mask = cmems_mask;
        err_mask(isfinite(m_err))=NaN;
        
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

        cmems_sla = ncread(filename,'cmems_sla');
        for sla_i=1:size(cmems_sla,1)
            for sla_j=1:size(cmems_sla,2)
                cmems_sla_mean(sla_i,sla_j)=mean(cmems_sla(sla_i,sla_j,:),'omitnan');
            end
        end
        interped_ssh=ncread(filename,'interped_ssh');
        for sla_i=1:size(cmems_sla,1)
            for sla_j=1:size(cmems_sla,2)
                interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
            end
        end
        
        
        interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
        clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
        for t=1:length(inputyear)
            interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
        end
        cmems_sla_divided=reshape(cmems_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
        clim_cmems_sla=mean(cmems_sla_divided,4,'omitnan');
        for t=1:length(inputyear)
            cmems_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(cmems_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_cmems_sla);
        end
        
        
%         plot(squeeze(interped_sla_filtered(20,20,:)))
%         hold on 
%         plot(squeeze(cmems_sla_filtered(20,20,:)))
%         hold off



% start-------------------- lowpass SLA RMS plot
        nyears=[1:5];
        nc_varname_prefixes={'interped_sla_',};
        nc_titlename_prefixes={'lowpass'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                cmems_varname=['cmems_', num2str(nyear),'y_lowpass'];
                pngname=strcat(outfile, '_', testname,'_',regionname, '_err_sla_RMS_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2)
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
                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        eval(['model_data(:,:)', '=ncread(filename,', '''', nc_varname, '''', ',[1,1,varind],[inf,inf,1]', ');']);
                        eval(['data(:,:)', '=ncread(filename,','''', cmems_varname, '''',',[1,1,varind],[inf,inf,1]',');']);

%                         data(:,:) = cmems_sla(:,:,varind);
%                         model_data(:,:) = interped_sla(:,:,varind);
                        if(varind > nyear*6 && varind <= length(inputyear)*12-nyear*6)
                            if (exist('sq_diff' , 'var') ~= 1)
                                sq_diff=zeros(size(data));
                            end
                            sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                        else
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                        end
                    end
                end

%                     eval([nc_varname, '=ncread(filename,','''', nc_varname,'''',');']);
%                     eval([cmems_varname, '=',cmems_varname, '.*cmems_mask;']);

                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
%                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                    n_lp=(length(inputyear) * length(inputmonth) - nyear*12);
                    mean_rms=sqrt(sq_diff./n_lp);
                    mean_rms=mean_rms.*err_mask;
                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');
        
                    shading(gca,m_pcolor_shading_method);   
                
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat('LP-',num2str(nyear),'y-SLA', ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
    %                 caxis(colorbar_lev);
                    caxis([0 15]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    close all;
                    
                    matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_',nc_varname, '_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    save(matname, 'mean_rms', 'sq_diff_all');
                    
                    clear lon_rho mean_data sq_diff sq_diff_all
                    
                    
                end
            end
        end
% end-------------------- lowpass SLA RMS plot


% start-------------------- lowpass SLA RMS_corrected plot
        nyears=[1:5];
        nc_varname_prefixes={'corrected_interped_sla_'};
        nc_titlename_prefixes={'lowpass'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                cmems_varname=['cmems_', num2str(nyear),'y_lowpass'];
                pngname=strcat(outfile, '_', testname,'_',regionname, '_err_sla_RMS_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2)
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
                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        eval(['model_data(:,:)', '=ncread(filename,', '''', nc_varname, '''', ',[1,1,varind],[inf,inf,1]', ');']);
                        eval(['data(:,:)', '=ncread(filename,','''', cmems_varname, '''',',[1,1,varind],[inf,inf,1]',');']);

%                         data(:,:) = cmems_sla(:,:,varind);
%                         model_data(:,:) = interped_sla(:,:,varind);
                        if(varind > nyear*6 && varind <= length(inputyear)*12-nyear*6)
                            if (exist('sq_diff' , 'var') ~= 1)
                                sq_diff=zeros(size(data));
                            end
                            sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                        else
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                        end
                    end
                end

%                     eval([nc_varname, '=ncread(filename,','''', nc_varname,'''',');']);
%                     eval([cmems_varname, '=',cmems_varname, '.*cmems_mask;']);

                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
%                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                    n_lp=(length(inputyear) * length(inputmonth) - nyear*12);
                    mean_rms=sqrt(sq_diff./n_lp);
                    mean_rms=mean_rms.*err_mask;
                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');
        
                    shading(gca,m_pcolor_shading_method);   
                
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat('LP-',num2str(nyear),'y-SLA', ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
    %                 caxis(colorbar_lev);
                    caxis([0 15]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    close all;
                    
                    matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_',nc_varname, '_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    save(matname, 'mean_rms', 'sq_diff_all');
                    
                    clear lon_rho mean_data sq_diff sq_diff_all
                    
                    
                end
            end
        end
% end-------------------- lowpass SLA RMS_corrected plot


% start-------------------- lowpass SLA RMS_corrected plot (diff)
        nyears=[2,5];
        nc_varname_prefixes={'corrected_interped_sla_'};
        nc_varname_prefixes2={'interped_sla_'};
        nc_titlename_prefixes={'lowpass-c'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                cmems_varname=['cmems_', num2str(nyear),'y_lowpass'];
                pngname=strcat(outfile, '_', 'diff','_',regionname, '_err_sla_RMS_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2)
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
                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                    '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                    '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
   
                        eval(['model_data(:,:)', '=ncread(diffname,', '''', nc_varname2, '''', ',[1,1,varind],[inf,inf,1]', ');']);
                        eval(['model_data2(:,:)', '=ncread(diffname2,', '''', nc_varname, '''', ',[1,1,varind],[inf,inf,1]', ');']);                        
                        eval(['data(:,:)', '=ncread(filename,','''', cmems_varname, '''',',[1,1,varind],[inf,inf,1]',');']);

%                         data(:,:) = cmems_sla(:,:,varind);
%                         model_data(:,:) = interped_sla(:,:,varind);
                        if(varind > nyear*6 && varind <= length(inputyear)*12-nyear*6)
                            if (exist('sq_diff' , 'var') ~= 1)
                                sq_diff=zeros(size(data));
                                sq_diff2=zeros(size(data));
                            end
                            sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                            
                            sq_diff2=sq_diff2 + ((model_data2)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff2_all(:,:,varind)=((model_data2)*100-(data)*100).^2;
                        else
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                            sq_diff2_all(:,:,varind)=((model_data2)*100-(data)*100).^2;
                        end
                    end
                end

%                     eval([nc_varname, '=ncread(filename,','''', nc_varname,'''',');']);
%                     eval([cmems_varname, '=',cmems_varname, '.*cmems_mask;']);

                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
%                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                    n_lp=(length(inputyear) * length(inputmonth) - nyear*12);
                    mean_rms=sqrt(sq_diff./n_lp);
                    mean_rms=mean_rms.*err_mask;
                    mean_rms2=sqrt(sq_diff2./n_lp);
                    mean_rms2=mean_rms2.*err_mask;
                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms'-mean_rms2');
                    
                    mean_rms_diff=mean_rms-mean_rms2;
                    shading(gca,m_pcolor_shading_method);   
                
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat('LP-',num2str(nyear),'y-SLA', ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:)-mean_rms2(:),'omitnan'),1)));  
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
    %                 caxis(colorbar_lev);
                    caxis([-2 2]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    close all;
                    
                    matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_diff_',nc_varname, '_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    save(matname, 'mean_rms_diff');
                    
                    clear lon_rho mean_data sq_diff sq_diff_all
                    
                    
                end
            end
        end
% end-------------------- lowpass SLA RMS_corrected plot

% start-------------------- lowpass SLA RMS_corrected plot (diff, normalized)
        nyears=[2,5];
        nc_varname_prefixes={'corrected_interped_sla_'};
        nc_varname_prefixes2={'interped_sla_'};
        nc_titlename_prefixes={'lowpass-c'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                cmems_varname=['cmems_', num2str(nyear),'y_lowpass'];
                pngname=strcat(outfile, '_', 'diff','_',regionname, '_err_sla_RMS_diff_norm', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%                 if (exist(pngname , 'file') ~= 2)
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
                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                    '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                    '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
   
                        eval(['model_data(:,:)', '=ncread(diffname,', '''', nc_varname2, '''', ',[1,1,varind],[inf,inf,1]', ');']);
                        eval(['model_data2(:,:)', '=ncread(diffname2,', '''', nc_varname, '''', ',[1,1,varind],[inf,inf,1]', ');']);                        
                        eval(['data(:,:)', '=ncread(filename,','''', cmems_varname, '''',',[1,1,varind],[inf,inf,1]',');']);

%                         data(:,:) = cmems_sla(:,:,varind);
%                         model_data(:,:) = interped_sla(:,:,varind);
                        if(varind > nyear*6 && varind <= length(inputyear)*12-nyear*6)
                            if (exist('sq_diff' , 'var') ~= 1)
                                sq_diff=zeros(size(data));
                                sq_diff2=zeros(size(data));
                            end
                            sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                            
                            sq_diff2=sq_diff2 + ((model_data2)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff2_all(:,:,varind)=((model_data2)*100-(data)*100).^2;
                        else
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                            sq_diff2_all(:,:,varind)=((model_data2)*100-(data)*100).^2;
                        end
                    end
                end

%                     eval([nc_varname, '=ncread(filename,','''', nc_varname,'''',');']);
%                     eval([cmems_varname, '=',cmems_varname, '.*cmems_mask;']);

                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
%                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                    n_lp=(length(inputyear) * length(inputmonth) - nyear*12);
                    mean_rms=sqrt(sq_diff./n_lp);
                    mean_rms=mean_rms.*err_mask;
                    mean_rms2=sqrt(sq_diff2./n_lp);
                    mean_rms2=mean_rms2.*err_mask;
                    
                    mean_rms_diff=mean_rms-mean_rms2;
                    std_mean_rms_diff=std(mean_rms_diff(:), 'omitnan');
                    n_mean_rms_diff=mean_rms_diff/std_mean_rms_diff;
                                        
                    m_pcolor(cut_lon_rho',cut_lat_rho',n_mean_rms_diff');
                    
                    shading(gca,m_pcolor_shading_method);   
                
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat('LP-',num2str(nyear),'y-SLA', ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', ...
                        'M=', num2str(round(mean(n_mean_rms_diff(:),'omitnan'),2)));  
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
    %                 caxis(colorbar_lev);
                    caxis([-6 6]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    close all;
                    
                    matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_diff_',nc_varname, '_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    save(matname, 'mean_rms_diff');
                    
                    clear lon_rho mean_data sq_diff sq_diff_all
                    
                    
%                 end
            end
        end
% end-------------------- lowpass SLA RMS_corrected plot (diff, normalized)



% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
% end-------------------- make timedata for time series  


% start-------------------- SLA_RMS plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_err_rms_', 'sla','_',num2str(min(inputyear),'%04i'), ...
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
% %                         data_info = ncinfo(filename, varname); 
% %                         data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
% %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% %                         model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
% %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         data(:,:) = cmems_sla(:,:,varind);
%                         model_data(:,:) = interped_sla(:,:,varind);
% %                         if (exist('mean_data' , 'var') ~= 1)
%                             mean_data=zeros(size(data));
%                             mean_model_data=mean_data;
%                           end
%                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
%                         mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
%                     end
%                 end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
%                 mean_data= mean_data .* mask_model;
%                 mean_model_data = mean_model_data .* mask_model;
%                 cmems_msl=mean(mean_data(:), 'omitnan');
%                 model_msl=mean(mean_model_data(:), 'omitnan');
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
%                         data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        data(:,:) = cmems_sla(:,:,varind);
                        model_data(:,:) = interped_sla(:,:,varind);
                        if (exist('sq_diff' , 'var') ~= 1)
                            sq_diff=zeros(size(data));
                        end
%                             sq_diff=sq_diff + ((model_data-model_msl)*100-(data-cmems_msl)*100).^2;  %% 100 -> m to cm
                            sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                            
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_rms=sqrt(sq_diff./(length(inputyear) * length(inputmonth)));
                mean_rms=mean_rms.*err_mask;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('SLA', ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([0 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    save(matname, 'mean_rms', 'sq_diff_all');
                
                clear lon_rho mean_data sq_diff sq_diff_all
            end
% end-------------------- SLA_RMS plot

% start-------------------- SLA_RMS_corrected plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_err_rms_corrected_', 'sla','_',num2str(min(inputyear),'%04i'), ...
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

                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

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
%                         data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        data(:,:) = cmems_sla(:,:,varind);
%                         model_data(:,:) = interped_sla(:,:,varind);
                        model_data(:,:)=ncread(filename,'corrected_interped_sla',[lon_min(1) lat_min(1) varind], ...
                                                    [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('sq_diff' , 'var') ~= 1)
                            sq_diff=zeros(size(data));
                        end
%                             sq_diff=sq_diff + ((model_data-model_msl)*100-(data-cmems_msl)*100).^2;  %% 100 -> m to cm
                            sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                            
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_rms=sqrt(sq_diff./(length(inputyear) * length(inputmonth)));
                mean_rms=mean_rms.*err_mask;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('SLA', ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([0 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_corrected_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    save(matname, 'mean_rms', 'sq_diff_all');
                
                clear lon_rho mean_data sq_diff sq_diff_all
            end
% end-------------------- SLA_RMS_corrected plot


% start-------------------- SLA_filtered_RMS plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_err_rms_', 'sla','_filtered_',num2str(min(inputyear),'%04i'), ...
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
% %                         data_info = ncinfo(filename, varname); 
% %                         data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
% %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% %                         model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
% %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         data(:,:) = cmems_sla(:,:,varind);
%                         model_data(:,:) = interped_sla(:,:,varind);
% %                         if (exist('mean_data' , 'var') ~= 1)
%                             mean_data=zeros(size(data));
%                             mean_model_data=mean_data;
%                           end
%                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
%                         mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
%                     end
%                 end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
%                 mean_data= mean_data .* mask_model;
%                 mean_model_data = mean_model_data .* mask_model;
%                 cmems_msl=mean(mean_data(:), 'omitnan');
%                 model_msl=mean(mean_model_data(:), 'omitnan');
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
%                         data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        data(:,:) = cmems_sla_filtered(:,:,varind);
                        model_data(:,:) = interped_sla_filtered(:,:,varind);
                        if (exist('sq_diff' , 'var') ~= 1)
                            sq_diff=zeros(size(data));
                        end
%                             sq_diff=sq_diff + ((model_data-model_msl)*100-(data-cmems_msl)*100).^2;  %% 100 -> m to cm
                            sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
                            sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_rms=sqrt(sq_diff./(length(inputyear) * length(inputmonth)));
                mean_rms=mean_rms.*err_mask;
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('SLA', ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([0 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                    '_rms_filtered_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                    save(matname, 'mean_rms', 'sq_diff_all');
                clear lon_rho mean_data sq_diff sq_diff_all
            end
% end-------------------- SLA_filtered_RMS plot


% 
% % start-------------------- SLA_RMS_seasonal plot
%             figdir2=[figrawdir,'CLIM\'];
%             if (exist(strcat(figdir2) , 'dir') ~= 7)
%                 mkdir(strcat(figdir2));
%             end 
%             outfile = strcat(figdir2,regionname);
%             for monthij=1:length(inputmonth)
%                 tempmonth=inputmonth(monthij);
%                 pngname=strcat(outfile, '_', testname,'_',regionname, '_',num2str(tempmonth, '%02i'),'_err_rms_', 'sla','_',num2str(min(inputyear),'%04i'), ...
%                     '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%                 if (exist(pngname , 'file') ~= 2)        
%                     run(param_script);
% 
%                             if (exist('lon_min' , 'var') ~= 1)
%                                 lon_cmems=ncread(filename, 'lon_cmems');
%                                 lat_cmems=ncread(filename, 'lat_cmems');
%                                 [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
%                                 [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
%                                 cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                                 cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                             end
% 
% 
%                     mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%                     mask_model(mask_model==0)=NaN;
% 
% 
%                     for yearij=1:length(inputyear)
%                         tempyear=inputyear(yearij);
%                         yearstr=num2str(tempyear, '%04i');
%                         for monthij2=1:length(inputmonth)
%                             tempmonth2=inputmonth(monthij2);
%                             monthstr=num2str(tempmonth2, '%02i');
%                             varind=((yearij-1)*12)+monthij2
%                             if (exist('lon_min' , 'var') ~= 1)
%                                 lon_cmems=ncread(filename, 'lon_cmems');
%                                 lat_cmems=ncread(filename, 'lat_cmems');
%                                 [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
%                                 [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
%                                 cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                                 cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                             end
% 
%                             data(:,:) = cmems_sla(:,:,varind);
%                             model_data(:,:) = interped_sla(:,:,varind);
%                             if (exist('sq_diff' , 'var') ~= 1)
%                                 sq_diff=zeros(size(data));
%                             end
%     %                             sq_diff=sq_diff + ((model_data-model_msl)*100-(data-cmems_msl)*100).^2;  %% 100 -> m to cm
%                                 sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
%                                 sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
%                         end
%                     end
% 
% 
% %                         if (exist('mrms_seasonal' , 'var') ~= 1)
%         %                     size_var=size(sq_diff_all_rcm_sat,1);
%                             size_x=size(sq_diff_all,1);
%                             size_y=size(sq_diff_all,2);
%                             size_t=size(sq_diff_all,3);
%                             sq_diff_all_interped_divided = ...
%                                 reshape(sq_diff_all, [size_x, size_y, size_t/12. 12]);
%                             for xij=1:size_x
%                                 for yij=1:size_y
%                                     mrms_interped_seasonal(xij,yij,tempmonth)=sqrt(mean(sq_diff_all_interped_divided(xij,yij,:,tempmonth),3,'omitnan'));
%                                 end
%                             end
%                             mrms_interped_seasonal(:,:,tempmonth)=mrms_interped_seasonal(:,:,tempmonth).*cmems_mask;
%     %                         mrms_interped_seasonal_spatial=mean(sq_diff_all_interped_divided,1,'omitnan');
%     %                         mrms_interped_seasonal_divided=reshape(mrms_interped_seasonal_spatial, [size_t/12, 12]);
% 
% %                             mrms_interped_seasonal(1,tempmonth)=sqrt(mean(mrms_interped_seasonal_divided(:,tempmonth),1,'omitnan'));
% 
% %                         end
% 
% 
%                     m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%                     hold on;    
% % 
% %                     mean_rms=sqrt(sq_diff./(length(inputyear) * length(inputmonth)));
% %                     mean_rms=mean_rms.*err_mask;
%                     
%                     m_pcolor(cut_lon_rho',cut_lat_rho',mrms_interped_seasonal(:,:,tempmonth)');
%                     
%                     temp_mrms=squeeze(mrms_interped_seasonal(:,:,tempmonth));
% 
%                     shading(gca,m_pcolor_shading_method);   
% 
%                     m_gshhs_i('color',m_gshhs_line_color)  
%                     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%                     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
%                     titlename = strcat('SLA', ' rms,', calendarname{tempmonth}(1:3),',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(temp_mrms(:),'omitnan'),1)));  %% + glacier contribution
%                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%                     % set colorbar 
%                     h = colorbar;
%                     colormap(jet);
%                     set(h,'fontsize',colorbar_fontsize);
%                     title(h,'(cm)','fontsize',colorbar_title_fontsize);
%     %                 caxis(colorbar_lev);
%                     caxis([0 15]);
% 
% 
%                     set(gcf, 'PaperUnits', 'points');
%                     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%                     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%                     saveas(gcf,pngname,'tif');
%                     close all;
% %                     matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
% %                         '_rms_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
% %                         save(matname, 'mean_rms', 'sq_diff_all');
% 
%                     clear lon_rho mean_data sq_diff sq_diff_all
%             end
%             end
% % end-------------------- SLA_RMS_seasonal plot
% 
% % start-------------------- lowpass SLA RMS_seasonal plot
%         nyears=[1,2,3,4,5];
%         nc_varname_prefixes={'interped_sla_',};
%         nc_titlename_prefixes={'lowpass'};
%         
%         for nyeari=1:length(nyears)
%             nyear=nyears(nyeari);
%             for nc_varnameij=1:length(nc_varname_prefixes)
%                 nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
%                 nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
%                 cmems_varname=['cmems_', num2str(nyear),'y_lowpass'];
%                 for monthij=1:length(inputmonth)
%                     tempmonth=inputmonth(monthij);
%                     pngname=strcat(outfile, '_', testname,'_',regionname, '_', num2str(tempmonth, '%02i'), '_err_sla_RMS_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
%                         '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%                     if (exist(pngname , 'file') ~= 2)
%                     for yearij=1:length(inputyear)
%                         tempyear=inputyear(yearij);
%                         yearstr=num2str(tempyear, '%04i');
%                         for monthij2=1:length(inputmonth)
%                             tempmonth2=inputmonth(monthij2);
%                             monthstr=num2str(tempmonth2, '%02i');
%                             varind=((yearij-1)*12)+monthij2
%                             if (exist('lon_min' , 'var') ~= 1)
%                                 lon_cmems=ncread(filename, 'lon_cmems');
%                                 lat_cmems=ncread(filename, 'lat_cmems');
%                                 [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
%                                 [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
%                                 cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                                 cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                             end
%                             mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%                             mask_model(mask_model==0)=NaN;
%                             eval(['model_data(:,:)', '=ncread(filename,', '''', nc_varname, '''', ',[1,1,varind],[inf,inf,1]', ');']);
%                             eval(['data(:,:)', '=ncread(filename,','''', cmems_varname, '''',',[1,1,varind],[inf,inf,1]',');']);
% 
%     %                         data(:,:) = cmems_sla(:,:,varind);
%     %                         model_data(:,:) = interped_sla(:,:,varind);
%                             if(varind > nyear*6 && varind < length(inputyear)*12-nyear*6)
%                                 if (exist('sq_diff' , 'var') ~= 1)
%                                     sq_diff=zeros(size(data));
%                                 end
%                                 sq_diff=sq_diff + ((model_data)*100-(data)*100).^2;  %% 100 -> m to cm
%                                 sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
%                             else
%                                 sq_diff_all(:,:,varind)=((model_data)*100-(data)*100).^2;
%                             end
%                         end
%                     end
% 
%     %                     eval([nc_varname, '=ncread(filename,','''', nc_varname,'''',');']);
%     %                     eval([cmems_varname, '=',cmems_varname, '.*cmems_mask;']);
% 
%                         lon_cmems=ncread(filename,'lon_cmems');
%                         lat_cmems=ncread(filename,'lat_cmems');
%                         size_x=size(sq_diff_all,1);
%                         size_y=size(sq_diff_all,2);
%                         size_t=size(sq_diff_all,3);
%                         sq_diff_all_interped_divided = ...
%                             reshape(sq_diff_all, [size_x, size_y, size_t/12, 12]);
%                         for xij=1:size_x
%                             for yij=1:size_y
%                                 mrms_interped_seasonal(xij,yij,tempmonth)=sqrt(mean(sq_diff_all_interped_divided(xij,yij,:,tempmonth),3,'omitnan'));
%                             end
%                         end
%                         mrms_interped_seasonal(:,:,tempmonth)=mrms_interped_seasonal(:,:,tempmonth).*cmems_mask;
%                         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%                         hold on;
%     %                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
%     %                     n_lp=(length(inputyear) * length(inputmonth) - nyear*12);
%     %                     mean_rms=sqrt(sq_diff./n_lp);
%     %                     mean_rms=mean_rms.*err_mask;
% 
% 
%                         m_pcolor(cut_lon_rho',cut_lat_rho',mrms_interped_seasonal(:,:,tempmonth)');
% 
%                         temp_mrms=squeeze(mrms_interped_seasonal(:,:,tempmonth));
% 
%                         shading(gca,m_pcolor_shading_method);   
% 
%                         m_gshhs_i('color',m_gshhs_line_color)  
%                         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%                         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
%                         titlename = strcat('LP-',num2str(nyear),'y-SLA,',calendarname{tempmonth}(1:3), ', rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(temp_mrms(:),'omitnan'),1)));  
%                         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%                         % set colorbar 
%                         h = colorbar;
%                         colormap(jet);
%                         set(h,'fontsize',colorbar_fontsize);
%                         title(h,'(cm)','fontsize',colorbar_title_fontsize);
%         %                 caxis(colorbar_lev);
%                         caxis([0 15]);
% 
%                         set(gcf, 'PaperUnits', 'points');
%                         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%                         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%                         saveas(gcf,pngname,'tif');
%                         close all;
% 
%                         clear lon_rho mean_data sq_diff sq_diff_all
% 
%                     end
%                 end
%             end
%         end
% % end-------------------- lowpass SLA RMS_seasonal plot  

    end
end
