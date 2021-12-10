% %  Updated 07-Dec-2021 by Yong-Yub Kim, 

close all; clear all;  clc;
warning off;

% years_group=[1990, 2000, 2010];
years_group=[2015, 2050];
% years_group=[2020, 2030, 2040, 2050];
seasons_group={'all', 'summer', 'winter'};

for years_groupi=1:length(years_group)
    for seasons_groupi=1:length(seasons_group)

        GCM_info.name = {'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
        GCM_info.abbs = {'GCM-CNE', 'GCM-ECV', 'GCM-ACC', 'GCM-CNH', 'GCM-CMC'};

        % GCM_info.name = {'ACCESS-CM2'};

        % GCM_info.name = {'ens03'};

        % GCM_info.region ={'NWP', 'YS', 'AKP2'}
        % GCM_info.region ={'NWP', 'AKP4'};
        % GCM_info.region ={'NWP', 'AKP4'};
        GCM_info.region ={'AKP4'};

        % GCM_info.region ={'YS'};
        % GCM_info.vars = {'SSH'};
        % GCM_info.vars = {'SST', 'SSS', 'SSH'};
        GCM_info.vars = {'SSH'};
        % GCM_info.vars = {'BT'};

        % GCM_info.region ={'NWP'}
        for testnameind2=1:length(GCM_info.name)
            for regionind2=1:length(GCM_info.region)
                close all;
                clearvars '*' -except regionind2 testnameind2 ...
                    GCM_info.region GCM_info.name GCM_info.vars GCM_info.abbs years_group seasons_group
                % % % 
                tmp.fs=filesep;

               % %     set dropbox path
                if (strcmp(computer,'PCWIN64'))
                    tmp.dropboxpath = 'C:\Users\User\Dropbox';
                else
                    tmp.dropboxpath = '/home/kimyy/Dropbox';
                end
                addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
                [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
                addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
                    tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
                    'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));

                tmp.shadlev = [0 35];
                tmp.rms_shadlev = [0 4];
                tmp.trendlev = [-10 10];  %% trend lev
                tmp.abstrendlev =[4 7];
                tmp.reltrendlev =[-5 5];
                tmp.conlev  = 0:5:35;
                tmp.meanplotlev =[-0.3 0.3];
                tmp.trendplotlev = [3 7];
                tmp.sshlev =[-0.7 1.3];
                tmp.sshdifflev = [40 70];

                % for snu_desktopd
%                 testname=GCM_info.name{testnameind2};    % % need to change
                GCM_info.abb=GCM_info.abbs{testnameind2};
        %         inputyear1 = [2020]; % % put year which you want to plot [year year ...]
        %         inputyear1 = [1993:2014]; % % put year which you want to plot [year year ...]
        %         inputyear1 = [1995:2014]; % % put year which you want to plot [year year ...]
        %         inputyear1 = [1985:2014]; % % put year which you want to plot [year year ...]
        %         inputyear1 = [1993:1993]; % % put year which you want to plot [year year ...]
        %         inputyear1 = [2015:2015]; % % put year which you want to plot [year year ...]
               GCM_info.years = years_group(years_groupi);  

                season=seasons_group{seasons_groupi};
                switch(season)
                    case 'all'
                        RCM_info.months =1:12;  
                    case 'summer'
                        RCM_info.months =8;  
                    case 'winter'
                        RCM_info.months =2;  
                end
                
                if max(inputyear1)<=2014
                    GCM_info.scenname ='historical';
                else
                    GCM_info.scenname ='ssp585';
                end

                tmp.regionname=GCM_info.region{regionind2};
                [GCM_grid.refpolygon, GCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

                [cmaps.bwrmap, tmp.error_status] = Func_0009_get_colormaps('bwr', tmp.dropboxpath);
                cmaps.wrmap = cmaps.bwrmap(51:100,:);
                [cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);

                % % for windows
%                 figrawdir =strcat('D:\MEPL\project\SSH\6th_year\figure\cmip6\',testname,filesep,regionname,filesep); % % where figure files will be saved
                tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\new_fig_param_kyy_CMIP6_', regionname, '.m'];
    %             dirs.cmip6dir = strcat('G:\Data\Model\cmip6\'); % % where data files are
                dirs.cmip6dir = strcat('D:\Data\Model\CMIP6\NWP\'); % % where data files are
                dirs.figrawdir =strcat('D:\MEPL\project\SSH\6th_year\figure\nwp_1_20\'); % % where figure files will be saved
                
%                 tmp.testname=testname;
%                 RCM_info.years=inputyear1;

                run(tmp.param_script);
                
                
        % start-------------------- earlier decadal current plot
                pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_uv_',num2str(min(inputyear1),'%04i'), ...
                    '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
                dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'vec', tmp.fs, ...
                    num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
                if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
                    mkdir(strcat(dirs.figdir));
                end 
                tmp.tifname=strcat(dirs.figdir, tmp.testname, '_clim_uv_',num2str(min(RCM_info.years),'%04i'), ...
                    '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg
        %         if (exist(pngname , 'file') ~= 2)        
                    for yearij=1:length(inputyear1)
                        tempyear=inputyear1(yearij);
                        yearstr=num2str(tempyear, '%04i');
                        for monthij=1:length(inputmonth)
                            tempmonth=inputmonth(monthij);
                            monthstr=num2str(tempmonth, '%02i');
                            
                            varname='uo';
                            filedir = strcat(dirs.cmip6dir, varname, filesep, GCM_info.scenname, '\Omon\', testname, filesep); % % where data files are
                            ufilename=[filedir, filesep, varname, '_Omon_', GCM_info.scenname, '_', testname, '_', yearstr, '.nc'];
                            tind_u=tempmonth;
                            
                            varname='vo';
                            filedir = strcat(dirs.cmip6dir, varname, filesep, GCM_info.scenname, '\Omon\', testname, filesep); % % where data files are
                            vfilename=[filedir, filesep, varname, '_Omon_', GCM_info.scenname, '_', testname, '_', yearstr, '.nc'];
                            tind_v=tempmonth;
                            if (exist('lon_rho' , 'var') ~= 1)
                                lonfilename=[filedir, filesep, 'lon_Omon_', GCM_info.scenname, '_', testname, '.nc'];
                                latfilename=[filedir, filesep, 'lat_Omon_', GCM_info.scenname, '_', testname, '.nc'];
                                switch(testname)
                                    case{'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
                                        lonname='lon';
                                        latname='lat';
                                        xname = 'x';
                                        yname = 'y';
                                    case{'EC-Earth3-Veg', 'ACCESS-CM2', 'CMCC-ESM2'}
                                        lonname='longitude';
                                        latname='latitude';
                                        xname = 'i';
                                        yname = 'j';
                                end
                                lon_rho=ncread(lonfilename, lonname);
                                lat_rho=ncread(latfilename, latname);
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
                            
                            data_info = ncinfo(ufilename, 'uo'); 
                            u = ncread(ufilename,'uo',[lon_min(1) lat_min(1) 1 tind_u], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                            v = ncread(vfilename,'vo',[lon_min(1) lat_min(1) 1 tind_v], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
        
                            if (exist('mean_u' , 'var') ~= 1)
                                mean_u=zeros(size(u));
                                mean_v=zeros(size(v));
                            end
                            mean_u=mean_u + (u / (length(inputyear1) * length(inputmonth)));
                            mean_v=mean_v + (v / (length(inputyear1) * length(inputmonth)));
                        end
                    end
        %             u_rho = u2rho_2d(mean_u')';
        %             v_rho = v2rho_2d(mean_v')';
                    u_rho = mean_u;
                    v_rho = mean_v;
                    
                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho=u_rho.*mask_model;
                    v_rho=v_rho.*mask_model;
                    
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval;
        %                 ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
        %                 ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2));
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
        
                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        
                    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
        
                    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    if (min(inputyear1) == max(inputyear1))
                        titlename = strcat('UV mean, ',abbname,',(',num2str(min(inputyear1),'%04i'),') ');  %% + glacier contribution
                    else
                        titlename = strcat('UV mean, ',abbname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
                    end
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
        
                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
                    saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
                    close all;
                    clear lon_rho mean_u ref_vec_x_range
        %         end
        % end-------------------- earlier decadal current plot


        % start-------------------- earlier decadal SST, SSS, SSH plot
                    for varind2=1:length(GCM_info.vars)
                        variable=GCM_info.vars{varind2};
                        pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', variable,'_',num2str(min(inputyear1),'%04i'), ...
                            '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg

                        dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, variable, tmp.fs, ...
                            num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
                        if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
                            mkdir(strcat(dirs.figdir));
                        end 
                        tmp.tifname=strcat(dirs.figdir, tmp.testname, '_clim_', variable, '_',num2str(min(RCM_info.years),'%04i'), ...
                             '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg

        %                 if (exist(pngname , 'file') ~= 2)        
                            run(param_script);
                            for yearij=1:length(inputyear1)
                                tempyear=inputyear1(yearij);
                                yearstr=num2str(tempyear, '%04i');
                                for monthij=1:length(inputmonth)
                                    tempmonth=inputmonth(monthij);
                                    monthstr=num2str(tempmonth, '%02i');
                                    filedir = strcat(dirs.cmip6dir, varname, filesep, GCM_info.scenname, '\Omon\', testname, filesep); % % where data files are
                                    filename=[filedir, filesep, varname, '_Omon_', GCM_info.scenname, '_', testname, '_', yearstr, '.nc'];                                                    
                                    tind=tempmonth;

                                    if (exist('lon_rho' , 'var') ~= 1)
                                        lonfilename=[filedir, filesep, 'lon_Omon_', GCM_info.scenname, '_', testname, '.nc'];
                                        latfilename=[filedir, filesep, 'lat_Omon_', GCM_info.scenname, '_', testname, '.nc'];
                                        switch(testname)
                                            case{'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
                                                lonname='lon';
                                                latname='lat';
                                            case{'EC-Earth3-Veg', 'ACCESS-CM2', 'CMCC-ESM2'}
                                                lonname='longitude';
                                                latname='latitude';
                                        end
                                        lon_rho=ncread(lonfilename, lonname);
                                        lat_rho=ncread(latfilename, latname);
                                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                    end
                                    data_info = ncinfo(filename, varname); 

                                    if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
        % %                                 zos = ncread(vfilename,'zos',[lon_min(1) lat_min(1) 1 tind_v], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
        % %                                 zosmask=zeros(size(zos));
        % %                                 zosmask(isfinite(zos))=1;
                                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                                    elseif (strcmp(variable,'SSH')==1)
                                        data = ncread(filename,varname,[lon_min(1) lat_min(1) tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                                    end
                                    if (exist('mean_data' , 'var') ~= 1)
                                        mean_data=zeros(size(data));
                                    end
                                    mean_data=mean_data + (data / (length(inputyear1) * length(inputmonth)));
                                end
                            end
                            if(strcmp(variable, 'SST'))
                                mean_data=mean_data;
                            elseif (strcmp(variable,'SSH')==1)
        %                         mean_data=mean_data-mean(mean_data(:),'omitnan');
                                load(['D:\Data\Model\CMIP6\zos_correction\GCM_corr_', GCM_info.scenname, '.mat']);
                                if strcmp(GCM_info.scenname, 'historical')
                                    GCM_correction=GCM_corr(testnameind2).ts((min(inputyear1)-1985)*12+1:(max(inputyear1)-1985)*12+12);
                                else
                                    GCM_correction=GCM_corr(testnameind2).ts((min(inputyear1)-2015)*12+1:(max(inputyear1)-2015)*12+12);
                                end
                                GCM_corr_reshape=reshape(GCM_correction, [12, length(GCM_correction)/12]);
                                GCM_corr_tgt=GCM_corr_reshape(inputmonth,:);
                                GCM_corr_mean=mean(GCM_corr_tgt(:));
                                mean_data=mean_data + GCM_corr_mean;
        %                         mean(mean_data(:), 'omitnan')
                                [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho)
        %                         -5.5537; -5.4880; 0.5504; -5.5819; -5.8024; mean
        %                         ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor(testname);
                                ssh_correction_for_fig = Func_0017_SSH_correction_for_CMIP6_RMSE(testname);
                                mean_data=mean_data-ssh_correction_for_fig;
                            end
                            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                            mask_model(mask_model==0)=NaN;
                            mean_data=mean_data.*mask_model;
                            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                            hold on;

                            [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho)
                            m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                            shading(gca,m_pcolor_shading_method);   

                            if(strcmp(variable, 'SST'))
                                [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                                        'linewidth', 1.5, 'linestyle', '-');
                                    clabel(C,h2,'FontSize',13,'Color','k', ...
                                        'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                            end

                            m_gshhs_i('color',m_gshhs_line_color)  
                            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                            if (min(inputyear1) == max(inputyear1))
                                titlename = strcat(variable, ' mean, ',abbname,',(',num2str(min(inputyear1),'%04i'),') ');  %% + glacier contribution
                            else
                                titlename = strcat(variable, ' mean, ',abbname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
                            end
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
                            saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
                            saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
                            close all;
                            clear lon_rho mean_data
        %                 end
                    end
        % end-------------------- earlier decadal SST plot
        % % % % 
        % 

% % %         % start-------------------- earlier YSBCW plot
% % %                 if (strcmp(regionname, 'YS')==1)
% % %                     pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', 'YSBCW','_',num2str(min(inputyear1),'%04i'), ...
% % %                         '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
% % %                     if (exist(pngname , 'file') ~= 2)       
% % %                         variable='BT';
% % %                         run(param_script);
% % %                         for yearij=1:length(inputyear1)
% % %                             tempyear=inputyear1(yearij);
% % %                             yearstr=num2str(tempyear, '%04i');
% % %                             tempmonth=8;
% % %                             monthstr=num2str(tempmonth, '%02i');
% % % 
% % %                             filedir = strcat(dirs.cmip6dir, varname, filesep, GCM_info.scenname, '\Omon\', testname, filesep); % % where data files are
% % %                             filename=[filedir, filesep, varname, '_Omon_', GCM_info.scenname, '_', testname, '_', yearstr, '.nc'];                                                    
% % %                             tind=tempmonth;
% % % 
% % %                              if (exist('lon_rho' , 'var') ~= 1)
% % %                                 lonfilename=[filedir, filesep, 'lon_Omon_', GCM_info.scenname, '_', testname, '.nc'];
% % %                                 latfilename=[filedir, filesep, 'lat_Omon_', GCM_info.scenname, '_', testname, '.nc'];
% % %                                 switch(testname)
% % %                                     case{'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
% % %                                         lonname='lon';
% % %                                         latname='lat';
% % %                                         xname = 'x';
% % %                                         yname = 'y';
% % %                                     case{'EC-Earth3-Veg', 'ACCESS-CM2', 'CMCC-ESM2'}
% % %                                         lonname='longitude';
% % %                                         latname='latitude';
% % %                                         xname = 'i';
% % %                                         yname = 'j';
% % %                                 end
% % %                                 lon_rho=ncread(lonfilename, lonname);
% % %                                 lat_rho=ncread(latfilename, latname);
% % %                                 [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
% % %                                 cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % %                                 cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % %                             end
% % %                             data_info = ncinfo(filename, varname); 
% % % 
% % %                             data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
% % %                             numx=size(data,1);
% % %                             numy=size(data,2);
% % %                             for lonind=1:numx
% % %                                 for latind=1:numy
% % %                                     vert_temp=squeeze(data(lonind,latind,:,1));
% % %                                     nanind=find(isnan(vert_temp));
% % %                                     botind=min(nanind)-1;
% % %                                     if (botind>0)
% % %                                         bottomtemp(lonind,latind)=vert_temp(botind);
% % %                                     else
% % %                                         bottomtemp(lonind,latind)=NaN;
% % %                                     end
% % %                                     clear nanind botind
% % %                                 end
% % %                             end
% % %                             if (exist('mean_data' , 'var') ~= 1)
% % %                                 mean_data=zeros(size(bottomtemp));
% % %                             end
% % %                             mean_data=mean_data + (bottomtemp / length(inputyear1));
% % %                         end
% % % 
% % %                         m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % %                         hold on;
% % % 
% % % 
% % %                         m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
% % %                         shading(gca,m_pcolor_shading_method);   
% % % 
% % %                         [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [6 8 10], m_contour_color, 'linewidth', m_contour_linewidth);
% % %                         clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
% % %                             'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
% % % 
% % %                         m_gshhs_i('color',m_gshhs_line_color)  
% % %                         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % 
% % %                         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % %                         titlename = strcat('YSBCW-Aug', ' mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
% % %                         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % 
% % %                         % set colorbar 
% % %                         h = colorbar;
% % %                         colormap(jet);
% % %                         set(h,'fontsize',colorbar_fontsize);
% % %                         title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
% % %                         caxis(colorbar_lev);
% % % 
% % %                         set(gcf, 'PaperUnits', 'points');
% % %                         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % %                         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % %                         saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
% % %                         close all;
% % %                         clear lon_rho mean_data
% % %                     end
% % %                 end
% % %         % end-------------------- earlier YSBCW plot
            end
        end
    end
end