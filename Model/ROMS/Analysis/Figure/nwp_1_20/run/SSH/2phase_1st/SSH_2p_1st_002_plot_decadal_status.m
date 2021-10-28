% %  Updated 09-Oct-2021 by Yong-Yub Kim, structure


close all; clear all;  clc;
warning off;

% if(isempty(gcp('nocreate')))
%     parpool(4);
% end

% % % configuration of RCM
RCM_info.name={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
% RCM_info.name={'test2107', 'test2108', 'test2109', 'test2110', 'test2111'};
% RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
RCM_info.abbs = {'RCM-CNRM', 'RCM-EC-Veg', 'RCM-ACC', 'RCM-CNRM-HR', 'RCM-CMCC'};
RCM_info.model = 'nwp_1_20';
% RCM_info.dataroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';
RCM_info.region = {'AKP4'};
% RCM_info.region = {'NWP', 'AKP4'};
% RCM_info.vars = {'SST', 'SSH', 'SSS', 'Uwind', 'Vwind'};
% RCM_info.vars = {'Uwind', 'Vwind'};

RCM_info.vars = {'SSH'};
% RCM_info.years = 1985:2014;  
RCM_info.years = 1993:2014;  
% RCM_info.years = 1995:2014;  

RCM_info.months =1:12;  
RCM_grid.dl = 1/20;
RCM_grid.gridname = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v'};

% RCM_info.testnames = {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};


% all_region2 ={'NWP', 'YS', 'AKP2'}

% all_region2 ={'NWP'};

% RCM_info.vars = {'SST', 'SSH', 'SSS'};
% RCM_info.vars = {'SSH'};

% all_region2 ={'NWP'}
for testnameind2=1:length(RCM_info.name)
    for regionind2=1:length(RCM_info.region)
        close all;
        clearvars '*' -except RCM_info RCM_grid testnameind2 regionind2
        tmp.fs=filesep;
        % % % 
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

        % for snu_desktop
        tmp.testname=RCM_info.name{testnameind2};   % % need to change
        tmp.abb=RCM_info.abbs{testnameind2};    % % need to change

%         RCM_info.years = [1985:2014]; % % put year which you want to plot [year year ...]
%           RCM_info.months = [1:12]; % % put month which you want to plot [month month ...]

        flags.fig_name{1}='earlier decadal current plot';
        flags.fig_name{2}='later decadal current plot';
        flags.fig_name{3}='earlier decadal SST, SSS plot';
        flags.fig_name{4}='later decadal SST, SSS plot';
        flags.fig_name{5}='earlier decadal YSBCW plot';
        flags.fig_name{6}='later decadal YSBCW plot';
        
        for flagi=1:6
            fig_flags{flagi,2}=0;
        end
        flags.fig_switch(1)=0;  %1 or 2
        flags.fig_switch(2)=0; 
        flags.fig_switch(3)=2;
        flags.fig_switch(4)=0;
        flags.fig_switch(5)=0;
        flags.fig_switch(6)=0;
        
        tmp.variable ='zeta';
%         run('nwp_polygon_point.m');
        tmp.regionname=RCM_info.region{regionind2};
        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);
       
        [cmaps.bwrmap, tmp.error_status] = Func_0009_get_colormaps('bwr', tmp.dropboxpath);
        cmaps.wrmap = cmaps.bwrmap(51:100,:);
        [cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);
        
        dirs.figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\'); % % where figure files will be saved
        tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param2_kyy_', tmp.regionname, '.m'];
        dirs.filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\backup_surf\', tmp.testname, '\run\'); % % where data files are          
        dirs.matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\mean\');
        dirs.griddir = strcat('D:\Data\Model\ROMS\nwp_1_20\backup_surf\'); % % where grid data files are            
        
        
        for gridi=1:length(RCM_grid.gridname)
            RCM_grid.(['filename_', RCM_grid.gridname{gridi}])=[dirs.griddir, 'NWP_pck_ocean_', RCM_grid.gridname{gridi}, '_NWP.nc'];
        end
        
        run(tmp.param_script);
        
        if (exist(strcat(dirs.matdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.matdir));
        end 

        if flags.fig_switch(1) > 0
            flags.fig_tmp = flags.fig_switch(1);
            SSH_2p_1st_002_sub_001_current_plot;
        end
        
        if flags.fig_switch(3) > 0
            flags.fig_tmp = flags.fig_switch(1);
            SSH_2p_1st_002_sub_003_surface_variable_plot;
        end



% start-------------------- earlier YSBCW plot
        fig_flag=fig_flags{5,2};
        while (fig_flag)
            if (strcmp(regionname, 'YS')==1)
                tifname=strcat(outfile, '_', testname,'_',regionname, '_clim_', 'YSBCW','_',num2str(min(RCM_info.years),'%04i'), ...
                    '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(tifname , 'file') ~= 2 || fig_flag==2)      
                    variable='BT';
                    run(param_script);
                    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'), '.mat'];
                    if (exist(matname , 'file') ~= 2 || fig_flag==2)
                        for yearij=1:length(RCM_info.years)
                            tmp.tempyear=RCM_info.years(yearij);
                            tmp.yearstr=num2str(tmp.tempyear, '%04i');
                            tmp.tempmonth=8;
                            tmp.monthstr=num2str(tmp.tempmonth, '%02i');
                            filename=[filedir, tmp.yearstr, filesep,'pck_', testname, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];

                            if (exist('lon_rho' , 'var') ~= 1)
                                lon_rho=ncread(filename, 'lon_rho');
                                lat_rho=ncread(filename, 'lat_rho');
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
                            data_info = ncinfo(filename, varname); 

                            if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
                                data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                            elseif (strcmp(variable,'SSH')==1)
                                data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                            elseif (strcmp(variable,'BT')==1)
                                data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                            end
                            if (exist('mean_data' , 'var') ~= 1)
                                mean_data=zeros(size(data));
                            end
                            mean_data=mean_data + (data / length(RCM_info.years));
                        end
                        save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
                    else
                        load(matname);
                    end

                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;

                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                    shading(gca,m_pcolor_shading_method);   

                    [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [6 8 10], m_contour_color, 'linewidth', m_contour_linewidth);
                    clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                        'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);

                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat('YSBCW-Aug', ' mean, ',testname,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');  %% + glacier contribution
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
                    saveas(gcf,tifname,'tif');
                    RemoveWhiteSpace([], 'file', tifname);
                    close all;
                    clear lon_rho mean_data
                end
                
            end
            fig_flag=0;
        end
       
% start-------------------- later YSBCW plot
        fig_flag=fig_flags{6,2};
        while (fig_flag)
            if (strcmp(regionname, 'YS')==1)
                tifname=strcat(outfile, '_', testname,'_',regionname, '_clim_', 'YSBCW','_',num2str(min(inputyear2),'%04i'), ...
                    '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(tifname , 'file') ~= 2 || fig_flag==2)       
                    var='BT';
                    run(param_script);
                    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
                    if (exist(matname , 'file') ~= 2|| fig_flag==2)
                        for yearij=1:length(inputyear2)
                            tmp.tempyear=inputyear2(yearij);
                            tmp.yearstr=num2str(tmp.tempyear, '%04i');
                            tmp.tempmonth=8;
                            tmp.monthstr=num2str(tmp.tempmonth, '%02i');
                            filename=[filedir, tmp.yearstr, filesep,'pck_', testname, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];

                            if (exist('lon_rho' , 'var') ~= 1)
                                lon_rho=ncread(filename, 'lon_rho');
                                lat_rho=ncread(filename, 'lat_rho');
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
                            data_info = ncinfo(filename, varname); 

                            if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
                                data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                            elseif (strcmp(variable,'SSH')==1)
                                data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                            elseif (strcmp(variable,'BT')==1)
                                data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                            end
                            if (exist('mean_data' , 'var') ~= 1)
                                mean_data=zeros(size(data));
                            end
                            mean_data=mean_data + (data / length(inputyear2));
                        end
                        save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
                    else
                        load(matname);
                    end

                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;

                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                    shading(gca,m_pcolor_shading_method);   

                    [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [6 8 10], m_contour_color, 'linewidth', m_contour_linewidth);
                    clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                        'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);

                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat('YSBCW-Aug', ' mean, ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
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
                    saveas(gcf,tifname,'tif');
                    RemoveWhiteSpace([], 'file', tifname);
                    close all;
                    clear lon_rho mean_data
                end
            end
            fig_flag=0;
        end

% % start-------------------- earlier decadal current + SST plot
%         tifname=strcat(outfile, '_', testname,'_',regionname, '_clim_uv_SST_',num2str(min(RCM_info.years),'%04i'), ...
%             '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg
% %         if (exist(tifname , 'file') ~= 2)        
%             for yearij=1:length(RCM_info.years)
%                 tmp.tempyear=RCM_info.years(yearij);
%                 tmp.yearstr=num2str(tmp.tempyear, '%04i');
%                 for monthij=1:length(RCM_info.months)
%                     tmp.tempmonth=RCM_info.months(monthij);
%                     tmp.monthstr=num2str(tmp.tempmonth, '%02i');
%                     filename=[filedir, tmp.yearstr, '\', testname, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
% 
%                     if (exist('lon_rho' , 'var') ~= 1)
%                         lon_rho=ncread(filename, 'lon_rho');
%                         lat_rho=ncread(filename, 'lat_rho');
%                         lon_u=ncread(filename, 'lon_u');
%                         lat_u=ncread(filename, 'lat_u');
%                         lon_v=ncread(filename, 'lon_v');
%                         lat_v=ncread(filename, 'lat_v');
%                         [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
%                         [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(1/20, lonlat(1:4), lon_u, lat_u);
%                         [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(1/20, lonlat(1:4), lon_v, lat_v);
%                         cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                         cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                     end
%                     data_info = ncinfo(filename, 'u'); 
% %                     u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) data_info.Size(3) 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
% %                     v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) data_info.Size(3) 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%                     u = ncread(filename,'u',[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1) lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%                     v = ncread(filename,'v',[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1) 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
% 
%                     if (exist('mean_u' , 'var') ~= 1)
%                         mean_u=zeros(size(u));
%                         mean_v=zeros(size(v));
%                     end
%                     mean_u=mean_u + (u / (length(RCM_info.years) * length(RCM_info.months)));
%                     mean_v=mean_v + (v / (length(RCM_info.years) * length(RCM_info.months)));
%                 end
%             end
%             u_rho = u2rho_2d(mean_u')';
%             v_rho = v2rho_2d(mean_v')';
%             if (exist('tmp.ref_vec_x_range' , 'var') ~= 1)
%                 tmp.ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
%                 tmp.ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
%                 tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(m_quiver_x_interval/2)) : round(tmp.ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
%                 tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(m_quiver_y_interval/2)) : round(tmp.ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
%             end
%             u_rho(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=m_quiver_ref_u_value;
%             v_rho(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=m_quiver_ref_v_value;     
%             
%             var='SST';
%             run(param_script);
%             for yearij=1:length(RCM_info.years)
%                 tmp.tempyear=RCM_info.years(yearij);
%                 tmp.yearstr=num2str(tmp.tempyear, '%04i');
%                 for monthij=1:length(RCM_info.months)
%                     tmp.tempmonth=RCM_info.months(monthij);
%                     tmp.monthstr=num2str(tmp.tempmonth, '%02i');
%                     filename=[filedir, tmp.yearstr, '\', testname, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
% 
%                     if (exist('lon_rho' , 'var') ~= 1)
%                         lon_rho=ncread(filename, 'lon_rho');
%                         lat_rho=ncread(filename, 'lat_rho');
%                         [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
%                         cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                         cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                     end
%                     data_info = ncinfo(filename, varname); 
% 
%                     if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
%                         data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
%                     elseif (strcmp(variable,'SSH')==1)
%                         data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%                     elseif (strcmp(variable,'BT')==1)
%                         data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
%                     end
%                     if (exist('mean_data' , 'var') ~= 1)
%                         mean_data=zeros(size(data));
%                     end
%                     mean_data=mean_data + (data / (length(RCM_info.years) * length(RCM_info.months)));
%                 end
%             end
%             mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%                 mask_model(mask_model==0)=NaN;
%             mean_data=mean_data.*mask_model;
%             u_rho=u_rho.*mask_model;
%             v_rho=v_rho.*mask_model;    
%             
%             m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             hold on;
%             m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
%                 shading(gca,m_pcolor_shading_method);   
%               % set colorbar 
%                 h = colorbar;
%                 colormap(jet);
%                 set(h,'fontsize',colorbar_fontsize);
%                 title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
%                 caxis([10 35]);
%         
%             m_gshhs_i('color',m_gshhs_line_color)  
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%             uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
%                             cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
%                             u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
%                             v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
%                             'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
% 
%             m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% %             titlename = strcat('UV mean, ',testname,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');  %% + glacier contribution
% %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
%                       
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%             saveas(gcf,tifname,'tif');
%             hold off
%             close all;
%             clear lon_rho mean_u tmp.ref_vec_x_range
% %         end
% % end-------------------- earlier decadal current + SST plot
    
    end
end


