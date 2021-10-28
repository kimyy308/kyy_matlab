% %  Updated 10-Oct-2021 by Yong-Yub Kim, structure

% start-------------------- earlier decadal SST, SSS plot
            for varind2=1:length(RCM_info.vars)
                tmp.variable=RCM_info.vars{varind2};
                dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
                    num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
                if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
                    mkdir(strcat(dirs.figdir));
                end 
            
                tmp.tifname=strcat(dirs.figdir, tmp.testname, '_clim_', tmp.variable, '_', ...
                    num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'),'.tif'); %% ~_year_month.jpg
                if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(3)==2)      
                    run(tmp.param_script);
                    tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
                        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'), '.mat'];
                    if (exist(tmp.matname , 'file') ~= 2 || flags.fig_switch(3)==2)
                        for yearij=1:length(RCM_info.years)
                            tmp.tempyear=RCM_info.years(yearij);
                            tmp.yearstr=num2str(tmp.tempyear, '%04i');
                            for monthij=1:length(RCM_info.months)
                                tmp.tempmonth=RCM_info.months(monthij);
                                tmp.monthstr=num2str(tmp.tempmonth, '%02i');
                                
                                switch tmp.testname
                                    case {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
                                        tmp.filename=[dirs.filedir, param.varname, tmp.fs, tmp.yearstr, filesep, ...
                                            'pck_', tmp.testname, '_', param.varname, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                                    otherwise
                                        tmp.filename=[dirs.filedir, param.varname, tmp.fs, tmp.yearstr, filesep, ...
                                            'NWP_pck_', tmp.testname, '_', param.varname, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                                end
                                if (isfield(RCM_grid, 'lon_rho') ~= 1)
                                    for gridi=1:length(RCM_grid.gridname)
                                        RCM_grid.(RCM_grid.gridname{gridi})=ncread(RCM_grid.(['filename_', RCM_grid.gridname{gridi}]), RCM_grid.gridname{gridi});
                                    end

                                    [RCM_grid.lon_min, RCM_grid.lon_max, RCM_grid.lat_min, RCM_grid.lat_max] = ...
                                        findind_Y(1/20, RCM_grid.domain(1:4), RCM_grid.lon_rho, RCM_grid.lat_rho);
                                    cut_lon_rho = ...
                                        RCM_grid.lon_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                                    cut_lat_rho = ...
                                        RCM_grid.lat_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                                end
                                tmp.data_info = ncinfo(tmp.filename, param.varname); 

                                if (strcmp(tmp.variable,'SST')==1 || strcmp(tmp.variable,'SSS')==1)
                                    tmp.data = ncread(tmp.filename,param.varname,[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                                        [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                                elseif (strcmp(tmp.variable,'SSH')==1 || strcmp(tmp.variable,'Uwind')==1 || strcmp(tmp.variable,'Vwind')==1)
                                    tmp.data = ncread(tmp.filename,param.varname,[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                                        [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%                                 elseif (strcmp(tmp.variable,'BT')==1)
%                                     tmp.data = ncread(tmp.filename,param.varname,[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                                end
                                if (exist('mean_data', 'var') ~= 1)
                                    mean_data=zeros(size(tmp.data));
                                end
                                mean_data=mean_data + ...
                                    (tmp.data / (length(RCM_info.years) * length(RCM_info.months)));
                            end
                        end
                        save(tmp.matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
                    else
                        load(tmp.matname);
                    end
                    
%                     [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho)
                    
                    if(strcmp(tmp.variable, 'SST'))
                        mean_data=mean_data;
                    elseif (strcmp(tmp.variable,'SSH')==1)
%                         ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor(tmp.testname);
                        ssh_correction_for_fig = Func_0017_SSH_correction_for_CMIP6_RMSE(tmp.testname);
                        mean_data=mean_data-ssh_correction_for_fig;
                    end

                    RCM_grid.mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
                    RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
                    mean_data=mean_data.*RCM_grid.mask_model;
                
                    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
                    hold on;

                    [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho);
                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                    shading(gca,param.m_pcolor_shading_method);   
                    if(strcmp(tmp.variable, 'SST'))
%                         [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [2, 4, 6, 8, 10], 'color','k', ...
%                                 'linewidth', 1.5, 'linestyle', '-');
                        [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [5, 10, 15, 20, 25, 30], 'color','k', ...
                                'linewidth', 1.5, 'linestyle', '-');
                            clabel(C,h2,'FontSize',13,'Color','k', ...
                                'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                    end
                    
                    m_gshhs_i('color',param.m_gshhs_line_color)  
                    m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
                    if min(RCM_info.years) == max(RCM_info.years)
                        tmp.titlename = strcat(tmp.variable, ', ',tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
                    else
                        tmp.titlename = strcat(tmp.variable, ', ',tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
                    end
                    title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
%                     if (strcmp(tmp.variable,'SSH')==1)
%                         colormap(flip(cool));
%                     else
                        colormap(jet);
%                     end

                    set(h,'fontsize',param.colorbar_fontsize);
                    title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
                    caxis(param.colorbar_lev);
                    
                    disp(['M = ', num2str(tmp.m_value)]);
                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
                    set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
                    saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
                    close all;
                    clear mean_data
                    RCM_grid=rmfield(RCM_grid, 'lon_rho');
                end
            end

% % % % start-------------------- later decadal SST, SSS plot
% % %         fig_flag=fig_flags{4,2};
% % %         while (fig_flag)
% % %             for varind2=1:length(RCM_info.vars)
% % %                 variable=RCM_info.vars{varind2};
% % %                 tifname=strcat(outfile, '_', testname,'_',tmp.regionname, '_clim_', variable,'_',num2str(min(inputyear2),'%04i'), ...
% % %                     '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
% % %                 if (exist(tifname , 'file') ~= 2 || fig_flag==2)      
% % %                     run(param_script);
% % %                     matname = [matdir, testname, '_', tmp.regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
% % %                     if (exist(matname , 'file') ~= 2 || fig_flag==2)
% % %                         for yearij=1:length(inputyear2)
% % %                             tmp.tempyear=inputyear2(yearij);
% % %                             tmp.yearstr=num2str(tmp.tempyear, '%04i');
% % %                             for monthij=1:length(RCM_info.months)
% % %                                 tmp.tempmonth=RCM_info.months(monthij);
% % %                                 tmp.monthstr=num2str(tmp.tempmonth, '%02i');
% % %                                 filename=[filedir, tmp.yearstr, filesep,'pck_', testname, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
% % % 
% % %                                 if (exist('lon_rho' , 'var') ~= 1)
% % %                                     lon_rho=ncread(filename, 'lon_rho');
% % %                                     lat_rho=ncread(filename, 'lat_rho');
% % %                                     [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
% % %                                     cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % %                                     cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % %                                 end
% % %                                 data_info = ncinfo(filename, varname); 
% % % 
% % %                                 if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
% % %                                     data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
% % %                                 elseif (strcmp(variable,'SSH')==1)
% % %                                     data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
% % %                                 elseif (strcmp(variable,'BT')==1)
% % %                                     data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
% % %                                 end
% % %                                 if (exist('mean_data' , 'var') ~= 1)
% % %                                     mean_data=zeros(size(data));
% % %                                 end
% % %                                 mean_data=mean_data + (data / (length(inputyear2) * length(RCM_info.months)));
% % %                             end
% % %                         end
% % %                         save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
% % %                     else
% % %                         load(matname);
% % %                     end
% % %                     mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
% % %                     mask_model(mask_model==0)=NaN;
% % %                     mean_data=mean_data.*mask_model;
% % %                     
% % %                     m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % %                     hold on;
% % % 
% % %                     m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
% % %                     shading(gca,m_pcolor_shading_method);   
% % % 
% % %                     m_gshhs_i('color',m_gshhs_line_color)  
% % %                     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % 
% % %                     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % %                     titlename = strcat(variable, ' mean, ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
% % %                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % 
% % %                     % set colorbar 
% % %                     h = colorbar;
% % %                     colormap(jet);
% % %                     set(h,'fontsize',colorbar_fontsize);
% % %                     title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
% % %                     caxis(colorbar_lev);
% % %                     
% % %                     disp(['M = ', num2str(mean(mean_data(:),'omitnan'))]);
% % %                     set(gcf, 'PaperUnits', 'points');
% % %                     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % %                     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % %                     saveas(gcf,tifname,'tif');
% % %                     RemoveWhiteSpace([], 'file', tifname);
% % %                     close all;
% % %                     clear lon_rho mean_data
% % %                 end
% % %             end
% % %             fig_flag=0;
% % %         end