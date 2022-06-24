% %  Updated 10-Oct-2021 by Yong-Yub Kim, structure
% %  Updated 10-Nov-2021 by Yong-Yub Kim, colorbar

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
        '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.tif'); %% ~_year_month.jpg
    if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(3)==2)      
        run(tmp.param_script);
        tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
        '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
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
                            tmp.uwindfilename=[dirs.filedir, 'Uwind', tmp.fs, tmp.yearstr, filesep, ...
                                'NWP_pck_', tmp.testname, '_', 'Uwind', '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                            tmp.vwindfilename=[dirs.filedir, 'Vwind', tmp.fs, tmp.yearstr, filesep, ...
                                'NWP_pck_', tmp.testname, '_', 'Vwind', '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                            tmp.sustrfilename=[dirs.filedir, 'sustr', tmp.fs, tmp.yearstr, filesep, ...
                                'NWP_pck_', tmp.testname, '_', 'sustr', '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                            tmp.svstrfilename=[dirs.filedir, 'svstr', tmp.fs, tmp.yearstr, filesep, ...
                                'NWP_pck_', tmp.testname, '_', 'svstr', '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
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
                        cut_lon_psi = ...
                            RCM_grid.lon_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1)-1, RCM_grid.lat_min(1):RCM_grid.lat_max(1)-1);
                        cut_lat_psi = ...
                            RCM_grid.lat_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1)-1, RCM_grid.lat_min(1):RCM_grid.lat_max(1)-1);
                        RCM_grid.xdist=1./RCM_grid.pm;
                        RCM_grid.ydist=1./RCM_grid.pn;
                        cut_xdist= ...
                            RCM_grid.xdist(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                        cut_ydist= ...
                            RCM_grid.ydist(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                    end
%                     tmp.data_info = ncinfo(tmp.filename, param.varname); 

                    if (strcmp(tmp.variable,'SST')==1 || strcmp(tmp.variable,'SSS')==1)
                        tmp.data = ncread(tmp.filename,param.varname,[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    elseif (strcmp(tmp.variable,'SSH')==1 || strcmp(tmp.variable,'Uwind')==1 || strcmp(tmp.variable,'Vwind')==1 || strcmp(tmp.variable,'shflux')==1)
                        tmp.data = ncread(tmp.filename,param.varname,[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                    elseif (strcmp(tmp.variable,'u')==1)
                       tmp.data = ncread(tmp.filename,param.varname,[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1) RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                       tmp.data = u2rho_2d(tmp.data')';
                    elseif (strcmp(tmp.variable,'v')==1)
                       tmp.data = ncread(tmp.filename,param.varname,[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1) 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                       tmp.data = v2rho_2d(tmp.data')';
                   elseif (strcmp(tmp.variable,'wcurl')==1)
                        tmp.uwind = ncread(tmp.uwindfilename,'Uwind',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        tmp.vwind = ncread(tmp.vwindfilename,'Vwind',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        tmp.data=get_curl(tmp.uwind, tmp.vwind, cut_xdist, cut_ydist);
                    elseif (strcmp(tmp.variable,'wstrcurl')==1)
                        tmp.sustr = ncread(tmp.sustrfilename,'sustr',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1) RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        tmp.svstr = ncread(tmp.svstrfilename,'svstr',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        tmp.wstrcurl=get_curl(tmp.sustr, tmp.svstr, cut_xdist, cut_ydist);
                        tmp.data=griddata(cut_lon_psi,cut_lat_psi,tmp.wstrcurl,cut_lon_rho,cut_lat_rho);
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

        [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho);

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
            tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
        else
            tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
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
        switch tmp.variable
            case {'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'u', 'v', }
                caxis(param.colorbar_lev);
            case {'wstrcurl', 'wcurl'}
%                 caxis([-max(abs(mean_data(:))), max(abs(mean_data(:)))]);
                caxis(param.colorbar_lev);
                colormap(cmaps.byrmap);
        end

        disp(['M = ', num2str(tmp.m_value)]);
        m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        clear mean_data
        if (isfield(RCM_grid, 'lon_rho') == 1)
            RCM_grid=rmfield(RCM_grid, 'lon_rho');
        end
    end
end

function wstrcurl=get_curl(uwstr, vwstr, xdist, ydist)
    loncount=size(uwstr,1);
    latcount=size(uwstr,2);
    if (size(uwstr,1) == size(vwstr,1))
        wstrcurl=NaN(loncount, latcount);
        for i=1:loncount-2
            for j=1:latcount-2
                wstrcurl(i+1,j+1) = (vwstr(i+2,j+1)-vwstr(i,j+1)) ...
                    / (xdist(i+1, j+1)*2) - ...
                    (uwstr(i+1,j+2)-uwstr(i+1,j)) ...
                    / (ydist(i+1, j+1)*2);
            end
        end
        wstrcurl(1,2:latcount-2)=wstrcurl(2,2:latcount-2);
        wstrcurl(loncount, 2:latcount-2) =wstrcurl(loncount-1,2:latcount-2);
        wstrcurl(2:loncount-2,1)=wstrcurl(2:loncount-2,2);
        wstrcurl(2:loncount-2,latcount)=wstrcurl(2:loncount-2,latcount-1);
        wstrcurl(1,1)=wstrcurl(2,2);
        wstrcurl(1,latcount)=wstrcurl(2,latcount-1);
        wstrcurl(loncount,1)=wstrcurl(loncount-1,2);
        wstrcurl(loncount,latcount)=wstrcurl(loncount-1,latcount-1);
    elseif (size(uwstr,1) ~= size(vwstr,1))
        uloncount=size(uwstr,1);
        ulatcount=size(uwstr,2);
        vloncount=size(vwstr,1);
        vlatcount=size(vwstr,2);
        wstrcurl=NaN(uloncount, vlatcount);
        for i=1:uloncount-1
            for j=1:vlatcount-1
                wstrcurl(i,j) = (vwstr(i+1,j)-vwstr(i,j)) ...
                    / (xdist(i, j)*2) - ...
                    (uwstr(i,j+1)-uwstr(i,j)) ...
                    / (ydist(i, j)*2);  % N m-2 m-1
            end
        end
    end
end