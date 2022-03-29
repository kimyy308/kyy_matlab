% %  Updated 22-Feb-2022 by Yong-Yub Kim, make

for varind2=1:length(RCM_info.vert_vars)
%% variable, dir, figname set    
    tmp.variable=RCM_info.vert_vars{varind2};
    dirs.figdir=[dirs.figrawdir,'vertical', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs, ...
        num2str(RCM_info.vert_section(1)), '_', ...
        num2str(RCM_info.vert_section(2)), '_', ...
        num2str(RCM_info.vert_section(3)), '_', ...
        num2str(RCM_info.vert_section(4)), '_', ...
        num2str(RCM_info.vert_section(5)), '_', ...
        num2str(RCM_info.vert_section(6)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 
    tmp.tifname=strcat(dirs.figdir, tmp.testname, '_clim_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
%% loop and calculation    
    if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(1)==2)      
        run(tmp.param_script);
        tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
            '_vert_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
        '_', RCM_info.season, '.mat'];
%         if (exist(tmp.matname , 'file') ~= 2 || flags.fig_switch(1)==2)
            for yearij=1:length(RCM_info.years)
                tmp.tempyear=RCM_info.years(yearij);
                tmp.yearstr=num2str(tmp.tempyear, '%04i');
                for monthij=1:length(RCM_info.months)
                    tmp.tempmonth=RCM_info.months(monthij);
                    tmp.monthstr=num2str(tmp.tempmonth, '%02i');

                    tmp.filename=[dirs.vert_filedir,  tmp.fs, tmp.yearstr, filesep, ...
                        tmp.testname,  '_monthly_std_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                    
                    RCM_vert_grid.vert_section = RCM_info.vert_section;
                    if (isfield(RCM_vert_grid, 'lon_rho') ~= 1)
%                         for gridi=1:length(RCM_vert_grid.gridname)
%                             RCM_vert_grid.(RCM_vert_grid.gridname{gridi})=ncread(RCM_vert_grid.(['filename_', RCM_vert_grid.gridname{gridi}]), RCM_vert_grid.gridname{gridi});
%                         end
                        RCM_vert_grid.lon_rho = ncread(tmp.filename, 'lon_rho');
                        RCM_vert_grid.lat_rho = ncread(tmp.filename, 'lat_rho');
                        RCM_vert_grid.depth = ncread(tmp.filename, 'depth');
                        
                        
                        [RCM_vert_grid.lon_min, RCM_vert_grid.lon_max, RCM_vert_grid.lat_min, RCM_vert_grid.lat_max] = ...
                            findind_Y(1/20, RCM_info.vert_section(1:4), RCM_vert_grid.lon_rho, RCM_vert_grid.lat_rho);
                        RCM_vert_grid.lon_min = RCM_vert_grid.lon_min+1;
                        RCM_vert_grid.lat_min = RCM_vert_grid.lat_min+1;
                        RCM_vert_grid.lon_max = RCM_vert_grid.lon_max-1;
                        RCM_vert_grid.lat_max = RCM_vert_grid.lat_max-1;
                        RCM_vert_grid.lon_count = RCM_vert_grid.lon_max - RCM_vert_grid.lon_min + 1;
                        RCM_vert_grid.lat_count = RCM_vert_grid.lat_max - RCM_vert_grid.lat_min + 1;
                        if (RCM_vert_grid.lon_count > RCM_vert_grid.lat_count)
                            RCM_vert_grid.lonlatflag=1;
                            RCM_vert_grid.lat_max=RCM_vert_grid.lat_min;
                            RCM_vert_grid.lat_count=1;
                        else
                            RCM_vert_grid.lonlatflag=2;
                            RCM_vert_grid.lon_max=RCM_vert_grid.lon_min;
                            RCM_vert_grid.lon_count=1;
                        end
                        
                        tmp.depth=-RCM_vert_grid.depth;
                        tmp.depth=abs(tmp.depth-RCM_info.vert_section(5));
                        RCM_vert_grid.dep_min = find(tmp.depth==min(tmp.depth));
                        tmp.depth=-RCM_vert_grid.depth;
                        tmp.depth=abs(tmp.depth-RCM_info.vert_section(6));
                        RCM_vert_grid.dep_max = find(tmp.depth==min(tmp.depth));
                        RCM_vert_grid.dep_count= RCM_vert_grid.dep_min - RCM_vert_grid.dep_max +1;
                        RCM_vert_grid.thick = diff([0; (RCM_vert_grid.depth(2:end)+RCM_vert_grid.depth(1:end-1))/2; 3750]);
                        
%                         RCM_info.vert_section = [129, 131, 37, 37, -200, 0]; 
                        
                        
                        RCM_vert_grid.cut_depth = ...
                            -RCM_vert_grid.depth(RCM_vert_grid.dep_max : RCM_vert_grid.dep_min);
                        RCM_vert_grid.cut_thick = ...
                            RCM_vert_grid.thick(RCM_vert_grid.dep_max : RCM_vert_grid.dep_min);
                        
                        if RCM_vert_grid.lonlatflag==1
                            RCM_vert_grid.cut_lon_rho = ...
                                RCM_vert_grid.lon_rho(RCM_vert_grid.lon_min(1):RCM_vert_grid.lon_max(1), RCM_vert_grid.lat_min(1):RCM_vert_grid.lat_max(1));
                            RCM_vert_grid.cut_lat_rho = ...
                                RCM_vert_grid.lat_rho(RCM_vert_grid.lon_min(1):RCM_vert_grid.lon_max(1), RCM_vert_grid.lat_min(1):RCM_vert_grid.lat_max(1));
                            RCM_vert_grid.cut_lon_rho2 = repmat(RCM_vert_grid.cut_lon_rho, 1, RCM_vert_grid.dep_count, 1);
                            RCM_vert_grid.cut_lat_rho2 = repmat(RCM_vert_grid.cut_lat_rho, 1, RCM_vert_grid.dep_count, 1);
                            RCM_vert_grid.cut_depth2 = repmat(RCM_vert_grid.cut_depth, 1, RCM_vert_grid.lon_count, 1)';
                            RCM_vert_grid.cut_thick2 = repmat(RCM_vert_grid.cut_thick, 1, RCM_vert_grid.lon_count, 1)';
                            RCM_vert_grid.lonlat = RCM_vert_grid.cut_lon_rho2;
                        elseif RCM_vert_grid.lonlatflag==2
                            RCM_vert_grid.cut_lon_rho = ...
                                RCM_vert_grid.lon_rho(RCM_vert_grid.lon_min(1):RCM_vert_grid.lon_max(1), RCM_vert_grid.lat_min(1):RCM_vert_grid.lat_max(1))';
                            RCM_vert_grid.cut_lat_rho = ...
                                RCM_vert_grid.lat_rho(RCM_vert_grid.lon_min(1):RCM_vert_grid.lon_max(1), RCM_vert_grid.lat_min(1):RCM_vert_grid.lat_max(1))';
                            RCM_vert_grid.cut_lon_rho2 = repmat(RCM_vert_grid.cut_lon_rho, 1, RCM_vert_grid.dep_count, 1);
                            RCM_vert_grid.cut_lat_rho2 = repmat(RCM_vert_grid.cut_lat_rho, 1, RCM_vert_grid.dep_count, 1);
                            RCM_vert_grid.cut_depth2 = repmat(RCM_vert_grid.cut_depth, 1, RCM_vert_grid.lat_count, 1)';
                            RCM_vert_grid.cut_thick2 = repmat(RCM_vert_grid.cut_thick, 1, RCM_vert_grid.lat_count, 1)';
                            RCM_vert_grid.lonlat = RCM_vert_grid.cut_lat_rho2;
                        end
                    end
%                     tmp.data_info = ncinfo(tmp.filename, param.varname); 

                    switch tmp.variable
                        case 'rho'
                            tmp.error_status=1;
                        otherwise
                            tmp.data = ncread(tmp.filename,param.varname,[RCM_vert_grid.lon_min(1) RCM_vert_grid.lat_min(1) RCM_vert_grid.dep_max 1], ...
                            [RCM_vert_grid.lon_max(1)-RCM_vert_grid.lon_min(1)+1 RCM_vert_grid.lat_max(1)-RCM_vert_grid.lat_min(1)+1 ...
                            RCM_vert_grid.dep_min-RCM_vert_grid.dep_max+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    end
%                     tmp.zeta = ncread(tmp.filename,param.varname,[RCM_vert_grid.lon_min(1) RCM_vert_grid.lat_min(1) RCM_vert_grid.dep_max 1], ...
%                             [RCM_vert_grid.lon_max(1)-RCM_vert_grid.lon_min(1)+1 RCM_vert_grid.lat_max(1)-RCM_vert_grid.lat_min(1)+1 ...
%                             RCM_vert_grid.dep_min-RCM_vert_grid.dep_max+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    
                    if (exist('mean_data', 'var') ~= 1)
                        mean_data=zeros(size(tmp.data));
                    end
                    mean_data=mean_data + ...
                        (tmp.data / (length(RCM_info.years) * length(RCM_info.months)));
                end  %% month
            end  %% year
%             save(tmp.matname, 'mean_data', 'RCM_vert_grid', 'RCM_info');
%         else
%             load(tmp.matname, 'mean_data', 'RCM_vert_grid');
%         end  %% mat flag

        [tmp.m_value, tmp.error_status] = ...
            Func_0021_get_vert_area_weighted_mean(mean_data, RCM_vert_grid.cut_lon_rho2, RCM_vert_grid.cut_lat_rho2, ...
            RCM_vert_grid.cut_thick2, RCM_vert_grid.lonlatflag);
        hold on;
        pcolor(RCM_vert_grid.lonlat,RCM_vert_grid.cut_depth2, squeeze(mean_data));
        
        shading(gca,param.m_pcolor_shading_method);   
        if(strcmp(tmp.variable, 'vert_temp'))
            [C,h2]=contour(RCM_vert_grid.lonlat, RCM_vert_grid.cut_depth2, squeeze(mean_data), [5, 10, 15, 20, 25, 30], 'color','k', ...
                    'linewidth', 1.5, 'linestyle', '-');
                clabel(C,h2,'FontSize',13,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        end
        set(gca, 'box', 'on', 'fontsize', param.m_grid_fontsize);
%         m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
        if min(RCM_info.years) == max(RCM_info.years)
            tmp.titlename = strcat(tmp.variable(6:end), ', ', RCM_info.season(1:4), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
        else
            tmp.titlename = strcat(tmp.variable(6:end), ', ', RCM_info.season(1:4), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
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
            case {'vert_temp', 'vert_salt', 'vert_rho'}
                caxis(param.colorbar_lev);
            case {'vert_u', 'vert_v'}
                caxis([-max(abs(mean_data(:))), max(abs(mean_data(:)))]);
                colormap(cmaps.byrmap);
            case{'vert_temp_nogrd', 'vert_salt_nogrd', 'vert_rho_nogrd'}
                colormap(jet)
%                 caxis([min(mean_data(:)) max(mean_data(:))])
%                 max(abs(mean_data(:)))-mean(abs(mean_data(:)), 'omitnan')
        end

        disp(['M = ', num2str(tmp.m_value)]);
        param.text_pos_x = (RCM_info.vert_section(2)-RCM_info.vert_section(1))/20+RCM_info.vert_section(1);
        param.text_pos_y = (RCM_info.vert_section(6)-RCM_info.vert_section(5))/20+RCM_info.vert_section(5);

        text(param.text_pos_x, param.text_pos_y, ['M = ', num2str(round(tmp.m_value,2))], 'FontSize', param.m_quiver_ref_text_fontsize); 

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        clear mean_data
%         RCM_vert_grid=rmfield(RCM_vert_grid, 'lon_rho');
        clear RCM_vert_grid;

    end  %% fig flag
end  %% var

function wstrcurl=get_curl(uwstr, vwstr, xdist, ydist)
    loncount=size(uwstr,1);
    latcount=size(uwstr,2);
    if (size(uwstr,1) == size(vwstr,1))
        wstrcurl=NaN(loncount, latcount);
        for i=1:loncount-2
            for j=1:latcount-2
                wstrcurl(i+1,j+1) = (vwstr(i+2,j+1)-vwstr(i,j+1)) ...
                    / (xdist(i+1, j+1)*2) + ...
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
                    / (xdist(i, j)*2) + ...
                    (uwstr(i,j+1)-uwstr(i,j)) ...
                    / (ydist(i, j)*2);  % N m-2 m-1
            end
        end
    end
end