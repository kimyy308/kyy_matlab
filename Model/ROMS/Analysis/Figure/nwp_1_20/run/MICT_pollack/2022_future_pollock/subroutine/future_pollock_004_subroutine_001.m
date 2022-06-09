% %  Updated 26-Apr-2021 by Yong-Yub Kim, 

% start-------------------- earlier decadal SST, SSS plot
for varind2=1:length(RCM_info.vars)
%% set variable name & figure directory
    tmp.variable=RCM_info.vars{varind2};
    dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 

%% set figure file name
    tmp.tifname=strcat(dirs.figdir, tmp.testname, '_ts_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
    if (exist(tmp.tifname , 'file') ~= 2 || fig_flag==2)      
        run(tmp.param_script);
%% set data file name (mat) for writing
        tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
        '_', RCM_info.season, '.mat'];
        if (exist(tmp.matname , 'file') ~= 2 || fig_flag==2)
            for yearij=1:length(RCM_info.years)  %% yearly loop
                tmp.tempyear=RCM_info.years(yearij);
                tmp.yearstr=num2str(tmp.tempyear, '%04i');
                for monthij=1:length(RCM_info.months)  %% monthly loop
                    tmp.tempmonth=RCM_info.months(monthij);
                    tmp.monthstr=num2str(tmp.tempmonth, '%02i');
%% set nc file name for reading
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
                    
%% get grid information                    
                    if (isfield(RCM_grid, 'lon_rho') ~= 1)
                        for gridi=1:length(RCM_grid.gridname)
                            RCM_grid.(RCM_grid.gridname{gridi})=ncread(RCM_grid.(['filename_', RCM_grid.gridname{gridi}]), RCM_grid.gridname{gridi});
                        end

                        [RCM_grid.lon_min, RCM_grid.lon_max, RCM_grid.lat_min, RCM_grid.lat_max] = ...
                            findind_Y(1/20, RCM_grid.domain(1:4), RCM_grid.lon_rho, RCM_grid.lat_rho);
                        RCM_grid.cut_lon_rho = ...
                            RCM_grid.lon_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                        RCM_grid.cut_lat_rho = ...
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

%% get data
                    if (strcmp(tmp.variable,'SST')==1 || strcmp(tmp.variable,'SSS')==1)
                        tmp.data = ncread(tmp.filename,param.varname,[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    elseif (strcmp(tmp.variable,'SSH')==1 || strcmp(tmp.variable,'Uwind')==1 || strcmp(tmp.variable,'Vwind')==1 ...
                            || strcmp(tmp.variable,'shflux')==1 || strcmp(tmp.variable,'swrad')==1)
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
                        tmp.data=griddata(cut_lon_psi,cut_lat_psi,tmp.wstrcurl,RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho);
                    end
                    tmp.data(tmp.data > 10^6)=NaN;
                    
%                     if (isfield(RCM_data, 'mean') ~= 1)
%                         RCM_data.mean=zeros(size(tmp.data));
%                     end
                    if yearij==1 & monthij==1
                        RCM_data.mean=zeros(size(tmp.data));
                    end

                    RCM_data.mean=RCM_data.mean + ...
                        (tmp.data / (length(RCM_info.years) * length(RCM_info.months)));
                    
%                     RCM_data.all(:,:,(yearij-1)*12+monthij)=tmp.data;
                    RCM_data.all(:,:,yearij,monthij)=tmp.data;
                end %% monthly loop
            end  %% yearly loop
            tmp.data_all=reshape(RCM_data.all, [size(tmp.data,1), size(tmp.data,2), length(RCM_info.years)*length(RCM_info.months)]);
            RCM_data.yearly_all=mean(RCM_data.all,4);
            for i=1:size(tmp.data,1)
                for j =1:size(tmp.data,2)
                    RCM_data.std(i,j)=std(tmp.data_all(i,j,:));
                    RCM_data.std_yearly(i,j)=std(RCM_data.yearly_all(i,j,:));
                end
            end
            [RCM_data.m_value, error_status] = Func_0011_get_area_weighted_mean(mean(tmp.data_all,3), RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
            [RCM_data.yearly_m_value, error_status] = Func_0011_get_area_weighted_mean(mean(RCM_data.yearly_all,3), RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
            [RCM_data.std_m_value, error_status] = Func_0011_get_area_weighted_mean(RCM_data.std, RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
            [RCM_data.std_yearly_m_value, error_status] = Func_0011_get_area_weighted_mean(RCM_data.std_yearly, RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
            for t=1:size(tmp.data_all,3)
                [RCM_data.spa_mean(t), error_status]= Func_0011_get_area_weighted_mean(tmp.data_all(:,:,t), RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
            end
            for t=1:size(RCM_data.yearly_all,3)
                [RCM_data.spa_yearly_mean(t), error_status]= Func_0011_get_area_weighted_mean(RCM_data.yearly_all(:,:,t), RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
            end
%             a=1
            save(tmp.matname, 'RCM_data', 'RCM_grid', 'RCM_info');
            disp('future_pollock_004_subroutine_001');
        else
            load(tmp.matname, 'RCM_data', 'RCM_grid');
        end 
% % % 
% % %         
% % %         
% % % %         plot(RCM_info.years,spa_RCM_data.mean); set(gca, 'fontsize', 20)
% % % %         pcolor(RCM_grid.cut_lon_rho', RCM_grid.cut_lat_rho', comb_mean_yearly_data(:,:,end)'); shading flat; colorbar; set(gca, 'fontsize', 20)
% % %         
% % %         if(strcmp(tmp.variable, 'SST'))
% % %             RCM_data.mean=RCM_data.mean;
% % %         elseif (strcmp(tmp.variable,'SSH')==1)
% % % %                         ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor(tmp.testname);
% % %             ssh_correction_for_fig = Func_0017_SSH_correction_for_CMIP6_RMSE(tmp.testname);
% % %             RCM_data.mean=RCM_data.mean-ssh_correction_for_fig;
% % %         end
% % % 
% % %         RCM_grid.mask_model = double(inpolygon(RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
% % %         RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
% % %         RCM_data.mean=RCM_data.mean.*RCM_grid.mask_model;
% % % 
% % %         m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
% % %         hold on;
% % % 
% % %         [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data.mean, RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
% % %         m_pcolor(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho',RCM_data.mean');
% % %         shading(gca,param.m_pcolor_shading_method);   
% % %         if(strcmp(tmp.variable, 'SST'))
% % % %                         [C,h2]=m_contour(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho', RCM_data.mean', [2, 4, 6, 8, 10], 'color','k', ...
% % % %                                 'linewidth', 1.5, 'linestyle', '-');
% % %             [C,h2]=m_contour(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho', RCM_data.mean', [5, 10, 15, 20, 25, 30], 'color','k', ...
% % %                     'linewidth', 1.5, 'linestyle', '-');
% % %                 clabel(C,h2,'FontSize',13,'Color','k', ...
% % %                     'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
% % %         end
% % % 
% % %         
% % %         m_gshhs_i('color',param.m_gshhs_line_color)  
% % %         m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land
% % % 
% % %         m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
% % %         if min(RCM_info.years) == max(RCM_info.years)
% % %             tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
% % %         else
% % %             tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
% % %         end
% % %         title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title
% % % 
% % %         % set colorbar 
% % %         h = colorbar;
% % % %                     if (strcmp(tmp.variable,'SSH')==1)
% % % %                         colormap(flip(cool));
% % % %                     else
% % %             colormap(jet);
% % % %                     end
% % % 
% % %         set(h,'fontsize',param.colorbar_fontsize);
% % %         title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
% % %         switch tmp.variable
% % %             case {'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'u', 'v', }
% % %                 caxis(param.colorbar_lev);
% % %             case {'wstrcurl', 'wcurl'}
% % % %                 caxis([-max(abs(RCM_data.mean(:))), max(abs(RCM_data.mean(:)))]);
% % %                 caxis(param.colorbar_lev);
% % %                 colormap(cmaps.byrmap);
% % %         end
% % % 
% % %         disp(['M = ', num2str(tmp.m_value)]);
% % % %         m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 
% % % 
% % %         set(gcf, 'PaperUnits', 'points');
% % %         set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
% % %         set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
% % %         saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
% % %         close all;
        clear RCM_data.mean
        RCM_grid=rmfield(RCM_grid, 'lon_rho');
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