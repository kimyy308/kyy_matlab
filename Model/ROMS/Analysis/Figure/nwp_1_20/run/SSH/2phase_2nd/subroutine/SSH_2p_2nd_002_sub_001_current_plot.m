% %  Updated 09-Oct-2021 by Yong-Yub Kim, structure


% start-------------------- earlier decadal current plot

dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'vec', tmp.fs, ...
    num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
    mkdir(strcat(dirs.figdir));
end 
%             outfile = strcat(figdir,regionname,);

tmp.tifname=strcat(dirs.figdir, tmp.testname, '_clim_uv_',num2str(min(RCM_info.years),'%04i'), ...
    '_',num2str(max(RCM_info.years),'%04i'), ...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.tif'); %% ~_year_month.jpg
tmp.variable='UV';
if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(1)==2)      
    tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable, ...
        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'), ...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.mat'];
    if (exist(tmp.matname , 'file') ~= 2 || flags.fig_switch(1)==2)
        for yearij=1:length(RCM_info.years)
            tmp.tempyear=RCM_info.years(yearij);
            tmp.yearstr=num2str(tmp.tempyear, '%04i');
            for monthij=1:length(RCM_info.months)
                tmp.tempmonth=RCM_info.months(monthij);
                tmp.monthstr=num2str(tmp.tempmonth, '%02i');
%                         D:\Data\Model\ROMS\nwp_1_20\backup_surf\test2102\run\v\1985\pck_test2102_v_monthly_1985_01.nc
                switch tmp.testname
                    case {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
                        tmp.ufilename=[dirs.filedir, 'u', tmp.fs, tmp.yearstr, filesep,'pck_', tmp.testname, '_u_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                        tmp.vfilename=[dirs.filedir, 'v', tmp.fs, tmp.yearstr, filesep,'pck_', tmp.testname, '_v_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                    otherwise
                        tmp.ufilename=[dirs.filedir, 'u', tmp.fs, tmp.yearstr, filesep,'NWP_pck_', tmp.testname, '_u_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                        tmp.vfilename=[dirs.filedir, 'v', tmp.fs, tmp.yearstr, filesep,'NWP_pck_', tmp.testname, '_v_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                end
                
               
                if (isfield(RCM_grid, 'lon_rho') ~= 1)
                    for gridi=1:length(RCM_grid.gridname)
                        RCM_grid.(RCM_grid.gridname{gridi})=ncread(RCM_grid.(['filename_', RCM_grid.gridname{gridi}]), RCM_grid.gridname{gridi});
                    end
%                             lat_rho=ncread(filename, 'lat_rho');
%                             lon_u=ncread(filename, 'lon_u');
%                             lat_u=ncread(filename, 'lat_u');
%                             lon_v=ncread(filename, 'lon_v');
%                             lat_v=ncread(filename, 'lat_v');

                    [RCM_grid.lon_min, RCM_grid.lon_max, RCM_grid.lat_min, RCM_grid.lat_max] = ...
                        findind_Y(1/20, RCM_grid.domain(1:4), RCM_grid.lon_rho, RCM_grid.lat_rho);
%                             [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(1/20, lonlat(1:4), lon_u, lat_u);
%                             [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(1/20, lonlat(1:4), lon_v, lat_v);
                    cut_lon_rho = ...
                        RCM_grid.lon_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                    cut_lat_rho = ...
                        RCM_grid.lat_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                end
                tmp.data_info = ncinfo(tmp.ufilename, 'u'); 
                tmp.u = ncread(tmp.ufilename,'u', ...
                    [RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                    [RCM_grid.lon_max(1)-RCM_grid.lon_min(1) RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                tmp.v = ncread(tmp.vfilename,'v', ...
                    [RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                    [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1) 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                if (isfield(tmp, 'mean_u') ~= 1)
                    tmp.mean_u=zeros(size(tmp.u));
                    tmp.mean_v=zeros(size(tmp.v));
                end
                tmp.mean_u=tmp.mean_u + (tmp.u / (length(RCM_info.years) * length(RCM_info.months)));
                tmp.mean_v=tmp.mean_v + (tmp.v / (length(RCM_info.years) * length(RCM_info.months)));
            end
        end
        u_rho = u2rho_2d(tmp.mean_u')';
        v_rho = v2rho_2d(tmp.mean_v')';
        if (isfield(tmp, 'ref_vec_x_range') ~= 1)
            tmp.ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location) ...
                == min(abs(cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location)));
            tmp.ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location) ...
                == min(abs(cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location)))+param.m_quiver_y_interval*2;
            switch tmp.regionname
                case 'AKP4'
                    tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2)) : ...
                        round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2))+param.m_quiver_x_interval;
                    tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2)) : ...
                        round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2))+param.m_quiver_y_interval;
                otherwise
                    tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2)) : ...
                        round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2))+param.m_quiver_x_interval*3;
                    tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2)) : ...
                        round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2))+param.m_quiver_y_interval*3;
            end
        end
        u_rho(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_u_value;
        v_rho(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_v_value;     
        save(tmp.matname, 'u_rho','v_rho', 'cut_lon_rho', 'cut_lat_rho');
    else
        load(tmp.matname);
    end

    RCM_grid.mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
    RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
    u_rho=u_rho(1:size(RCM_grid.mask_model,1),:).*RCM_grid.mask_model;
    v_rho=v_rho(:,1:size(RCM_grid.mask_model,2)).*RCM_grid.mask_model;  

    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
    hold on;

    m_gshhs_i('color',param.m_gshhs_line_color)  
    m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

    uvplot=m_quiver(cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                    cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                    u_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                    v_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                    'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);

    m_text(param.m_quiver_ref_text_x_location, param.m_quiver_ref_text_y_location, param.m_quiver_ref_text, 'FontSize', param.m_quiver_ref_text_fontsize); 
    m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
    if min(RCM_info.years) == max(RCM_info.years)
        param.titlename = strcat('UV, ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),') ');  %% + glacier contribution
    else
        param.titlename = strcat('UV, ', RCM_info.season(1:3), ', ', tmp.abb, ',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');  %% + glacier contribution
    end
    title(param.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
    set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
    saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
%                 system(['magick ', tifname, ' -trim ', tifname]);

    close all;
%             clear RCM_grid.lon_rho tmp.mean_u tmp.ref_vec_x_range
    if (isfield(RCM_grid, 'lon_rho') == 1)
        RCM_grid=rmfield(RCM_grid, 'lon_rho');
    end
    if (isfield(tmp, 'mean_u') == 1)
        tmp=rmfield(tmp, 'mean_u');
    end
    if (isfield(tmp, 'ref_vec_x_range') == 1)
        tmp=rmfield(tmp, 'ref_vec_x_range');
    end
end

% % % % start-------------------- later decadal current plot
% % %         fig_flag=fig_flags{2,2};
% % %         while (fig_flag)
% % %             tifname=strcat(outfile, '_', testname,'_', regionname, '_clim_uv_',num2str(min(inputyear2),'%04i'), ...
% % %                 '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
% % %             variable='UV';        
% % %             if (exist(tifname , 'file') ~= 2 || fig_flag==2)      
% % %                 matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
% % %                 if (exist(matname , 'file') ~= 2 || fig_flag==2)
% % %                     for yearij=1:length(inputyear2)
% % %                         tmp.tempyear=inputyear2(yearij);
% % %                         tmp.yearstr=num2str(tmp.tempyear, '%04i');
% % %                         for monthij=1:length(RCM_info.months)
% % %                             tmp.tempmonth=RCM_info.months(monthij);
% % %                             tmp.monthstr=num2str(tmp.tempmonth, '%02i');
% % %                             filename=[filedir, tmp.yearstr, filesep,'pck_', testname, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
% % % 
% % %                             if (exist('lon_rho' , 'var') ~= 1)
% % %                                 lon_rho=ncread(filename, 'lon_rho');
% % %                                 lat_rho=ncread(filename, 'lat_rho');
% % %                                 lon_u=ncread(filename, 'lon_u');
% % %                                 lat_u=ncread(filename, 'lat_u');
% % %                                 lon_v=ncread(filename, 'lon_v');
% % %                                 lat_v=ncread(filename, 'lat_v');
% % %                                 [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
% % %                                 [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(1/20, lonlat(1:4), lon_u, lat_u);
% % %                                 [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(1/20, lonlat(1:4), lon_v, lat_v);
% % %                                 cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % %                                 cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % %                             end
% % %                             data_info = ncinfo(filename, 'u'); 
% % %                             u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) data_info.Size(3) 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
% % %                             v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) data_info.Size(3) 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
% % % 
% % %                             if (exist('mean_u' , 'var') ~= 1)
% % %                                 mean_u=zeros(size(u));
% % %                                 mean_v=zeros(size(v));
% % %                             end
% % %                             mean_u=mean_u + (u / (length(inputyear2) * length(RCM_info.months)));
% % %                             mean_v=mean_v + (v / (length(inputyear2) * length(RCM_info.months)));
% % %                         end
% % %                     end
% % %                     u_rho = u2rho_2d(mean_u')';
% % %                     v_rho = v2rho_2d(mean_v')';
% % %                     if (exist('tmp.ref_vec_x_range' , 'var') ~= 1)
% % %                         tmp.ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
% % %                         tmp.ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
% % %                         tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(m_quiver_x_interval/2)) : round(tmp.ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
% % %                         tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(m_quiver_y_interval/2)) : round(tmp.ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
% % %                     end
% % %                     u_rho(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=m_quiver_ref_u_value;
% % %                     v_rho(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=m_quiver_ref_v_value;     
% % %                     save(matname, 'u_rho','v_rho', 'cut_lon_rho', 'cut_lat_rho');
% % %                 else
% % %                     load(matname);
% % %                 end
% % %                 
% % %                 mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
% % %                 mask_model(mask_model==0)=NaN;
% % %                 u_rho=u_rho(1:size(mask_model,1),:).*mask_model;
% % %                 v_rho=v_rho(:,1:size(mask_model,2)).*mask_model;  
% % %                 
% % %                 m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % %                 hold on;
% % %                 m_gshhs_i('color',m_gshhs_line_color)  
% % %                 m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % 
% % %                 uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
% % %                                 cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
% % %                                 u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
% % %                                 v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
% % %                                 'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
% % % 
% % %                 m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
% % %                 m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % %                 titlename = strcat('UV mean, ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
% % %                 title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % 
% % %                 set(gcf, 'PaperUnits', 'points');
% % %                 set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % %                 set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % %                 saveas(gcf,tifname,'tif');
% % %                 RemoveWhiteSpace([], 'file', tifname);
% % %                 close all;
% % %                 clear lon_rho mean_u tmp.ref_vec_x_range
% % %             end
% % %             fig_flag=0;
% % %         end