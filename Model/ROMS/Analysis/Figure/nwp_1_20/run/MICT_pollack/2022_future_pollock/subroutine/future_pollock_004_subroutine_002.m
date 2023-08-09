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
    tmp.tifname=strcat(dirs.figdir, tmp.testname, '_surf_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
    if (exist(tmp.tifname , 'file') ~= 2 || fig_flag==2)      
        run(tmp.param_script);
%% set data file name (mat) for writing
        tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
        '_', RCM_info.season, '.mat'];
        switch tmp.testname
            case 'ENS4_hist'
               tmp.tset= {'test2117', 'test2118','test2119','test2120'};
               for titi=1:4
                    tmp.ttname=tmp.tset{titi};
                    dirs.ttmatdir = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/', tmp.ttname, '/run/mean/');
                    tmp.matname = [dirs.ttmatdir, tmp.ttname, '_', tmp.regionname, '_', tmp.variable,...
                        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                    if (exist(tmp.matname , 'file') ~= 2)
                        disp('please get data from subroutine_004_001 first')
                    else
                        load(tmp.matname, 'RCM_data', 'RCM_grid');
                    end 
                    tmp.tmean(:,:,titi)=RCM_data.mean;
               end
               RCM_data.mean=mean(tmp.tmean,3);
           case 'ENS4_fut'
               tmp.tset= {'test2127', 'test2128','test2129','test2130'};
               for titi=1:4
                    tmp.ttname=tmp.tset{titi};
                    dirs.ttmatdir = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/', tmp.ttname, '/run/mean/');
                    tmp.matname = [dirs.ttmatdir, tmp.ttname, '_', tmp.regionname, '_', tmp.variable,...
                        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                    if (exist(tmp.matname , 'file') ~= 2)
                        disp('please get data from subroutine_004_001 first')
                    else
                        load(tmp.matname, 'RCM_data', 'RCM_grid');
                    end 
                    tmp.tmean(:,:,titi)=RCM_data.mean;
               end
               RCM_data.mean=mean(tmp.tmean,3);
            otherwise
                if (exist(tmp.matname , 'file') ~= 2)
                    disp('please get data from subroutine_004_001 first')
                else
                    load(tmp.matname, 'RCM_data', 'RCM_grid');
                end 
        end
        
        if(strcmp(tmp.variable, 'SST'))
            RCM_data.mean=RCM_data.mean;
        elseif (strcmp(tmp.variable,'SSH')==1)
            ssh_correction_for_fig = Func_0017_SSH_correction_for_CMIP6_RMSE(tmp.testname);
            RCM_data.mean=RCM_data.mean-ssh_correction_for_fig;
        end

        RCM_grid.mask_model = double(inpolygon(RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
        RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
        RCM_data.mean=RCM_data.mean.*RCM_grid.mask_model;
        
        %% READ OBS, RMSE, CORR
        OBS_info.filename=['/Volumes/kyy_raid/Data/Observation/OISST/monthly_kimyy/', 'avhrr_only_monthly_v2_1995-2014_JanFeb.nc'];
        OBS_grid.lon=ncread(OBS_info.filename, 'lon');
        OBS_grid.lat=ncread(OBS_info.filename, 'lat');
        OBS_data.raw=ncread(OBS_info.filename, 'temp');
         [OBS_grid.lat_rho, OBS_grid.lon_rho] = meshgrid(OBS_grid.lat, OBS_grid.lon);

        [OBS_grid.refpolygon, OBS_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

        [OBS_grid.lon_min, OBS_grid.lon_max, OBS_grid.lat_min, OBS_grid.lat_max] = ...
                                    findind_Y(1/20, OBS_grid.domain(1:4), OBS_grid.lon_rho, OBS_grid.lat_rho);
        OBS_grid.cut_lon_rho = ...
            OBS_grid.lon_rho(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1));
        OBS_grid.cut_lat_rho = ...
            OBS_grid.lat_rho(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1));
        OBS_data.cut_raw=OBS_data.raw(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1),:);
        
        OBS_grid.mask_model = double(inpolygon(OBS_grid.cut_lon_rho,OBS_grid.cut_lat_rho,OBS_grid.refpolygon(:,1),OBS_grid.refpolygon(:,2)));
        OBS_grid.mask_model(OBS_grid.mask_model==0)=NaN;
        OBS_data.mean=mean(OBS_data.cut_raw, 3).*OBS_grid.mask_model;

        for i=1:size(OBS_data.cut_raw,3)
            RCM_data.interp_all(:,:,i)=griddata(RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho, RCM_data.all(:,:,ceil(i/2),2-mod(i,2)), ...
                OBS_grid.cut_lon_rho, OBS_grid.cut_lat_rho);
            RCM_data.sqe(:,:,i)=(RCM_data.interp_all(:,:,i)-OBS_data.cut_raw(:,:,i)).^2;
        end
        RCM_data.rmse=sqrt(mean(RCM_data.sqe,3));
        RCM_data.interp_mean=mean(RCM_data.interp_all,3);
        
        %% mean of RMSE
        [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data.rmse, OBS_grid.cut_lon_rho, OBS_grid.cut_lat_rho);

        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[127, 130, 37, 41],OBS_grid.cut_lon_rho,OBS_grid.cut_lat_rho); % southern EKB
        EKB_lon=RCM_grid.lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        EKB_lat=RCM_grid.lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        EKB_data=RCM_data.rmse(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [EKB_mean, error_status] = Func_0011_get_area_weighted_mean(EKB_data, EKB_lon, EKB_lat);

        %% pattern corr
        tmp.mask=RCM_data.interp_mean.*OBS_data.mean;
        tmp.a=RCM_data.interp_mean(isfinite(tmp.mask));
        tmp.b=OBS_data.mean(isfinite(tmp.mask));
        RCM_data.pattern_corr=corrcoef(tmp.a(:), tmp.b(:));

        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
        hold on;

        [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data.mean, RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
        m_pcolor(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho',RCM_data.mean');
        shading(gca,param.m_pcolor_shading_method);   
        if(strcmp(tmp.variable, 'SST'))
            [C,h2]=m_contour(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho', RCM_data.mean', [2, 5, 10], 'color','k', ...
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
            case {'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'u', 'v', 'swrad', 'shflux'}
                caxis(param.colorbar_lev);
            case {'wstrcurl', 'wcurl' }
%                 caxis([-max(abs(RCM_data.mean(:))), max(abs(RCM_data.mean(:)))]);
                caxis(param.colorbar_lev);
                colormap(cmaps.byrmap);
        end

        disp(['M = ', num2str(tmp.m_value)]);
%         m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        clear RCM_data.mean
        RCM_grid=rmfield(RCM_grid, 'lon_rho');
    end
end