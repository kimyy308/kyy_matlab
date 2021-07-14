close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region ={'ECS(CRD)'}

for regionind=1:length(all_region)
    clearvars '*' -except regionind all_region

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
%         addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
        addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end


    shadlev = [30 35];
    windlev = [-10 10];
    rms_shadlev = [0 4];
    bias_shadlev = [-4 4];
    conlev  = shadlev(1) : 5 : shadlev(2);
    dl=1/10;
    % for snu_desktop
    testname='test08'   % % need to change
    inputyear = [2009:2018]; % % put year which you want to plot [year year ...]
    inputmonth = [8]; % % put month which you want to plot [month month ...]
    

    
    varname ='salt';  %% reference variable -> temperature
    run('nwp_polygon_point.m');
    regionname=all_region{regionind}
    switch(regionname)
        case('NWP') %% North western Pacific
            lonlat = [115, 164, 15, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('ES') %% East Sea
            refpolygon=espolygon;
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        case('ECS(CRD)') %% for East China Sea Changjiang river discharge path
            refpolygon=ecs_crdpolygon;
        case('EKB') %% North western Pacific
            lonlat = [127, 131, 37, 42];  %% East Korea Bay
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
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

    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\EAST\3rd_year\figure\nwp_1_10\',testname,'\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_ECS_CWD.m';
        filedir = strcat('G:\nwp_1_10\output\test08\run\'); % % where data files are
        avhrrdir='E:\Data\Observation\OISST\monthly_kimyy\';
    elseif (strcmp(system_name,'GLNXA64'))
    end

    run(param_script);
    ind=1; % initialize time index of comb_data
    ind2=1;
    for yearij = 1:length(inputyear)
        tempyear = inputyear(yearij); % current year in yealy loop
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij); % current month in monthly loop
            tic; % getting runtime, start
            
            if inputmonth(monthij)==1
                daystart = 1;
                dayend = 31;
            else
                daystart = sum(eomday(tempyear,1:tempmonth-1))+1;
                dayend = daystart + eomday(tempyear,tempmonth)-1;
                disp([num2str(yearij), 'y_',num2str(monthij),'m, from avg file ', num2str(daystart,'%04i'), ' to ', num2str(dayend, '%04i') ])
            end
            
            % read model grid
            grdname = strcat(filedir, num2str(inputyear(yearij),'%04i'), '\', 'ocean_avg_', num2str(daystart,'%04i'), '.nc'); % get name of gridfile
            if (exist('lon')==0)
                modelinfo=ncinfo(grdname);
                for i=1:length(modelinfo.Dimensions)
                    if(strcmp(modelinfo.Dimensions(i).Name,'xi_rho')==1)
                        length_lon=modelinfo.Dimensions(i).Length;  % get length of longitude
                    elseif(strcmp(modelinfo.Dimensions(i).Name,'eta_rho')==1)
                        length_lat=modelinfo.Dimensions(i).Length;  % get length of latitude
                    elseif(strcmp(modelinfo.Dimensions(i).Name,'s_rho')==1)
                        length_depth=modelinfo.Dimensions(i).Length;  % get length of depth
                    end
                end
                lon_whole = ncread(grdname,'lon_rho',[1 1],[length_lon, length_lat]);
                lat_whole = ncread(grdname,'lat_rho',[1 1],[length_lon, length_lat]);
                
                [lon_min, lon_max, lat_min, lat_max] = findind_Y(dl, lonlat, lon_whole', lat_whole');  %% lon_whole' and lat_whole' should have [y x] dimension

                lon = ncread(grdname,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                lat = ncread(grdname,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_model(1:length_lon,1:length_lat)=1;
                    otherwise
                        mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                end
            end
            
%             for dayij = daystart : dayend 
            for dayij = daystart-7 : dayend  %% for aquiring 7 - day ago winds to lag-correlation 
                tempday = dayij - daystart + 1; % current date in daily loop
                filename = strcat(filedir, num2str(inputyear(yearij),'%04i'), '\', 'ocean_avg_', num2str(dayij,'%04i'), '.nc');
                
                lon_length_sliced = lon_max(1)-lon_min(1)+1;
                lat_length_sliced = lat_max(1)-lat_min(1)+1;
                data = ncread(filename,varname,[lon_min(1), lat_min(1), length_depth, 1], [lon_length_sliced, lat_length_sliced, 1, 1]);  % get surface data
                uwind = ncread(filename,'Uwind',[lon_min(1), lat_min(1), 1], [lon_length_sliced, lat_length_sliced, 1]);  % get surface data
                vwind = ncread(filename,'Vwind',[lon_min(1), lat_min(1), 1], [lon_length_sliced, lat_length_sliced, 1]);  % get surface data
                u = ncread(filename,'u',[lon_min(1)+1, lat_min(1), length_depth, 1], [lon_length_sliced-1, lat_length_sliced, 1, 1]);  % get surface data
                v = ncread(filename,'v',[lon_min(1), lat_min(1)+1, length_depth, 1], [lon_length_sliced, lat_length_sliced-1, 1, 1]);  % get surface data
                
%                 same as u2rho_2d in romstools
                u_rho(2:lon_length_sliced-1,:) = 0.5 * (u(1:lon_length_sliced-2,:) + u(2:lon_length_sliced-1,:));
                u_rho(1,:)=u_rho(2,:);
                u_rho(lon_length_sliced,:)=u_rho(lon_length_sliced-1,:);
%                 same as v2rho_2d in romstools
                v_rho(:,2:lat_length_sliced-1) = 0.5 * (v(:,1:lat_length_sliced-2) + v(:,2:lat_length_sliced-1));
                v_rho(:,1) = v_rho(:,2);
                v_rho(:,lat_length_sliced) = v_rho(:,lat_length_sliced-1);
                
                if (exist('ref_vec_x_range' , 'var') ~= 1)
                    ref_vec_x_ind = find(abs(lon(:,1)-m_quiver_ref_text_x_location) == min(abs(lon(:,1)-m_quiver_ref_text_x_location)));
                    ref_vec_y_ind = find(abs(lat(1,:)-m_quiver_ref_text_y_location) == min(abs(lat(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                    ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval-1;
                    ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval-1;
                end
                u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;    
                uwind(ref_vec_x_range,ref_vec_y_range)=5.0;
                vwind(ref_vec_x_range,ref_vec_y_range)=0.00001;    
                
                if (dayij >=daystart)
                    comb_data(:,:,ind) = data;
                    comb_uwind(:,:,ind) = uwind;
                    comb_vwind(:,:,ind) = vwind;
                    comb_u(:,:,ind) = u_rho;
                    comb_v(:,:,ind) = v_rho;
                    comb_year(ind) = tempyear;
                    comb_month(ind) = tempmonth;
                    comb_day(ind) = tempday;
                    ind = ind + 1;
                end
                
                comb_data2(:,:,ind2) = data;
                comb_uwind2(:,:,ind2) = uwind;
                comb_vwind2(:,:,ind2) = vwind;
                comb_u2(:,:,ind2) = u_rho;
                comb_v2(:,:,ind2) = v_rho;
                comb_year2(ind2) = tempyear;
                comb_month2(ind2) = tempmonth;
                comb_day2(ind2) = tempday;
                ind2 = ind2 + 1;
                
            end
            toc;
        end
    end


% % test figure for salinity
%     figure(1)
%     pcolor(lon',lat',comb_data(:,:,310)')
%     shading flat
%     colorbar
%     caxis(shadlev)
%     title('salt')
%     
%     figure(4)
%     pcolor(lon',lat',comb_vwind(:,:,310)')
%     shading flat
%     colorbar
% %     caxis(shadlev)
%     title('vwind')
    
for i=1:lon_length_sliced
    for j=1:lat_length_sliced
        [R, P] = corrcoef(comb_data(i,j,:),comb_vwind(i,j,:));
        corr_data_southerly(i,j)=R(1,2);
        pval_data_southerly(i,j)=P(1,2);
        [R, P] = corrcoef(comb_data(i,j,:),-comb_uwind(i,j,:));
        corr_data_easterly(i,j)=R(1,2);
        pval_data_easterly(i,j)=P(1,2);
    end
end

for lag=1:7
    for i=1:lon_length_sliced
        for j=1:lat_length_sliced
            for year=2009:2018
                comb_vwind_lag(i,j,(year-2009)*tempday+1:(year-2008)*tempday) = ...
                    comb_vwind2(i,j,(year-2009)*(tempday+7)+8-lag:(year-2009)*(tempday+7)+8-lag +tempday-1);
                comb_uwind_lag(i,j,(year-2009)*tempday+1:(year-2008)*tempday) = ...
                    comb_uwind2(i,j,(year-2009)*(tempday+7)+8-lag:(year-2009)*(tempday+7)+8-lag +tempday-1);
            end
            [R, P] = corrcoef(comb_data(i,j,:),comb_vwind_lag(i,j,:));
            lag_corr_data_southerly(i,j,lag)=R(1,2);
            lag_pval_data_southerly(i,j,lag)=P(1,2);
            [R, P] = corrcoef(comb_data(i,j,:),-comb_uwind_lag(i,j,:));
            lag_corr_data_easterly(i,j,lag)=R(1,2);
            lag_pval_data_easterly(i,j,lag)=P(1,2);
        end
    end
end


% 
% % test figure for corrcoef
%     figure(2)   
%     pcolor(lon',lat',corr_data_southerly(:,:)')
%     shading flat
%     colorbar
%     caxis([-0.5 0.5])
%     colormap(bwrmap)
%     title('southerly(wind) & salt, corr coef')
%     
%     figure(3)
%     pcolor(lon',lat',corr_data_easterly(:,:)')
%     shading flat
%     colorbar
%     caxis([-0.5 0.5])
%     colormap(bwrmap)
%     title('easterly(wind) & salt, corr coef')
% 
%     figure(5)
%     pcolor(lon',lat',(pval_data_southerly(:,:)')*100.0)
%     shading flat
%     colorbar
%     caxis([0 20])
%     colormap(bwrmap)
%     title('easterly(wind) & salt, probablity value')
  
%   
for i= 1: size(comb_data,3)
%         figdir=[figrawdir,'salt\'];
%         outfile = strcat(figdir,regionname);
%         if (exist(strcat(figdir) , 'dir') ~= 7)
%             mkdir(strcat(figdir));
%         end 
% 
%         close all;
%         % avhrr climatological trend plot
%         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%         hold on;
% 
%         m_pcolor(lon',lat',comb_data(:,:,i)');
%         shading(gca,m_pcolor_shading_method);
%         m_gshhs_i('color',m_gshhs_line_color);
%         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
% %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
%         titlename = strcat(regionname,', ',num2str(comb_day(i),'%02i'),calendarname{comb_month(i)},num2str(comb_year(i),'%04i'),', Salt');
%         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%         % set colorbar 
%         h = colorbar;
%         colormap(colormap_style);
%         set(h,'fontsize',colorbar_fontsize);
% %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%         caxis(shadlev);
% 
%         % set grid
%         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%         set(gcf, 'PaperUnits', 'points');
%         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%         
%         tifname=strcat(outfile, '_salt',num2str(min(shadlev)),'_',num2str(comb_year(i),'%04i'),'_', ...
%             num2str(comb_month(i),'%02i'),'_',num2str(comb_day(i),'%02i'),'.tif'); %% ~_year_month.tif
%         drawnow;
%         saveas(gcf,tifname,'tif');
% 
%         disp(' ')
%         disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%         disp(' ')
%         disp([' File path is : ',tifname])
%         disp(' ')
% 
%         hold off

% % % % % % % % % % % % % % % %         uwind
%         figdir=[figrawdir,'uwind\'];
%         outfile = strcat(figdir,regionname);
%         if (exist(strcat(figdir) , 'dir') ~= 7)
%             mkdir(strcat(figdir));
%         end 
% 
%         close all;
%         % avhrr climatological trend plot
%         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%         hold on;
% 
%         m_pcolor(lon',lat',comb_uwind(:,:,i)');
%         shading(gca,m_pcolor_shading_method);
%         m_gshhs_i('color',m_gshhs_line_color);
%         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
% %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
%         titlename = strcat(regionname,', ',num2str(comb_day(i),'%02i'),calendarname{comb_month(i)},num2str(comb_year(i),'%04i'),', Uwind');
%         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%         % set colorbar 
%         h = colorbar;
%         colormap(bwrmap);
%         set(h,'fontsize',colorbar_fontsize);
%         title(h,'m/s','fontsize',colorbar_title_fontsize);
%         caxis(windlev);
% 
%         % set grid
%         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%         set(gcf, 'PaperUnits', 'points');
%         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%         
%         tifname=strcat(outfile,'_uwind_',num2str(comb_year(i),'%04i'),'_', ...
%             num2str(comb_month(i),'%02i'),'_',num2str(comb_day(i),'%02i'),'.tif'); %% ~_year_month.tif
%         drawnow;
%         saveas(gcf,tifname,'tif');
% 
%         disp(' ')
%         disp(['clim_', num2str(tempmonth), '_', 'uwind', ' plot is created.'])
%         disp(' ')
%         disp([' File path is : ',tifname])
%         disp(' ')
% 
%         hold off
%         
%         % % % % % % % % % % % % % % %         vwind
%         figdir=[figrawdir,'vwind\'];
%         outfile = strcat(figdir,regionname);
%         if (exist(strcat(figdir) , 'dir') ~= 7)
%             mkdir(strcat(figdir));
%         end 
% 
%         close all;
%         % avhrr climatological trend plot
%         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%         hold on;
% 
%         m_pcolor(lon',lat',comb_vwind(:,:,i)');
%         shading(gca,m_pcolor_shading_method);
%         m_gshhs_i('color',m_gshhs_line_color);
%         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
% %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
%         titlename = strcat(regionname,', ',num2str(comb_day(i),'%02i'),calendarname{comb_month(i)},num2str(comb_year(i),'%04i'),', Vwind');
%         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%         % set colorbar 
%         h = colorbar;
%         colormap(bwrmap);
%         set(h,'fontsize',colorbar_fontsize);
%         title(h,'m/s','fontsize',colorbar_title_fontsize);
%         caxis(windlev);
% 
%         % set grid
%         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%         set(gcf, 'PaperUnits', 'points');
%         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%         
%         tifname=strcat(outfile,'_vwind_',num2str(comb_year(i),'%04i'),'_', ...
%             num2str(comb_month(i),'%02i'),'_',num2str(comb_day(i),'%02i'),'.tif'); %% ~_year_month.tif
%         drawnow;
%         saveas(gcf,tifname,'tif');
% 
%         disp(' ')
%         disp(['clim_', num2str(tempmonth), '_', 'vwind', ' plot is created.'])
%         disp(' ')
%         disp([' File path is : ',tifname])
%         disp(' ')
% 
%         hold off

%         % % % % % % % % % % % % % % %         current
%         figdir=[figrawdir,'current\'];
%         outfile = strcat(figdir,regionname);
%         if (exist(strcat(figdir) , 'dir') ~= 7)
%             mkdir(strcat(figdir));
%         end 
%         
%         close all;
%         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%         hold on;
% 
%         m_pcolor(lon',lat',comb_data(:,:,i)');
%         shading(gca,m_pcolor_shading_method);
%         m_gshhs_i('color',m_gshhs_line_color);
%         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
%                         lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
%                         comb_u(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,i) * m_quiver_vector_size, ...
%                         comb_v(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,i) * m_quiver_vector_size, ...
%                         'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
% 
%         titlename = strcat(regionname,', ',num2str(comb_day(i),'%02i'),calendarname{comb_month(i)},num2str(comb_year(i),'%04i'),', current');
%         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%         
%         % set colorbar 
%         h = colorbar;
%         colormap(jet_mod);
%         set(h,'fontsize',colorbar_fontsize);
% %         title(h,'m/s','fontsize',colorbar_title_fontsize);
%         caxis(shadlev);
%         
%         % set grid
%         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
%         m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
% 
%         set(gcf, 'PaperUnits', 'points');
%         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%         
%         tifname=strcat(outfile,'_current_',num2str(comb_year(i),'%04i'),'_', ...
%             num2str(comb_month(i),'%02i'),'_',num2str(comb_day(i),'%02i'),'.tif'); %% ~_year_month.tif
%         drawnow;
%         saveas(gcf,tifname,'tif');
% 
%         disp(' ')
%         disp(['clim_', num2str(tempmonth), '_', 'current', ' plot is created.'])
%         disp(' ')
%         disp([' File path is : ',tifname])
%         disp(' ')
% 
%         hold off
% 
%         
%                 % % % % % % % % % % % % % % %        wind
%         figdir=[figrawdir,'wind\'];
%         outfile = strcat(figdir,regionname);
%         if (exist(strcat(figdir) , 'dir') ~= 7)
%             mkdir(strcat(figdir));
%         end 
%         
%         close all;
%         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%         hold on;
% 
%         m_pcolor(lon',lat',comb_data(:,:,i)');
%         shading(gca,m_pcolor_shading_method);
%         m_gshhs_i('color',m_gshhs_line_color);
%         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
%                         lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
%                         comb_uwind(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,i) * m_quiver_vector_size/20.0, ...
%                         comb_vwind(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,i) * m_quiver_vector_size/20.0, ...
%                         'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
% 
%         titlename = strcat(regionname,', ',num2str(comb_day(i),'%02i'),calendarname{comb_month(i)},num2str(comb_year(i),'%04i'),', wind');
%         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%         
%         % set colorbar 
%         h = colorbar;
%         colormap(jet_mod);
%         set(h,'fontsize',colorbar_fontsize);
% %         title(h,'m/s','fontsize',colorbar_title_fontsize);
%         caxis(shadlev);
%         
%         % set grid
%         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
%         m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, [num2str(5),' m/s'], 'FontSize', m_quiver_ref_text_fontsize); 
% 
%         set(gcf, 'PaperUnits', 'points');
%         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%         
%         tifname=strcat(outfile,'_wind_',num2str(comb_year(i),'%04i'),'_', ...
%             num2str(comb_month(i),'%02i'),'_',num2str(comb_day(i),'%02i'),'.tif'); %% ~_year_month.tif
%         drawnow;
%         saveas(gcf,tifname,'tif');
% 
%         disp(' ')
%         disp(['clim_', num2str(tempmonth), '_', 'wind', ' plot is created.'])
%         disp(' ')
%         disp([' File path is : ',tifname])
%         disp(' ')
% 
%         hold off
end


% % correlation, pval, significant correlation plot
for grouping=1:1
    % %                 % % % % % % % % % % % % % % %      southerly & salt
    % % % % % % % % % % % % % % % %                 correlation coeffcient
    % 
    % figdir=[figrawdir,'analysis\'];
    % outfile = strcat(figdir,regionname);
    % if (exist(strcat(figdir) , 'dir') ~= 7)
    %     mkdir(strcat(figdir));
    % end 
    % 
    % close all;
    % m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    % hold on;
    % 
    % m_pcolor(lon',lat',corr_data_southerly(:,:)');
    % shading(gca,m_pcolor_shading_method);
    % m_gshhs_i('color',m_gshhs_line_color);
    % m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    % 
    % %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
    % titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),', southerly(wind) & salt, corr coef');
    % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    % 
    % % set colorbar 
    % h = colorbar;
    % colormap(bwrmap);
    % set(h,'fontsize',colorbar_fontsize);
    % %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
    % caxis([-0.5 0.5]);
    % 
    % % set grid
    % m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    % 
    % set(gcf, 'PaperUnits', 'points');
    % set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    % set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    % 
    % tifname=strcat(outfile, '_southerly_salt_corr_coef','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    %     num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    % drawnow;
    % saveas(gcf,tifname,'tif');
    % 
    % disp(' ')
    % disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
    % disp(' ')
    % disp([' File path is : ',tifname])
    % disp(' ')
    % 
    % hold off
    % 
    % 
    %         
    % %                 % % % % % % % % % % % % % % %      easterly & salt
    % % % % % % % % % % % % % % % %                 correlation coeffcient
    % 
    % figdir=[figrawdir,'analysis\'];
    % outfile = strcat(figdir,regionname);
    % if (exist(strcat(figdir) , 'dir') ~= 7)
    %     mkdir(strcat(figdir));
    % end 
    % 
    % close all;
    % m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    % hold on;
    % 
    % m_pcolor(lon',lat',corr_data_easterly(:,:)');
    % shading(gca,m_pcolor_shading_method);
    % m_gshhs_i('color',m_gshhs_line_color);
    % m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    % 
    % %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
    % titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),', easterly(wind) & salt, corr coef');
    % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    % 
    % % set colorbar 
    % h = colorbar;
    % colormap(bwrmap);
    % set(h,'fontsize',colorbar_fontsize);
    % %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
    % caxis([-0.5 0.5]);
    % 
    % % set grid
    % m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    % 
    % set(gcf, 'PaperUnits', 'points');
    % set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    % set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    % 
    % tifname=strcat(outfile, '_easterly_salt_corr_coef','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    %     num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    % drawnow;
    % saveas(gcf,tifname,'tif');
    % 
    % disp(' ')
    % disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
    % disp(' ')
    % disp([' File path is : ',tifname])
    % disp(' ')
    % 
    % hold off

    % 
    % %                 % % % % % % % % % % % % % % %      southerly & salt
    % % % % % % % % % % % % % % % %                 probability value (p-val),
    % % % % % % % % % % % % % % %  significant level --> generally 0.05< or 0.01<
    % 
    % figdir=[figrawdir,'analysis\'];
    % outfile = strcat(figdir,regionname);
    % if (exist(strcat(figdir) , 'dir') ~= 7)
    %     mkdir(strcat(figdir));
    % end 
    % 
    % close all;
    % m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    % hold on;
    % 
    % m_pcolor(lon',lat',pval_data_southerly(:,:)');
    % shading(gca,m_pcolor_shading_method);
    % m_gshhs_i('color',m_gshhs_line_color);
    % m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    % 
    % %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
    % titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),', southerly(wind) & salt, corr-pval');
    % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    % 
    % % set colorbar 
    % h = colorbar;
    % colormap(bwrmap);
    % set(h,'fontsize',colorbar_fontsize);
    % %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
    % caxis([0 0.1]);
    % 
    % % set grid
    % m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    % 
    % set(gcf, 'PaperUnits', 'points');
    % set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    % set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    % 
    % tifname=strcat(outfile, '_southerly_salt_corr_pval','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    %     num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    % drawnow;
    % saveas(gcf,tifname,'tif');
    % 
    % disp(' ')
    % disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
    % disp(' ')
    % disp([' File path is : ',tifname])
    % disp(' ')
    % 
    % hold off
    % 
    %     
    % %                 % % % % % % % % % % % % % % %      easterly & salt
    % % % % % % % % % % % % % % % %                 probability value (p-val),
    % % % % % % % % % % % % % % %  significant level --> generally 0.05< or 0.01<
    % 
    % figdir=[figrawdir,'analysis\'];
    % outfile = strcat(figdir,regionname);
    % if (exist(strcat(figdir) , 'dir') ~= 7)
    %     mkdir(strcat(figdir));
    % end 
    % 
    % close all;
    % m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    % hold on;
    % 
    % m_pcolor(lon',lat',pval_data_easterly(:,:)');
    % shading(gca,m_pcolor_shading_method);
    % m_gshhs_i('color',m_gshhs_line_color);
    % m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    % 
    % %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
    % titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),', easterly(wind) & salt, corr-pval');
    % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    % 
    % % set colorbar 
    % h = colorbar;
    % colormap(bwrmap);
    % set(h,'fontsize',colorbar_fontsize);
    % %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
    % caxis([0 0.1]);
    % 
    % % set grid
    % m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    % 
    % set(gcf, 'PaperUnits', 'points');
    % set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    % set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    % 
    % tifname=strcat(outfile, '_easterly_salt_corr_pval','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    %     num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    % drawnow;
    % saveas(gcf,tifname,'tif');
    % 
    % disp(' ')
    % disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
    % disp(' ')
    % disp([' File path is : ',tifname])
    % disp(' ')
    % 
    % hold off


    % 
    % %                 % % % % % % % % % % % % % % %      southerly & salt
    % % % % % % % %                 correlation coefficient
    % % % % % % % % % % % % % % % %                 significance level is higher than 95
    % % % % % % % % % % % % % % %  significant level --> generally 0.05< or 0.01<
    % 
    % figdir=[figrawdir,'analysis\'];
    % outfile = strcat(figdir,regionname);
    % if (exist(strcat(figdir) , 'dir') ~= 7)
    %     mkdir(strcat(figdir));
    % end 
    % 
    % close all;
    % m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    % hold on;
    % pval_mask_southerly=pval_data_southerly(:,:);
    % pval_mask_southerly(pval_data_southerly(:,:)<=0.05)=1.0;
    % pval_mask_southerly(pval_data_southerly(:,:)>0.05)=0;
    %     
    % m_pcolor(lon',lat',corr_data_southerly(:,:)' .* pval_mask_southerly');
    % shading(gca,m_pcolor_shading_method);
    % m_gshhs_i('color',m_gshhs_line_color);
    % m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    % 
    % %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
    % titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),', southerly(wind) & salt, corr, sig-lev 95%');
    % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    % 
    % % set colorbar 
    % h = colorbar;
    % colormap(bwrmap);
    % set(h,'fontsize',colorbar_fontsize);
    % %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
    % caxis([-1 1]);
    % 
    % % set grid
    % m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    % 
    % set(gcf, 'PaperUnits', 'points');
    % set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    % set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    % 
    % tifname=strcat(outfile, '_southerly_salt_corr_sig_95_','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    %     num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    % drawnow;
    % saveas(gcf,tifname,'tif');
    % 
    % disp(' ')
    % disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
    % disp(' ')
    % disp([' File path is : ',tifname])
    % disp(' ')
    % 
    % hold off
    % 
    %     
    % 
    % %                 % % % % % % % % % % % % % % %      easterly & salt
    % % % % % % % %                 correlation coefficient
    % % % % % % % % % % % % % % % %                 significance level is higher than 95
    % 
    % figdir=[figrawdir,'analysis\'];
    % outfile = strcat(figdir,regionname);
    % if (exist(strcat(figdir) , 'dir') ~= 7)
    %     mkdir(strcat(figdir));
    % end 
    % 
    % close all;
    % m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    % hold on;
    % pval_mask_easterly=pval_data_easterly(:,:);
    % pval_mask_easterly(pval_data_easterly(:,:)<=0.05)=1.0;
    % pval_mask_easterly(pval_data_easterly(:,:)>0.05)=0;
    %     
    % m_pcolor(lon',lat',corr_data_easterly(:,:)' .* pval_mask_easterly');
    % shading(gca,m_pcolor_shading_method);
    % m_gshhs_i('color',m_gshhs_line_color);
    % m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    % 
    % %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
    % titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),', easterly(wind) & salt, corr, sig-lev 95%');
    % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    % 
    % % set colorbar 
    % h = colorbar;
    % colormap(bwrmap);
    % set(h,'fontsize',colorbar_fontsize);
    % %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
    % caxis([-1 1]);
    % 
    % % set grid
    % m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    % 
    % set(gcf, 'PaperUnits', 'points');
    % set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    % set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    % 
    % tifname=strcat(outfile, '_easterly_salt_corr_sig_95_','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    %     num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    % drawnow;
    % saveas(gcf,tifname,'tif');
    % 
    % disp(' ')
    % disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
    % disp(' ')
    % disp([' File path is : ',tifname])
    % disp(' ')
    % 
    % hold off
end

% % % % % 
% % % % % lag correlation plot
% % % % % 
for lag=1:7
%     %                 % % % % % % % % % % % % % % %      southerly & salt
%     % % % % % % % % % % % % % % %                 lag! correlation coeffcient
% 
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
% 
%     m_pcolor(lon',lat',lag_corr_data_southerly(:,:,lag)');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%     titlename = strcat(regionname,', ',calendarname{comb_month(8)},  ...
%         num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ',-',num2str(lag),'day lag southerly(wind) & salt, corr coef');
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([-0.5 0.5]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%     tifname=strcat(outfile, '_southerly_salt_corr_coef','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off
% 
% 
% 
%     %                 % % % % % % % % % % % % % % %      easterly & salt
%     % % % % % % % % % % % % % % %                lag ! correlation coeffcient
% 
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
% 
%     m_pcolor(lon',lat',lag_corr_data_easterly(:,:,lag)');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%     titlename = strcat(regionname,', ',calendarname{comb_month(8)},  ...
%         num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ',-',num2str(lag),'day lag easterly(wind) & salt, corr coef');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([-0.5 0.5]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%     tifname=strcat(outfile, '_easterly_salt_corr_coef','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off
% 
% 
%     %                 % % % % % % % % % % % % % % %   lag   southerly & salt
%     % % % % % % % % % % % % % % %                 probability value (p-val),
%     % % % % % % % % % % % % % %  significant level --> generally 0.05< or 0.01<
% 
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
% 
%     m_pcolor(lon',lat',lag_pval_data_southerly(:,:,lag)');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
% 
%     titlename = strcat(regionname,', ',calendarname{comb_month(8)},  ...
%         num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ',-',num2str(lag),'day lag southerly(wind) & salt, pval');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([0 0.1]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%     tifname=strcat(outfile, '_southerly_salt_pval','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off
% 
% 
%     %                 % % % % % % % % % % % % % % %  lag    easterly & salt
%     % % % % % % % % % % % % % % %                 probability value (p-val),
%     % % % % % % % % % % % % % %  significant level --> generally 0.05< or 0.01<
% 
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
% 
%     m_pcolor(lon',lat',lag_pval_data_easterly(:,:,lag)');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%     titlename = strcat(regionname,', ',calendarname{comb_month(8)},  ...
%         num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ',-',num2str(lag),'day lag easterly(wind) & salt, pval');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([0 0.1]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%     tifname=strcat(outfile, '_easterly_salt_pval','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off



%     %                 % % % % % % % % % % % % % % %   lag   southerly & salt
%     % % % % % % %                 correlation coefficient
%     % % % % % % % % % % % % % % %                 significance level is higher than 95
%     % % % % % % % % % % % % % %  significant level --> generally 0.05< or 0.01<
% 
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
%     lag_pval_mask_southerly=lag_pval_data_southerly(:,:,lag);
%     lag_pval_mask_southerly(lag_pval_data_southerly(:,:,lag)<=0.05)=1.0;
%     lag_pval_mask_southerly(lag_pval_data_southerly(:,:,lag)>0.05)=0;
% 
%     m_pcolor(lon',lat',lag_corr_data_southerly(:,:,lag)' .* lag_pval_mask_southerly');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%     titlename = strcat(regionname,', ',calendarname{comb_month(8)},  ...
%         num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ',-',num2str(lag),'day lag southerly(wind) & salt, corr, sig-lev 95%');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([-0.5 0.5]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%     
%     tifname=strcat(outfile, '_southerly_salt_corr_sig_95','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off




%     %                 % % % % % % % % % % % % % % %   lag   easterly & salt
%     % % % % % % %                 correlation coefficient
%     % % % % % % % % % % % % % % %                 significance level is higher than 95
% 
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
%     lag_pval_mask_easterly=lag_pval_data_easterly(:,:,lag);
%     lag_pval_mask_easterly(lag_pval_data_easterly(:,:,lag)<=0.05)=1.0;
%     lag_pval_mask_easterly(lag_pval_data_easterly(:,:,lag)>0.05)=0;
% 
%     m_pcolor(lon',lat',lag_corr_data_easterly(:,:,lag)' .* lag_pval_mask_easterly');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%    titlename = strcat(regionname,', ',calendarname{comb_month(8)},  ...
%         num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ',-',num2str(lag),'day lag easterly(wind) & salt, corr, sig-lev 95%');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([-0.5 0.5]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%     tifname=strcat(outfile, '_easterly_salt_corr_sig_95','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off

end

    


% % % % % % % get lowest salt(30 or minvalue) latitude along longitude line 123 ~ 128
% % % % % % % and plot
figdir=[figrawdir,'analysis\'];
outfile = strcat(figdir,regionname);
if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 
for lati=123 : 128
    [lsindw, lsinde, lsinds, lsindn] = findind_Y(dl, [lati, lati, 28, 35], lon', lat');  %% lon_whole' and lat_whole' should have [y x] dimension
    lsind_lon=(lsindw+lsinde)/2;
    for t=1:size(comb_data,3)
        j=lsinds;
        while (j<lsindn)
            lat_lowsalt(lati-122,t)=NaN;
%             if(comb_data(lsind_lon,j,t) <= 30)
            if(comb_data(lsind_lon,j,t) <= min(comb_data(lsind_lon,lsinds:lsindn,t)))
                lat_lowsalt(lati-122,t)=lat(lsind_lon,j);
                j=lsindn-1;
            end
            j=j+1;
        end
        t_ax(t) = datenum([num2str(comb_year(t),'%04i'),'-',num2str(comb_month(t),'%02i'),'-',num2str(comb_day(t),'%02i')]);
    end
    year1_startind = datenum([num2str(comb_year(1),'%04i'),'-',num2str(1,'%02i'),'-',num2str(1,'%02i')]);
    year1_endind = datenum([num2str(comb_year(310),'%04i'),'-',num2str(12,'%02i'),'-',num2str(31,'%02i')]);

    % t_ax2=linspace(t_ax(1),t_ax(size(comb_data,3)),310);
    t_ax2=linspace(year1_startind,year1_endind,310);
 
    lowlat_plot=plot(t_ax2,lat_lowsalt(lati-122,:),'b');
    hold on
    for i=32:31:1+31*9
        line([t_ax2(i), t_ax2(i)],[min(lat_lowsalt(lati-122,:)), max(lat_lowsalt(lati-122,:))],'color','r')
    end
    datetick('x','yyyy')
    grid on
    grid minor
    hold off
    xlabel('year')
    ylabel('Salt==min Lowest Latitude (^oN)')
%     titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ', lowest lat (salt==30,',num2str(lati),'^oE)');
    titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
        ', lowest lat (salt==min,',num2str(lati),'^oE)');
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    set(lowlat_plot,'LineWidth',2);
    axis tight
    ylim([28 35]);
    drawnow
%     tifname=strcat(outfile, '_lowlat_30_',num2str(lati),'_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%               num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    tifname=strcat(outfile, '_lowlat_minsalt_',num2str(lati),'_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    saveas(gcf,tifname,'tif');
    close all;
end

% % lowlat_plot=plot(t_ax2,mean(lat_lowsalt,1,'omitnan'),'b');
% % hold on
% % for i=32:31:1+31*9
% %     line([t_ax2(i), t_ax2(i)],[min(lat_lowsalt(lati-122,:)), max(lat_lowsalt(lati-122,:))],'color','r')
% % end
% % datetick('x','yyyy')
% % grid on
% % grid minor
% % hold off
% % xlabel('year')
% % ylabel('Salt==30 Lowest Latitude (^oN)')
% % titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
% %     ', lowest lat (salt==30, latmean',')');
% % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % set(lowlat_plot,'LineWidth',2);
% % axis tight
% % ylim([28 35]);
% % drawnow
% % tifname=strcat(outfile, '_lowlat_30_latnanmean','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
% %             num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
% % saveas(gcf,tifname,'tif');
% % close all;


lowlat_plot=plot(t_ax2,mean(lat_lowsalt,1),'b');
hold on
for i=32:31:1+31*9
    line([t_ax2(i), t_ax2(i)],[min(lat_lowsalt(lati-122,:)), max(lat_lowsalt(lati-122,:))],'color','r')
end
datetick('x','yyyy')
grid on
grid minor
hold off
xlabel('year')
ylabel('Salt==30 Lowest Latitude (^oN)')
titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
    ', lowest lat (salt==min, latmean',')');
title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
set(lowlat_plot,'LineWidth',2);
axis tight
ylim([28 35]);
drawnow
tifname=strcat(outfile, '_lowlat_minsalt_latmean','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
                num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
saveas(gcf,tifname,'tif');
close all;


lat_lowsalt_mean=mean(lat_lowsalt,1);
[lon_min_lowsalt, lon_max_lowsalt, lat_min_lowsalt, lat_max_lowsalt] = findind_Y(dl, [123, 128, 31, 34], lon', lat');  %% lon' and lat' should have [y x] dimension
easterly_lowsalt2 = squeeze(mean(mean(-comb_uwind2(lon_min_lowsalt:lon_max_lowsalt,lat_min_lowsalt:lat_max_lowsalt,:),'omitnan'),'omitnan'));
southerly_lowsalt2 = squeeze(mean(mean(comb_vwind2(lon_min_lowsalt:lon_max_lowsalt,lat_min_lowsalt:lat_max_lowsalt,:),'omitnan'),'omitnan'));
lon_length_lowsalt=lon_max_lowsalt-lon_min_lowsalt+1;
lat_length_lowsalt=lat_max_lowsalt-lat_min_lowsalt+1;

% % % easterly and southerly with lowest salinity scatter plot, correlation
for lag=1:7
    figdir=[figrawdir,'analysis\'];
    outfile = strcat(figdir,regionname);
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
%     for i=1:lon_length_lowsalt
%         for j=1:lat_length_lowsalt
            for year=2009:2018
                southerly_lowsalt2_lag((year-2009)*tempday+1:(year-2008)*tempday) = ...
                    southerly_lowsalt2((year-2009)*(tempday+7)+8-lag:(year-2009)*(tempday+7)+8-lag +tempday-1);
                easterly_lowsalt2_lag((year-2009)*tempday+1:(year-2008)*tempday) = ...
                    easterly_lowsalt2((year-2009)*(tempday+7)+8-lag:(year-2009)*(tempday+7)+8-lag +tempday-1);
            end
%             
% %             [R, P] = corrcoef(comb_data(i,j,:),comb_vwind_lag(i,j,:));
% %             lag_corr_data_southerly(i,j,lag)=R(1,2);
% %             lag_pval_data_southerly(i,j,lag)=P(1,2);
% %             [R, P] = corrcoef(comb_data(i,j,:),-comb_uwind_lag(i,j,:));
% %             lag_corr_data_easterly(i,j,lag)=R(1,2);
% %             lag_pval_data_easterly(i,j,lag)=P(1,2);
%         end
%     end
% % % % easterly, lowlat scatter
    sc = scatter(lat_lowsalt_mean,easterly_lowsalt2_lag, 10, jet(length(lat_lowsalt_mean)), 'filled');
    lsl=lsline;
    lat_lowsalt_easterly_corr_coef = corrcoef(squeeze(lat_lowsalt_mean), squeeze(easterly_lowsalt2_lag), 'Alpha',0.05);
    set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
    yearchar1 = num2str(min(comb_year),'%04i');
    yearchar2 = num2str(max(comb_year),'%04i');
    titlename = strcat(calendarname{comb_month(8)}(1:3),  ...
        yearchar1(3:4),' to ', yearchar2(3:4), ...
        ',-',num2str(lag),' lag Easterly & minS Lat, corr');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    xlabel('Lowest Salt Latitude (^oN)','FontSize',vert_grid_fontsize,'fontweight','bold')
    ylabel('Easterly velocity (m/s)','FontSize',vert_grid_fontsize,'fontweight','bold')
    h = colorbar;
    colormap jet;
    caxis([min(inputyear), max(inputyear)]);
    set(h,'fontsize',colorbar_fontsize-3);
    title(h,'Year','fontsize',colorbar_title_fontsize);
    gcaaa=gca;
    textypos =  gcaaa.YLim(2)-(gcaaa.YLim(2)-gcaaa.YLim(1))/20;
    textxpos =  gcaaa.XLim(2)-(gcaaa.XLim(2)-gcaaa.XLim(1))/4.5;
    text1 = text(textxpos, textypos, ['r = ', num2str(lat_lowsalt_easterly_corr_coef(1,2),'%3.2f')], 'FontSize', colorbar_title_fontsize); 

    tifname=strcat(outfile, '_scatter_easterly_lowlat_corr','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
        num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    drawnow;
    saveas(gcf,tifname,'tif');
    close all;
    
    % % % % southerly, lowlat scatter
    sc = scatter(lat_lowsalt_mean,southerly_lowsalt2_lag, 10, jet(length(lat_lowsalt_mean)), 'filled');
    lsl=lsline;
    lat_lowsalt_southerly_corr_coef = corrcoef(squeeze(lat_lowsalt_mean), squeeze(southerly_lowsalt2_lag), 'Alpha',0.05);
    set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
    yearchar1 = num2str(min(comb_year),'%04i');
    yearchar2 = num2str(max(comb_year),'%04i');
    titlename = strcat(calendarname{comb_month(8)}(1:3),  ...
        yearchar1(3:4),' to ', yearchar2(3:4), ...
        ',-',num2str(lag),' lag Southerly & minS Lat, corr');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    xlabel('Lowest Salt Latitude (^oN)','FontSize',vert_grid_fontsize,'fontweight','bold')
    ylabel('Easterly velocity (m/s)','FontSize',vert_grid_fontsize,'fontweight','bold')
    h = colorbar;
    colormap jet;
    caxis([min(inputyear), max(inputyear)]);
    set(h,'fontsize',colorbar_fontsize-3);
    title(h,'Year','fontsize',colorbar_title_fontsize);
    gcaaa=gca;
    textypos =  gcaaa.YLim(2)-(gcaaa.YLim(2)-gcaaa.YLim(1))/20;
    textxpos =  gcaaa.XLim(2)-(gcaaa.XLim(2)-gcaaa.XLim(1))/4.5;
    text1 = text(textxpos, textypos, ['r = ', num2str(lat_lowsalt_southerly_corr_coef(1,2),'%3.2f')], 'FontSize', colorbar_title_fontsize); 

    tifname=strcat(outfile, '_scatter_southerly_lowlat_corr','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
        num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    drawnow;
    saveas(gcf,tifname,'tif');
    close all;
end

% sc=scatter(mean_isodepth, mean_surf_vel,marker_size,bwr_all, 'filled');
% scatter(lat_lowsalt_mean,southerly_lowsalt)
% sc = scatter(lat_lowsalt_mean,easterly_lowsalt, 4, 'filled');
% lsl=lsline;
% lat_lowsalt_easterly_corr_coef = corrcoef(squeeze(lat_lowsalt_mean), squeeze(easterly_lowsalt), 'Alpha',0.05);






% % test figure for lowest salt latitude plot
% plot(t_ax2,lat_lowsalt,'b');
% hold on
% for i=32:31:1+31*9
%     line([t_ax2(i), t_ax2(i)],[min(lat_lowsalt), max(lat_lowsalt)],'color','r')
% end
% datetick('x','yyyy')
% grid minor
% hold off




% % % % % %   get monthly data
for year=2009:2018
    riv_tr=ncread(['E:\Data\Model\ROMS\nwp_1_10\input\test08\roms_river_nwp_1_10_',num2str(year,'%04i'),'.nc'],'river_transport');
%     comb_riv_tr(12*(year-2009)+1:12*(year-2008))=riv_tr(1,:);
    comb_riv_tr(year-2008,:)=riv_tr(1,:);
    filename = [filedir,num2str(year,'%04i'),'\',testname,'_monthly_',num2str(year,'%04i'),'_08.nc']
    data = ncread(filename,varname,[lon_min(1), lat_min(1), length_depth, 1], [lon_length_sliced, lat_length_sliced, 1, 1]);  % get surface data
    comb_mondata(:,:,year-2008)=data;
    
    mon_uwind = ncread(filename,'Uwind',[lon_min(1), lat_min(1), 1], [lon_length_sliced, lat_length_sliced, 1]);  % get surface data
    mon_vwind = ncread(filename,'Vwind',[lon_min(1), lat_min(1), 1], [lon_length_sliced, lat_length_sliced, 1]);  % get surface data
    mon_u = ncread(filename,'u',[lon_min(1)+1, lat_min(1), length_depth, 1], [lon_length_sliced-1, lat_length_sliced, 1, 1]);  % get surface data
    mon_v = ncread(filename,'v',[lon_min(1), lat_min(1)+1, length_depth, 1], [lon_length_sliced, lat_length_sliced-1, 1, 1]);  % get surface data

%                 same as u2rho_2d in romstools
    mon_u_rho(2:lon_length_sliced-1,:) = 0.5 * (mon_u(1:lon_length_sliced-2,:) + mon_u(2:lon_length_sliced-1,:));
    mon_u_rho(1,:)=mon_u_rho(2,:);
    mon_u_rho(lon_length_sliced,:)=mon_u_rho(lon_length_sliced-1,:);
%                 same as v2rho_2d in romstools
    mon_v_rho(:,2:lat_length_sliced-1) = 0.5 * (mon_v(:,1:lat_length_sliced-2) + mon_v(:,2:lat_length_sliced-1));
    mon_v_rho(:,1) = mon_v_rho(:,2);
    mon_v_rho(:,lat_length_sliced) = mon_v_rho(:,lat_length_sliced-1);

    if (exist('ref_vec_x_range' , 'var') ~= 1)
        ref_vec_x_ind = find(abs(lon(:,1)-m_quiver_ref_text_x_location) == min(abs(lon(:,1)-m_quiver_ref_text_x_location)));
        ref_vec_y_ind = find(abs(lat(1,:)-m_quiver_ref_text_y_location) == min(abs(lat(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval-1;
        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval-1;
    end
    mon_u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
    mon_v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;    
    mon_uwind(ref_vec_x_range,ref_vec_y_range)=5.0;
    mon_vwind(ref_vec_x_range,ref_vec_y_range)=0.00001;    
    
    comb_mon_uwind(:,:,year-2008)=mon_uwind;
    comb_mon_vwind(:,:,year-2008)=mon_vwind;


% monthly current and wind plot (with salt)
    for grouping=1:1
        %             % % % % % % % % % % % % % % %         monthly current
        %         figdir=[figrawdir,'current\'];
        %         outfile = strcat(figdir,regionname);
        %         if (exist(strcat(figdir) , 'dir') ~= 7)
        %             mkdir(strcat(figdir));
        %         end 
        %         
        %         close all;
        %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        %         hold on;
        % 
        %         m_pcolor(lon',lat',comb_mondata(:,:,year-2008)');
        %         shading(gca,m_pcolor_shading_method);
        %         m_gshhs_i('color',m_gshhs_line_color);
        %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
        %                         lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
        %                         mon_u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
        %                         mon_v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
        %                         'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
        % 
        %         titlename = strcat(regionname,', ',calendarname{comb_month(i)},num2str(year,'%04i'),', current');
        %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
        % 
        %         
        %         % set colorbar 
        %         h = colorbar;
        %         colormap(jet_mod);
        %         set(h,'fontsize',colorbar_fontsize);
        % %         title(h,'m/s','fontsize',colorbar_title_fontsize);
        %         caxis(shadlev);
        %         
        %         % set grid
        %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
        %         m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
        % 
        %         set(gcf, 'PaperUnits', 'points');
        %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
        %         
        %         tifname=strcat(outfile,'_current_',num2str(year,'%04i'),'_', ...
        %             num2str(comb_month(i),'%02i'),'_','.tif'); %% ~_year_month.tif
        %         drawnow;
        %         saveas(gcf,tifname,'tif');
        % 
        %         disp(' ')
        %         disp(['clim_', num2str(tempmonth), '_', 'current', ' plot is created.'])
        %         disp(' ')
        %         disp([' File path is : ',tifname])
        %         disp(' ')
        % 
        %         hold off
        % 
        %         
        %                 % % % % % % % % % % % % % % %        monthly wind
        %         figdir=[figrawdir,'wind\'];
        %         outfile = strcat(figdir,regionname);
        %         if (exist(strcat(figdir) , 'dir') ~= 7)
        %             mkdir(strcat(figdir));
        %         end 
        %         
        %         close all;
        %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        %         hold on;
        % 
        %         m_pcolor(lon',lat',comb_mondata(:,:,year-2008)');
        %         shading(gca,m_pcolor_shading_method);
        %         m_gshhs_i('color',m_gshhs_line_color);
        %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
        %                         lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
        %                         mon_uwind(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size/20.0, ...
        %                         mon_vwind(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size/20.0, ...
        %                         'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
        % 
        %         titlename = strcat(regionname,', ',calendarname{comb_month(i)},num2str(year,'%04i'),', wind');
        %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
        % 
        %         
        %         % set colorbar 
        %         h = colorbar;
        %         colormap(jet_mod);
        %         set(h,'fontsize',colorbar_fontsize);
        % %         title(h,'m/s','fontsize',colorbar_title_fontsize);
        %         caxis(shadlev);
        %         
        %         % set grid
        %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
        %         m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, [num2str(5),' m/s'], 'FontSize', m_quiver_ref_text_fontsize); 
        % 
        %         set(gcf, 'PaperUnits', 'points');
        %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
        %         
        %         tifname=strcat(outfile,'_wind_',num2str(year,'%04i'),'_', ...
        %             num2str(comb_month(i),'%02i'),'_','.tif'); %% ~_year_month.tif
        %         drawnow;
        %         saveas(gcf,tifname,'tif');
        % 
        %         disp(' ')
        %         disp(['clim_', num2str(tempmonth), '_', 'wind', ' plot is created.'])
        %         disp(' ')
        %         disp([' File path is : ',tifname])
        %         disp(' ')
        % 
        %         hold off
    end
end

% % % % % %   get monthly data (July)
for year=2009:2018
    riv_tr=ncread(['E:\Data\Model\ROMS\nwp_1_10\input\test08\roms_river_nwp_1_10_',num2str(year,'%04i'),'.nc'],'river_transport');
%     comb_riv_tr(12*(year-2009)+1:12*(year-2008))=riv_tr(1,:);
    comb_riv_tr(year-2008,:)=riv_tr(1,:);
    filename = [filedir,num2str(year,'%04i'),'\',testname,'_monthly_',num2str(year,'%04i'),'_07.nc']
%     data = ncread(filename,varname,[lon_min(1), lat_min(1), length_depth, 1], [lon_length_sliced, lat_length_sliced, 1, 1]);  % get surface data
%     comb_mondata(:,:,year-2008)=data;
    
    mon_uwind_7 = ncread(filename,'Uwind',[lon_min(1), lat_min(1), 1], [lon_length_sliced, lat_length_sliced, 1]);  % get surface data
    mon_vwind_7 = ncread(filename,'Vwind',[lon_min(1), lat_min(1), 1], [lon_length_sliced, lat_length_sliced, 1]);  % get surface data
%     mon_u = ncread(filename,'u',[lon_min(1)+1, lat_min(1), length_depth, 1], [lon_length_sliced-1, lat_length_sliced, 1, 1]);  % get surface data
%     mon_v = ncread(filename,'v',[lon_min(1), lat_min(1)+1, length_depth, 1], [lon_length_sliced, lat_length_sliced-1, 1, 1]);  % get surface data

%                 same as u2rho_2d in romstools
%     mon_u_rho(2:lon_length_sliced-1,:) = 0.5 * (mon_u(1:lon_length_sliced-2,:) + mon_u(2:lon_length_sliced-1,:));
%     mon_u_rho(1,:)=mon_u_rho(2,:);
%     mon_u_rho(lon_length_sliced,:)=mon_u_rho(lon_length_sliced-1,:);
% %                 same as v2rho_2d in romstools
%     mon_v_rho(:,2:lat_length_sliced-1) = 0.5 * (mon_v(:,1:lat_length_sliced-2) + mon_v(:,2:lat_length_sliced-1));
%     mon_v_rho(:,1) = mon_v_rho(:,2);
%     mon_v_rho(:,lat_length_sliced) = mon_v_rho(:,lat_length_sliced-1);

    if (exist('ref_vec_x_range' , 'var') ~= 1)
        ref_vec_x_ind = find(abs(lon(:,1)-m_quiver_ref_text_x_location) == min(abs(lon(:,1)-m_quiver_ref_text_x_location)));
        ref_vec_y_ind = find(abs(lat(1,:)-m_quiver_ref_text_y_location) == min(abs(lat(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval-1;
        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval-1;
    end
%     mon_u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
%     mon_v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;    
    mon_uwind_7(ref_vec_x_range,ref_vec_y_range)=5.0;
    mon_vwind_7(ref_vec_x_range,ref_vec_y_range)=0.00001;    
    
    comb_mon_uwind_7(:,:,year-2008)=mon_uwind_7;
    comb_mon_vwind_7(:,:,year-2008)=mon_vwind_7;

end


% % % % % % % get lowest salt(30 or minvalue) latitude along longitude line
% % % % % % % 123 ~ 128 from monthly data
% % % % % % % and plot
figdir=[figrawdir,'analysis\'];
outfile = strcat(figdir,regionname);
if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 
for lati=123 : 128
    [lsindw, lsinde, lsinds, lsindn] = findind_Y(dl, [lati, lati, 28, 35], lon', lat');  %% lon_whole' and lat_whole' should have [y x] dimension
    lsind_lon=(lsindw+lsinde)/2;
    for t=1:size(comb_mondata,3)
        j=lsinds;
        while (j<lsindn)
            lat_lowsalt_mon(lati-122,t)=NaN;
%             if(comb_data(lsind_lon,j,t) <= 30)
            if(comb_mondata(lsind_lon,j,t) <= min(comb_mondata(lsind_lon,lsinds:lsindn,t)))
                lat_lowsalt_mon(lati-122,t)=lat(lsind_lon,j);
                j=lsindn-1;
            end
            j=j+1;
        end
        t_ax_mon(t) = datenum([num2str(comb_year(t*31),'%04i'),'-',num2str(comb_month(t*31),'%02i'),'-',num2str(comb_day(t*31),'%02i')]);
    end
    year1_startind = datenum([num2str(comb_year(1),'%04i'),'-',num2str(1,'%02i'),'-',num2str(1,'%02i')]);
    year1_endind = datenum([num2str(comb_year(310),'%04i'),'-',num2str(1,'%02i'),'-',num2str(1,'%02i')]);

    % t_ax2=linspace(t_ax(1),t_ax(size(comb_data,3)),310);
    t_ax2_mon=linspace(year1_startind,year1_endind,10);
 
    lowlat_mon_plot=plot(t_ax2_mon,lat_lowsalt_mon(lati-122,:),'b');
    hold on
    for i=1:10
        line([t_ax2_mon(i), t_ax2_mon(i)],[min(lat_lowsalt_mon(lati-122,:)), max(lat_lowsalt_mon(lati-122,:))],'color','r')
    end
    datetick('x','yyyy')
    grid on
    grid minor
    hold off
    xlabel('year')
    ylabel('Salt==min Lowest Latitude (^oN)')
%     titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ', lowest lat (salt==30,',num2str(lati),'^oE)');
    titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
        ', lowest lat-mon (salt==min,',num2str(lati),'^oE)');
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    set(lowlat_mon_plot,'LineWidth',2);
    axis tight
    ylim([28 35]);
    drawnow
%     tifname=strcat(outfile, '_lowlat_30_',num2str(lati),'_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%               num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    tifname=strcat(outfile, '_lowlat_mon_minsalt_',num2str(lati),'_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    saveas(gcf,tifname,'tif');
    close all;
end

% % % % mean lowlat, monthly plot
    lowlat_mon_plot=plot(t_ax2_mon,mean(lat_lowsalt_mon(:,:),1,'omitnan'),'b');
    hold on
%     for i=1:10
%         line([t_ax2_mon(i), t_ax2_mon(i)],[min(lat_lowsalt_mon(:,:)), max(lat_lowsalt_mon(:,:))],'color','r')
%     end
    datetick('x','yyyy')
    grid on
    grid minor
    hold off
    xlabel('year')
    ylabel('Salt==min Lowest Latitude (^oN)')
%     titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ', lowest lat (salt==30,',num2str(lati),'^oE)');
    titlename = strcat(regionname,', ',calendarname{comb_month(8)}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
        ', lowest lat-mon (salt==min,','mean^oE)');
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    set(lowlat_mon_plot,'LineWidth',2);
    set(gca, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize-2)
    axis tight
    ylim([31 34]);
    drawnow
%     tifname=strcat(outfile, '_lowlat_30_',num2str(lati),'_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%               num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    tifname=strcat(outfile, '_lowlat_mon_minsalt_','mean','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    saveas(gcf,tifname,'tif');
    close all;
    
    
lat_lowsalt_mon_mean=mean(lat_lowsalt_mon(:,:),1,'omitnan');
for mon=5:8
    sc = scatter(lat_lowsalt_mon_mean,comb_riv_tr(:,mon), 20, jet(length(lat_lowsalt_mon_mean)), 'filled');
    lsl=lsline;
    lat_lowsalt_riv_corr_coef = corrcoef(squeeze(lat_lowsalt_mon_mean), squeeze(comb_riv_tr(:,mon)), 'Alpha',0.05);
    set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
    yearchar1 = num2str(min(comb_year),'%04i');
    yearchar2 = num2str(max(comb_year),'%04i');
    titlename = strcat(calendarname{mon}(1:3),  ...
        yearchar1(3:4),' to ', yearchar2(3:4), ...
        ',',' riv tr & minS Lat');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    xlabel('Lowest Salt latitude (^oN)','FontSize',vert_grid_fontsize,'fontweight','bold')
    ylabel('Changjiang discharge (m^3/s)','FontSize',vert_grid_fontsize,'fontweight','bold')
    h = colorbar;
    colormap jet;
    caxis([min(inputyear), max(inputyear)]);
    set(h,'fontsize',colorbar_fontsize-3);
    title(h,'Year','fontsize',colorbar_title_fontsize);
    gcaaa=gca;
    textypos =  gcaaa.YLim(2)-(gcaaa.YLim(2)-gcaaa.YLim(1))/20;
    textxpos =  gcaaa.XLim(2)-(gcaaa.XLim(2)-gcaaa.XLim(1))/4.5;
    text1 = text(textxpos, textypos, ['r = ', num2str(lat_lowsalt_riv_corr_coef(1,2),'%3.2f')], 'FontSize', colorbar_title_fontsize); 

    tifname=strcat(outfile, '_scatter_riv_lowlat_corr','_',num2str(mon,'%02i'),'mon_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
        num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    drawnow;
    saveas(gcf,tifname,'tif');
    close all;
end

lat_lowsalt_mean=mean(lat_lowsalt,1);
[lon_min_lowsalt, lon_max_lowsalt, lat_min_lowsalt, lat_max_lowsalt] = findind_Y(dl, [123, 128, 31, 34], lon', lat');  %% lon' and lat' should have [y x] dimension
easterly_mon_lowsalt2 = squeeze(mean(mean(-comb_mon_uwind(lon_min_lowsalt:lon_max_lowsalt,lat_min_lowsalt:lat_max_lowsalt,:),'omitnan'),'omitnan'));
southerly_mon_lowsalt2 = squeeze(mean(mean(comb_mon_vwind(lon_min_lowsalt:lon_max_lowsalt,lat_min_lowsalt:lat_max_lowsalt,:),'omitnan'),'omitnan'));
lon_length_lowsalt=lon_max_lowsalt-lon_min_lowsalt+1;
lat_length_lowsalt=lat_max_lowsalt-lat_min_lowsalt+1;

% % % %  easterly and lowest salinity line
sc = scatter(lat_lowsalt_mon_mean,easterly_mon_lowsalt2, 20, jet(length(lat_lowsalt_mon_mean)), 'filled');
lsl=lsline;
lat_lowsalt_riv_corr_coef = corrcoef(squeeze(lat_lowsalt_mon_mean), squeeze(easterly_mon_lowsalt2), 'Alpha',0.05);
set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
yearchar1 = num2str(min(comb_year),'%04i');
yearchar2 = num2str(max(comb_year),'%04i');
titlename = strcat(calendarname{8}(1:3),  ...
    yearchar1(3:4),' to ', yearchar2(3:4), ...
    ',',' Easterly & minS Lat');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
xlabel('Lowest Salt latitude (^oN)','FontSize',vert_grid_fontsize,'fontweight','bold')
ylabel('Easterly Velocity (m/s)','FontSize',vert_grid_fontsize,'fontweight','bold')
h = colorbar;
colormap jet;
caxis([min(inputyear), max(inputyear)]);
set(h,'fontsize',colorbar_fontsize-3);
title(h,'Year','fontsize',colorbar_title_fontsize);
gcaaa=gca;
textypos =  gcaaa.YLim(2)-(gcaaa.YLim(2)-gcaaa.YLim(1))/20;
textxpos =  gcaaa.XLim(2)-(gcaaa.XLim(2)-gcaaa.XLim(1))/4.5;
text1 = text(textxpos, textypos, ['r = ', num2str(lat_lowsalt_riv_corr_coef(1,2),'%3.2f')], 'FontSize', colorbar_title_fontsize); 

tifname=strcat(outfile, '_scatter_easterly_lowlat_mon_corr','_',num2str(mon,'%02i'),'mon_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
drawnow;
saveas(gcf,tifname,'tif');
close all;

% % % %  southerly and lowest salinity line
sc = scatter(lat_lowsalt_mon_mean,southerly_mon_lowsalt2, 20, jet(length(lat_lowsalt_mon_mean)), 'filled');
lsl=lsline;
lat_lowsalt_riv_corr_coef = corrcoef(squeeze(lat_lowsalt_mon_mean), squeeze(southerly_mon_lowsalt2), 'Alpha',0.05);
set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
yearchar1 = num2str(min(comb_year),'%04i');
yearchar2 = num2str(max(comb_year),'%04i');
titlename = strcat(calendarname{8}(1:3),  ...
    yearchar1(3:4),' to ', yearchar2(3:4), ...
    ',',' Southerly & minS Lat');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
xlabel('Lowest Salt latitude (^oN)','FontSize',vert_grid_fontsize,'fontweight','bold')
ylabel('Southerly Velocity (m/s)','FontSize',vert_grid_fontsize,'fontweight','bold')
h = colorbar;
colormap jet;
caxis([min(inputyear), max(inputyear)]);
set(h,'fontsize',colorbar_fontsize-3);
title(h,'Year','fontsize',colorbar_title_fontsize);
gcaaa=gca;
textypos =  gcaaa.YLim(2)-(gcaaa.YLim(2)-gcaaa.YLim(1))/20;
textxpos =  gcaaa.XLim(2)-(gcaaa.XLim(2)-gcaaa.XLim(1))/4.5;
text1 = text(textxpos, textypos, ['r = ', num2str(lat_lowsalt_riv_corr_coef(1,2),'%3.2f')], 'FontSize', colorbar_title_fontsize); 

tifname=strcat(outfile, '_scatter_southerly_lowlat_mon_corr','_',num2str(mon,'%02i'),'mon_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
drawnow;
saveas(gcf,tifname,'tif');
close all;






easterly_mon_lowsalt2_7 = squeeze(mean(mean(-comb_mon_uwind_7(lon_min_lowsalt:lon_max_lowsalt,lat_min_lowsalt:lat_max_lowsalt,:),'omitnan'),'omitnan'));
southerly_mon_lowsalt2_7 = squeeze(mean(mean(comb_mon_vwind_7(lon_min_lowsalt:lon_max_lowsalt,lat_min_lowsalt:lat_max_lowsalt,:),'omitnan'),'omitnan'));

% % % %  easterly and lowest salinity line  (July wind)
sc = scatter(lat_lowsalt_mon_mean,easterly_mon_lowsalt2_7, 20, jet(length(lat_lowsalt_mon_mean)), 'filled');
lsl=lsline;
lat_lowsalt_riv_corr_coef = corrcoef(squeeze(lat_lowsalt_mon_mean), squeeze(easterly_mon_lowsalt2_7), 'Alpha',0.05);
set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
yearchar1 = num2str(min(comb_year),'%04i');
yearchar2 = num2str(max(comb_year),'%04i');
titlename = strcat(calendarname{8}(1:3),  ...
    yearchar1(3:4),' to ', yearchar2(3:4), ...
    ',',' Easterly(Jul) & minS Lat');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
xlabel('Lowest Salt latitude (^oN)','FontSize',vert_grid_fontsize,'fontweight','bold')
ylabel('Easterly Velocity (m/s)','FontSize',vert_grid_fontsize,'fontweight','bold')
h = colorbar;
colormap jet;
caxis([min(inputyear), max(inputyear)]);
set(h,'fontsize',colorbar_fontsize-3);
title(h,'Year','fontsize',colorbar_title_fontsize);
gcaaa=gca;
textypos =  gcaaa.YLim(2)-(gcaaa.YLim(2)-gcaaa.YLim(1))/20;
textxpos =  gcaaa.XLim(2)-(gcaaa.XLim(2)-gcaaa.XLim(1))/4.5;
text1 = text(textxpos, textypos, ['r = ', num2str(lat_lowsalt_riv_corr_coef(1,2),'%3.2f')], 'FontSize', colorbar_title_fontsize); 

tifname=strcat(outfile, '_scatter_easterly_jul_lowlat_mon_corr','_',num2str(mon,'%02i'),'mon_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
drawnow;
saveas(gcf,tifname,'tif');
close all;

% % % %  southerly and lowest salinity line (July wind)
sc = scatter(lat_lowsalt_mon_mean,southerly_mon_lowsalt2_7, 20, jet(length(lat_lowsalt_mon_mean)), 'filled');
lsl=lsline;
lat_lowsalt_riv_corr_coef = corrcoef(squeeze(lat_lowsalt_mon_mean), squeeze(southerly_mon_lowsalt2_7), 'Alpha',0.05);
set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
yearchar1 = num2str(min(comb_year),'%04i');
yearchar2 = num2str(max(comb_year),'%04i');
titlename = strcat(calendarname{8}(1:3),  ...
    yearchar1(3:4),' to ', yearchar2(3:4), ...
    ',',' Southerly(Jul) & minS Lat');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
xlabel('Lowest Salt latitude (^oN)','FontSize',vert_grid_fontsize,'fontweight','bold')
ylabel('Southerly Velocity (m/s)','FontSize',vert_grid_fontsize,'fontweight','bold')
h = colorbar;
colormap jet;
caxis([min(inputyear), max(inputyear)]);
set(h,'fontsize',colorbar_fontsize-3);
title(h,'Year','fontsize',colorbar_title_fontsize);
gcaaa=gca;
textypos =  gcaaa.YLim(2)-(gcaaa.YLim(2)-gcaaa.YLim(1))/20;
textxpos =  gcaaa.XLim(2)-(gcaaa.XLim(2)-gcaaa.XLim(1))/4.5;
text1 = text(textxpos, textypos, ['r = ', num2str(lat_lowsalt_riv_corr_coef(1,2),'%3.2f')], 'FontSize', colorbar_title_fontsize); 

tifname=strcat(outfile, '_scatter_southerly_jul_lowlat_mon_corr','_',num2str(mon,'%02i'),'mon_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
    num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
drawnow;
saveas(gcf,tifname,'tif');
close all;


% % % % % % % % depth plot
% % % % % % % % 
% nwp_depth = ncread(filename,'h');  % get surface data
% figdir=[figrawdir,'analysis\'];
% outfile = strcat(figdir,regionname);
% if (exist(strcat(figdir) , 'dir') ~= 7)
%     mkdir(strcat(figdir));
% end 
% 
% close all;
% m_proj(m_proj_name,'lon',[lon_whole(1,1)-0.5 lon_whole(end,end)+0.5],'lat',[lat_whole(1,1)-0.5 lat_whole(end,end)+0.5]);
% hold on;
% 
% m_pcolor(lon_whole',lat_whole',nwp_depth');
% shading(gca,m_pcolor_shading_method);
% m_gshhs_i('color',m_gshhs_line_color);
% m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
% %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
% titlename = strcat('North Western Pacific 1/10^o Model Depth');
% title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
% % set colorbar 
% h = colorbar;
% colormap(jet_mod);
% set(h,'fontsize',colorbar_fontsize);
% %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
% % caxis([-1.0 1.0]);
% 
% % set grid
% m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
% tifname=strcat(outfile, '_depth_','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%     num2str(mon,'%02i'),'.tif'); %% ~_year_month.tif
% drawnow;
% saveas(gcf,tifname,'tif');
% 
% disp(' ')
% disp(['depth', ' plot is created.'])
% disp(' ')
% disp([' File path is : ',tifname])
% disp(' ')
% 
% hold off
% 
% ecs_crd_depth= ncread(filename,'h',[lon_min(1), lat_min(1)], [lon_length_sliced, lat_length_sliced]);  % get surface data
% figdir=[figrawdir,'analysis\'];
% outfile = strcat(figdir,regionname);
% if (exist(strcat(figdir) , 'dir') ~= 7)
%     mkdir(strcat(figdir));
% end 
% 
% close all;
% m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% hold on;
% 
% m_pcolor(lon',lat',ecs_crd_depth');
% shading(gca,m_pcolor_shading_method);
% m_gshhs_i('color',m_gshhs_line_color);
% m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
% %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
% titlename = strcat(regionname,' 1/10^o Model Depth');
% title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
% % set colorbar 
% h = colorbar;
% colormap(jet_mod);
% set(h,'fontsize',colorbar_fontsize);
% %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
% caxis([0 100]);
% 
% % set grid
% m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
% tifname=strcat(outfile, '_ecs_crd_depth_','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%     num2str(mon,'%02i'),'.tif'); %% ~_year_month.tif
% drawnow;
% saveas(gcf,tifname,'tif');
% 
% disp(' ')
% disp(['depth', ' plot is created.'])
% disp(' ')
% disp([' File path is : ',tifname])
% disp(' ')
% 
% hold off




% comb_riv_tr
for mon=1:12
    for i=1:lon_length_sliced
        for j=1:lat_length_sliced
            [R, P] = corrcoef(comb_mondata(i,j,:),comb_riv_tr(:,mon));
            corr_data_transport(i,j,mon)=R(1,2);
            pval_data_transport(i,j,mon)=P(1,2);
        end
    end
end


%                 % % % % % % % % % % % % % % %      riv transport & salt
% % % % % % % % % % % % % % %                 correlation coeffcient
% % % % % % % % % % % % % % %                 pval

for mon = 1:12
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
% 
%     m_pcolor(lon',lat',corr_data_transport(:,:,mon)');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%     %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
%     titlename = strcat(regionname,', ', num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),',',calendarname{mon},' transport & salt(aug), corr coef');
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([-1.0 1.0]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%     tifname=strcat(outfile, '_transport_salt_corr_coef','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(mon,'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off
    
 
%     tr_plot=plot(2009:2018,comb_riv_tr(:,mon),'b');
%     hold on
% %     for i=32:31:1+31*9
% %         line([t_ax2(i), t_ax2(i)],[min(lat_lowsalt(lati-122,:)), max(lat_lowsalt(lati-122,:))],'color','r')
% %     end
% %     datetick('x','yyyy')
%     grid on
%     grid minor
%     hold off
%     xlabel('year')
%     ylabel('River transport')
%     titlename = strcat(regionname,', ',calendarname{mon}, num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
%         ', riv transport',')');
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
%     set(tr_plot,'LineWidth',2);
%     axis tight
%     ylim([10000 70000]);
%     drawnow
%     tifname=strcat(outfile, '_riv_transport','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%     num2str(mon,'%02i'),'.tif'); %% ~_year_month.tif
%     saveas(gcf,tifname,'tif');
%     close all;

% % % % % % % % % % % % tr & salt, pval
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
% 
%     m_pcolor(lon',lat',pval_data_transport(:,:,mon)');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%     %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
%     titlename = strcat(regionname,', ', num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),',',calendarname{mon},' transport & salt(aug), corr coef');
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([0 0.10]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%     tifname=strcat(outfile, '_transport_salt_pval','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(mon,'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off

%     
% % % % % % % % % % % % tr & salt, pval * coef
%     figdir=[figrawdir,'analysis\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
%     pval_mask=pval_data_transport(:,:,mon);
%     pval_mask(pval_data_transport(:,:,mon)<=0.05)=1.0;
%     pval_mask(pval_data_transport(:,:,mon)>0.05)=0;
% 
%     close all;
%     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
% 
%     m_pcolor(lon',lat',corr_data_transport(:,:,mon)' .* pval_mask');
%     shading(gca,m_pcolor_shading_method);
%     m_gshhs_i('color',m_gshhs_line_color);
%     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%     %         titlename = strcat(regionname,', ',calendarname{k},', Salt, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
%     titlename = strcat(regionname,', ', num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'),',',calendarname{mon},' transport & salt(aug), 95% significant corr coef');
%     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%     % set colorbar 
%     h = colorbar;
%     colormap(bwrmap);
%     set(h,'fontsize',colorbar_fontsize);
%     %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
%     caxis([-1 1]);
% 
%     % set grid
%     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%     tifname=strcat(outfile, '_transport_salt_95_sig_corr_coef','_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
%         num2str(mon,'%02i'),'.tif'); %% ~_year_month.tif
%     drawnow;
%     saveas(gcf,tifname,'tif');
% 
%     disp(' ')
%     disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
%     disp(' ')
%     disp([' File path is : ',tifname])
%     disp(' ')
% 
%     hold off



end


% comb_mon_uwind

for i=1:lon_length_sliced
    for j=1:lat_length_sliced
        [R, P] = corrcoef(comb_mondata(i,j,:),comb_mon_vwind(i,j,:));
        corr_data_southerly_mon(i,j)=R(1,2);
        pval_data_southerly_mon(i,j)=P(1,2);
        [R, P] = corrcoef(comb_mondata(i,j,:),-comb_mon_uwind(i,j,:));
        corr_data_easterly_mon(i,j)=R(1,2);
        pval_data_easterly_mon(i,j)=P(1,2);
    end
end

%                 % % % % % % % % % % % % % % %   monthly  easterly & salt
    % % % % % % %                 correlation coefficient
    % % % % % % % % % % % % % % %                 significance level is higher than 95

    figdir=[figrawdir,'analysis\'];
    outfile = strcat(figdir,regionname);
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 

    close all;
    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    pval_mask_easterly_mon=pval_data_easterly(:,:);
    pval_mask_easterly_mon(pval_data_easterly(:,:)<=0.05)=1.0;
    pval_mask_easterly_mon(pval_data_easterly(:,:)>0.05)=0;

    m_pcolor(lon',lat',corr_data_easterly_mon(:,:)' .* pval_mask_easterly_mon');
    shading(gca,m_pcolor_shading_method);
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

   titlename = strcat(regionname,', ',calendarname{comb_month(8)},  ...
        num2str(min(comb_year),'%04i'),' to ', num2str(max(comb_year),'%04i'), ...
        ',-','monthly easterly(wind) & salt, corr, sig-lev 95%');    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(bwrmap);
    set(h,'fontsize',colorbar_fontsize);
    %         title(h,'^oC/y','fontsize',colorbar_title_fontsize);
    caxis([-1 1]);

    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    tifname=strcat(outfile, '_mon_easterly_salt_corr_sig_95','_',num2str(lag,'%04i'),'day_lag_',num2str(min(comb_year),'%04i'),'_', num2str(max(comb_year)),'_', ...
        num2str(comb_month(i),'%02i'),'.tif'); %% ~_year_month.tif
    drawnow;
    saveas(gcf,tifname,'tif');

    disp(' ')
    disp(['clim_', num2str(tempmonth), '_', 'salt', ' plot is created.'])
    disp(' ')
    disp([' File path is : ',tifname])
    disp(' ')

    hold off
    
    close all;


end