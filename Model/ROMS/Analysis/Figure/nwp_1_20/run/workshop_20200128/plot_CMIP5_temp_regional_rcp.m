close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test57', 'test58','test59' ,'test60'};
all_testname2 = {'test57', 'test58', 'test59', 'test60'};

all_region2 ={'AKP2'}
% all_region2 ={'NWP'}

% all_region2 ={'YS'};

% all_var2 = {'SST', 'SSS', 'SSH', 'BT'};
all_var2 = {'SST'};

% all_region2 ={'NWP'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2
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
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
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
        abstrendlev =[4 7];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.3 0.3];
        trendplotlev = [3 7];
        sshlev =[-0.7 1.3];
        sshdifflev = [40 70];

        % for snu_desktop
        testname=all_testname2{testnameind2}    % % need to change
%         inputyear1 = [2006:2015]; % % put year which you want to plot [year year ...]
%         inputyear2 = [2066:2075]; % % put year which you want to plot [year year ...]
%         inputyear1 = [2006:2015]; % % put year which you want to plot [year year ...]
        diffyear = [2020, 2100];
        inputyear2 = [2040:2040]; % % put year which you want to plot [year year ...]
        uvyear =[2020];
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

%         varname ='zeta'
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
            case('SS') %% South Sea
                refpolygon=sspolygon;
            case('YS') %% Yellow Sea
                refpolygon=yspolygon;
            case('ECS') %% East China Sea
                refpolygon=ecspolygon;
            case('AKP') %% Around Korea Peninsula
                refpolygon=akppolygon;
            case('AKP2') %% Around Korea Peninsula
                refpolygon=akp2polygon;
            case('CA') %% Around Korea Peninsula
                refpolygon=capolygon;
            case('EKB') %% Around Korea Peninsula
                refpolygon=akp2polygon;
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

%         load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
%         filename = ['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
            filedir = strcat('G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        run(param_script);

        figdir=[figrawdir,'CLIM\workshop\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);


% start-------------------- decadal SST plot
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', var,'_',num2str(min(inputyear2),'%04i'), ...
                '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear2)
                    tempyear=inputyear2(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_rho=ncread(filename, 'lon_rho');
                            lat_rho=ncread(filename, 'lat_rho');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, varname); 
                        
                        if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        elseif (strcmp(var,'SSH')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        elseif (strcmp(var,'BT')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        end
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(inputyear2) * length(inputmonth)));
                    end
                end

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
            
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;

                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                % contour
                [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', conlev, m_contour_color, 'linewidth', m_contour_linewidth);
                clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' mean, ',testname,',(',num2str(max(inputyear2),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
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
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
        end
% end-------------------- decadal SST plot



% start-------------------- decadal SST diff plot
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_diff_', var,'_',num2str(min(diffyear),'%04i'), ...
                '_',num2str(max(diffyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                tempyear=diffyear(1);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');
                    filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                    if (exist('lon_rho' , 'var') ~= 1)
                        lon_rho=ncread(filename, 'lon_rho');
                        lat_rho=ncread(filename, 'lat_rho');
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end
                    data_info = ncinfo(filename, varname); 

                    if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    elseif (strcmp(var,'SSH')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                    elseif (strcmp(var,'BT')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    end
                    if (exist('mean_data' , 'var') ~= 1)
                        mean_data=zeros(size(data));
                    end
                    mean_data=mean_data + (data / (length(diffyear)/2 * length(inputmonth)));
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                mean_data_early=mean_data;
                clear mean_data
                tempyear=diffyear(2);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');
                    filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                    if (exist('lon_rho' , 'var') ~= 1)
                        lon_rho=ncread(filename, 'lon_rho');
                        lat_rho=ncread(filename, 'lat_rho');
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end
                    data_info = ncinfo(filename, varname); 

                    if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    elseif (strcmp(var,'SSH')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                    elseif (strcmp(var,'BT')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    end
                    if (exist('mean_data' , 'var') ~= 1)
                        mean_data=zeros(size(data));
                    end
                    mean_data=mean_data + (data / (length(diffyear)/2 * length(inputmonth)));
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                mean_data_later=mean_data;
                clear mean_data
                
                mean_data=mean_data_later - mean_data_early;
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
            
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;

                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' diff, ',testname,',(',num2str(max(diffyear),'%04i'),'-', num2str(min(diffyear),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
                title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
                caxis([-5 5]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
        end
% end-------------------- decadal SST diff plot



% start-------------------- decadal SST and current plot
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_uv_', var,'_',num2str(min(uvyear),'%04i'), ...
                '_',num2str(max(uvyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear2)
                    tempyear=uvyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_rho=ncread(filename, 'lon_rho');
                            lat_rho=ncread(filename, 'lat_rho');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, 'u'); 
                        u = ncread(filename,'u',[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1) lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        v = ncread(filename,'v',[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1) 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        
                        if (exist('mean_u' , 'var') ~= 1)
                            mean_u=zeros(size(u));
                            mean_v=zeros(size(v));
                        end
                        mean_u=mean_u + (u / (length(uvyear) * length(inputmonth)));
                        mean_v=mean_v + (v / (length(uvyear) * length(inputmonth)));
                    
                        data_info = ncinfo(filename, varname); 
                        
                        if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        elseif (strcmp(var,'SSH')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        elseif (strcmp(var,'BT')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        end
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(uvyear) * length(inputmonth)));
                    end
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                u_rho = u2rho_2d(mean_u')';
                v_rho = v2rho_2d(mean_v')';
                
                u_rho = u_rho .*mask_model;
                v_rho = v_rho .*mask_model;
                
                if (exist('ref_vec_x_range' , 'var') ~= 1)
                    ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                    ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                    ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                    ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                end
                u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
                
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                % contour
                [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', conlev, m_contour_color, 'linewidth', m_contour_linewidth);
                clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' mean, ',testname,',(',num2str(max(uvyear),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
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
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data mean_u
            end
        end
% end-------------------- decadal SST and current plot


        inputmonth = [2]; % % put month which you want to plot [month month ...]
% start-------------------- decadal SST plot winter
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_Feb_', var,'_',num2str(min(inputyear2),'%04i'), ...
                '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear2)
                    tempyear=inputyear2(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_rho=ncread(filename, 'lon_rho');
                            lat_rho=ncread(filename, 'lat_rho');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, varname); 
                        
                        if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        elseif (strcmp(var,'SSH')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        elseif (strcmp(var,'BT')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        end
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(inputyear2) * length(inputmonth)));
                    end
                end

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
            
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;

                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                % contour
                [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', conlev, m_contour_color, 'linewidth', m_contour_linewidth);
                clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' mean, ',testname,',(',num2str(max(inputyear2),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
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
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
        end
% end-------------------- decadal SST plot winter


% start-------------------- decadal SST diff plot winter
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_diff_Feb_', var,'_',num2str(min(diffyear),'%04i'), ...
                '_',num2str(max(diffyear),'%04i'), '.tif'); %% ~_year_month.jpg
%             if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                tempyear=diffyear(1);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');
                    filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                    if (exist('lon_rho' , 'var') ~= 1)
                        lon_rho=ncread(filename, 'lon_rho');
                        lat_rho=ncread(filename, 'lat_rho');
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end
                    data_info = ncinfo(filename, varname); 

                    if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    elseif (strcmp(var,'SSH')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                    elseif (strcmp(var,'BT')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    end
                    if (exist('mean_data' , 'var') ~= 1)
                        mean_data=zeros(size(data));
                    end
                    mean_data=mean_data + (data / (length(diffyear)/2 * length(inputmonth)));
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                mean_data_early=mean_data;
                clear mean_data
                tempyear=diffyear(2);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');
                    filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                    if (exist('lon_rho' , 'var') ~= 1)
                        lon_rho=ncread(filename, 'lon_rho');
                        lat_rho=ncread(filename, 'lat_rho');
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end
                    data_info = ncinfo(filename, varname); 

                    if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    elseif (strcmp(var,'SSH')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                    elseif (strcmp(var,'BT')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    end
                    if (exist('mean_data' , 'var') ~= 1)
                        mean_data=zeros(size(data));
                    end
                    mean_data=mean_data + (data / (length(diffyear)/2 * length(inputmonth)));
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                mean_data_later=mean_data;
                clear mean_data
                
                mean_data=mean_data_later - mean_data_early;
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
            
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;

                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' diff, ',testname,',(',num2str(max(diffyear),'%04i'),'-', num2str(min(diffyear),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
                title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
                caxis([-5 5]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
%             end
        end
% end-------------------- decadal SST diff plot winter

% start-------------------- decadal SST and current plot winter
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_Feb_uv_', var,'_',num2str(min(uvyear),'%04i'), ...
                '_',num2str(max(uvyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear2)
                    tempyear=uvyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_rho=ncread(filename, 'lon_rho');
                            lat_rho=ncread(filename, 'lat_rho');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, 'u'); 
                        u = ncread(filename,'u',[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1) lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        v = ncread(filename,'v',[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1) 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        
                        if (exist('mean_u' , 'var') ~= 1)
                            mean_u=zeros(size(u));
                            mean_v=zeros(size(v));
                        end
                        mean_u=mean_u + (u / (length(uvyear) * length(inputmonth)));
                        mean_v=mean_v + (v / (length(uvyear) * length(inputmonth)));
                    
                        data_info = ncinfo(filename, varname); 
                        
                        if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        elseif (strcmp(var,'SSH')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        elseif (strcmp(var,'BT')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        end
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(uvyear) * length(inputmonth)));
                    end
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                u_rho = u2rho_2d(mean_u')';
                v_rho = v2rho_2d(mean_v')';
                
                u_rho = u_rho .*mask_model;
                v_rho = v_rho .*mask_model;
                
                if (exist('ref_vec_x_range' , 'var') ~= 1)
                    ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                    ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                    ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                    ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                end
                u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
                
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                % contour
                [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', conlev, m_contour_color, 'linewidth', m_contour_linewidth);
                clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' mean, ',testname,',(',num2str(max(uvyear),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
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
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data mean_u
            end
        end
% end-------------------- decadal SST and current plot winter

        inputmonth = [8]; % % put month which you want to plot [month month ...]
% start-------------------- decadal SST plot summer
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_Aug_', var,'_',num2str(min(inputyear2),'%04i'), ...
                '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear2)
                    tempyear=inputyear2(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_rho=ncread(filename, 'lon_rho');
                            lat_rho=ncread(filename, 'lat_rho');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, varname); 
                        
                        if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        elseif (strcmp(var,'SSH')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        elseif (strcmp(var,'BT')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        end
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(inputyear2) * length(inputmonth)));
                    end
                end

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
            
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;

                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                % contour
                [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', conlev, m_contour_color, 'linewidth', m_contour_linewidth);
                clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' mean, ',testname,',(',num2str(max(inputyear2),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
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
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
        end
% end-------------------- decadal SST plot summer


% start-------------------- decadal SST diff plot summer
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_diff_Aug_', var,'_',num2str(min(diffyear),'%04i'), ...
                '_',num2str(max(diffyear),'%04i'), '.tif'); %% ~_year_month.jpg
%             if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                tempyear=diffyear(1);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');
                    filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                    if (exist('lon_rho' , 'var') ~= 1)
                        lon_rho=ncread(filename, 'lon_rho');
                        lat_rho=ncread(filename, 'lat_rho');
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end
                    data_info = ncinfo(filename, varname); 

                    if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    elseif (strcmp(var,'SSH')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                    elseif (strcmp(var,'BT')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    end
                    if (exist('mean_data' , 'var') ~= 1)
                        mean_data=zeros(size(data));
                    end
                    mean_data=mean_data + (data / (length(diffyear)/2 * length(inputmonth)));
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                mean_data_early=mean_data;
                clear mean_data
                tempyear=diffyear(2);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');
                    filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                    if (exist('lon_rho' , 'var') ~= 1)
                        lon_rho=ncread(filename, 'lon_rho');
                        lat_rho=ncread(filename, 'lat_rho');
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end
                    data_info = ncinfo(filename, varname); 

                    if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    elseif (strcmp(var,'SSH')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                    elseif (strcmp(var,'BT')==1)
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    end
                    if (exist('mean_data' , 'var') ~= 1)
                        mean_data=zeros(size(data));
                    end
                    mean_data=mean_data + (data / (length(diffyear)/2 * length(inputmonth)));
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                mean_data_later=mean_data;
                clear mean_data
                
                mean_data=mean_data_later - mean_data_early;
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
            
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;

                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' diff, ',testname,',(',num2str(max(diffyear),'%04i'),'-', num2str(min(diffyear),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
                title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
                caxis([-5 5]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
%             end
        end
% end-------------------- decadal SST diff plot summer


% start-------------------- decadal SST and current plot summer
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_Aug_uv_', var,'_',num2str(min(uvyear),'%04i'), ...
                '_',num2str(max(uvyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear2)
                    tempyear=uvyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_rho=ncread(filename, 'lon_rho');
                            lat_rho=ncread(filename, 'lat_rho');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, 'u'); 
                        u = ncread(filename,'u',[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1) lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        v = ncread(filename,'v',[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1) 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        
                        if (exist('mean_u' , 'var') ~= 1)
                            mean_u=zeros(size(u));
                            mean_v=zeros(size(v));
                        end
                        mean_u=mean_u + (u / (length(uvyear) * length(inputmonth)));
                        mean_v=mean_v + (v / (length(uvyear) * length(inputmonth)));
                    
                        data_info = ncinfo(filename, varname); 
                        
                        if (strcmp(var,'SST')==1 || strcmp(var,'SSS')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        elseif (strcmp(var,'SSH')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        elseif (strcmp(var,'BT')==1)
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        end
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(uvyear) * length(inputmonth)));
                    end
                end
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data=mean_data.*mask_model;
                
                u_rho = u2rho_2d(mean_u')';
                v_rho = v2rho_2d(mean_v')';
                
                u_rho = u_rho .*mask_model;
                v_rho = v_rho .*mask_model;
                
                if (exist('ref_vec_x_range' , 'var') ~= 1)
                    ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                    ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                    ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                    ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                end
                u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
                
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,m_pcolor_shading_method);   
                
                % contour
                [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', conlev, m_contour_color, 'linewidth', m_contour_linewidth);
                clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(var, ' mean, ',testname,',(',num2str(max(uvyear),'%04i'),') M=', num2str(round(mean(mean_data(:),'omitnan'),2)));  %% + glacier contribution
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
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data mean_u
            end
        end
% end-------------------- decadal SST and current plot summer
    end
end