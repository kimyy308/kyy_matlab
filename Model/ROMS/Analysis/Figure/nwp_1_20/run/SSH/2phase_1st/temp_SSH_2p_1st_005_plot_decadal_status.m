close all; clear all;  clc;
warning off;

% if(isempty(gcp('nocreate')))
%     parpool(4);
% end
all_testname2 = {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};

% all_region2 ={'NWP', 'YS', 'AKP2'}
% all_region2 ={'NWP', 'AKP4'};
% all_region2 ={'NWP', 'AKP4'};
% all_region2 ={'AKP4'};

all_region2 ={'YS'};

% all_var2 = {'SST', 'SSH', 'SSS'};
% all_var2 = {'SSH'};
all_var2 = {'SST'};
% all_var2 = {'SSS'};
% all_var2 = {'BT'};

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
            dropboxpath='C:\users\user\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
            addpath(genpath([dropboxpath '\source\matlab\function']));
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
%         inputyear1 = [1993:2014]; % % put year which you want to plot [year year ...]
        inputyear1 = [1985:2014]; % % put year which you want to plot [year year ...]
%         inputyear1 = [1993:2014]; % % put year which you want to plot [year year ...]
        
%         inputyear2 = [1985:2014]; % % put year which you want to plot [year year ...]
%         inputmonth = [1:12]; % % put month which you want to plot [month month ...]
%         inputmonth = [8]; % % put month which you want to plot [month month ...]
          inputmonth = [2]; % % put month which you want to plot [month month ...]

        for folding=1:1
            fig_flags{1,1}='earlier decadal current plot';
            fig_flags{2,1}='later decadal current plot';
            fig_flags{3,1}='earlier decadal SST, SSS plot';
            fig_flags{4,1}='later decadal SST, SSS plot';
            fig_flags{5,1}='earlier decadal YSBCW plot';
            fig_flags{6,1}='later decadal YSBCW plot';
        end
        for flagi=1:6
            fig_flags{flagi,2}=0;
        end
        fig_flags{1,2}=0;
        fig_flags{2,2}=0;
        fig_flags{3,2}=1;
        fig_flags{4,2}=0;
        fig_flags{5,2}=0;
        fig_flags{6,2}=0;
        
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
            case('AKP4') %% Around Korea Peninsula
                refpolygon=akp4polygon;
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

%         load(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
%         filename = ['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\user\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);
        
        switch testname
            case {'test61', 'test62', 'test63', 'test64'}
                drivename='H:\';
            case {'test57', 'test58', 'test59', 'test60'}
                drivename='I:\';
            case {'test65', 'test66', 'test67', 'test68'}
                drivename='G:\';
            case {'ens08', 'ens09', 'ens10'}
                drivename='E:\';
            case {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
                drivename='E:\';
        end
        
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             figrawdir =strcat('C:\Users\User\Desktop\6th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
%             filedir = strcat(drivename,'Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            filedir = strcat(drivename,'Data\Model\ROMS\nwp_1_20\output\', testname, '\run\packed_monthly\'); % % where data files are          
            matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\mean\');
            griddir = strcat('E:\Data\Model\ROMS\nwp_1_20\output\', testname, '\run\packed_monthly\'); % % where data files are            
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        run(param_script);
        
        if (exist(strcat(matdir) , 'dir') ~= 7)
            mkdir(strcat(matdir));
        end 
        
        figdir=[figrawdir,'CLIM\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        
%         gridname = strcat(griddir, num2str(tempyear,'%04i'), '\', ...
%                                 'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
        
        
% start-------------------- earlier decadal current plot
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_uv_',num2str(min(inputyear1),'%04i'), ...
                '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
            variable='UV';
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
                if (exist(matname , 'file') ~= 2 || fig_flag==2)
                    for yearij=1:length(inputyear1)
                        tempyear=inputyear1(yearij);
                        yearstr=num2str(tempyear, '%04i');
                        for monthij=1:length(inputmonth)
                            tempmonth=inputmonth(monthij);
                            monthstr=num2str(tempmonth, '%02i');
                            
                            filename=[filedir, yearstr, filesep,'pck_', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                            if (exist('lon_rho' , 'var') ~= 1)
                                lon_rho=ncread(filename, 'lon_rho');
                                lat_rho=ncread(filename, 'lat_rho');
                                lon_u=ncread(filename, 'lon_u');
                                lat_u=ncread(filename, 'lat_u');
                                lon_v=ncread(filename, 'lon_v');
                                lat_v=ncread(filename, 'lat_v');
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                                [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(1/20, lonlat(1:4), lon_u, lat_u);
                                [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(1/20, lonlat(1:4), lon_v, lat_v);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
                            data_info = ncinfo(filename, 'u'); 
                            u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) data_info.Size(3) 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                            v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) data_info.Size(3) 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)

                            if (exist('mean_u' , 'var') ~= 1)
                                mean_u=zeros(size(u));
                                mean_v=zeros(size(v));
                            end
                            mean_u=mean_u + (u / (length(inputyear1) * length(inputmonth)));
                            mean_v=mean_v + (v / (length(inputyear1) * length(inputmonth)));
                        end
                    end
                    u_rho = u2rho_2d(mean_u')';
                    v_rho = v2rho_2d(mean_v')';
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
                    save(matname, 'u_rho','v_rho', 'cut_lon_rho', 'cut_lat_rho');
                else
                    load(matname);
                end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                u_rho=u_rho(1:size(mask_model,1),:).*mask_model;
                v_rho=v_rho(:,1:size(mask_model,2)).*mask_model;  

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
                titlename = strcat('UV mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
%                 system(['magick ', pngname, ' -trim ', pngname]);
                
                close all;
                clear lon_rho mean_u ref_vec_x_range
            end
            fig_flag=0;
        end

% start-------------------- later decadal current plot
        fig_flag=fig_flags{2,2};
        while (fig_flag)
            pngname=strcat(outfile, '_', testname,'_', regionname, '_clim_uv_',num2str(min(inputyear2),'%04i'), ...
                '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
            variable='UV';        
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
                if (exist(matname , 'file') ~= 2 || fig_flag==2)
                    for yearij=1:length(inputyear2)
                        tempyear=inputyear2(yearij);
                        yearstr=num2str(tempyear, '%04i');
                        for monthij=1:length(inputmonth)
                            tempmonth=inputmonth(monthij);
                            monthstr=num2str(tempmonth, '%02i');
                            filename=[filedir, yearstr, filesep,'pck_', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                            if (exist('lon_rho' , 'var') ~= 1)
                                lon_rho=ncread(filename, 'lon_rho');
                                lat_rho=ncread(filename, 'lat_rho');
                                lon_u=ncread(filename, 'lon_u');
                                lat_u=ncread(filename, 'lat_u');
                                lon_v=ncread(filename, 'lon_v');
                                lat_v=ncread(filename, 'lat_v');
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                                [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(1/20, lonlat(1:4), lon_u, lat_u);
                                [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(1/20, lonlat(1:4), lon_v, lat_v);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
                            data_info = ncinfo(filename, 'u'); 
                            u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) data_info.Size(3) 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                            v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) data_info.Size(3) 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)

                            if (exist('mean_u' , 'var') ~= 1)
                                mean_u=zeros(size(u));
                                mean_v=zeros(size(v));
                            end
                            mean_u=mean_u + (u / (length(inputyear2) * length(inputmonth)));
                            mean_v=mean_v + (v / (length(inputyear2) * length(inputmonth)));
                        end
                    end
                    u_rho = u2rho_2d(mean_u')';
                    v_rho = v2rho_2d(mean_v')';
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
                    save(matname, 'u_rho','v_rho', 'cut_lon_rho', 'cut_lat_rho');
                else
                    load(matname);
                end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                u_rho=u_rho(1:size(mask_model,1),:).*mask_model;
                v_rho=v_rho(:,1:size(mask_model,2)).*mask_model;  
                
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
                titlename = strcat('UV mean, ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                RemoveWhiteSpace([], 'file', pngname);
                close all;
                clear lon_rho mean_u ref_vec_x_range
            end
            fig_flag=0;
        end

% start-------------------- earlier decadal SST, SSS plot
        fig_flag=fig_flags{3,2};
        while (fig_flag)
            for varind2=1:length(all_var2)
                variable=all_var2{varind2};
                pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', variable,'_',num2str(min(inputyear1),'%04i'), ...
                    '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
%                 if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                    run(param_script);
                    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
                    if (exist(matname , 'file') ~= 2 || fig_flag==2)
                        for yearij=1:length(inputyear1)
                            tempyear=inputyear1(yearij);
                            yearstr=num2str(tempyear, '%04i');
                            for monthij=1:length(inputmonth)
                                tempmonth=inputmonth(monthij);
                                monthstr=num2str(tempmonth, '%02i');
                                filename=[filedir, yearstr, filesep,'pck_', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

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
                                mean_data=mean_data + (data / (length(inputyear1) * length(inputmonth)));
                            end
                        end
                        save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
                    else
                        load(matname);
                    end
                    
%                     [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho)
                    
                    if(strcmp(variable, 'SST'))
                        mean_data=mean_data;
                    elseif (strcmp(variable,'SSH')==1)
                        ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor(testname);
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
                        [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [2, 4, 6, 8, 10], 'color','k', ...
                                'linewidth', 1.5, 'linestyle', '-');
%                         [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
%                                 'linewidth', 1.5, 'linestyle', '-');
                            clabel(C,h2,'FontSize',13,'Color','k', ...
                                'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                    end
                    
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat(variable, ' mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
%                     if (strcmp(variable,'SSH')==1)
%                         colormap(flip(cool));
%                     else
                        colormap(jet);
%                     end

                    set(h,'fontsize',colorbar_fontsize);
                    title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
                    caxis(colorbar_lev);
                    
                    disp(['M = ', num2str(mean(mean_data(:),'omitnan'))]);
                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
                    close all;
                    clear lon_rho mean_data
%                 end
            end
            fig_flag=0;
        end

% start-------------------- later decadal SST, SSS plot
        fig_flag=fig_flags{4,2};
        while (fig_flag)
            for varind2=1:length(all_var2)
                variable=all_var2{varind2};
                pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', variable,'_',num2str(min(inputyear2),'%04i'), ...
                    '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                    run(param_script);
                    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
                    if (exist(matname , 'file') ~= 2 || fig_flag==2)
                        for yearij=1:length(inputyear2)
                            tempyear=inputyear2(yearij);
                            yearstr=num2str(tempyear, '%04i');
                            for monthij=1:length(inputmonth)
                                tempmonth=inputmonth(monthij);
                                monthstr=num2str(tempmonth, '%02i');
                                filename=[filedir, yearstr, filesep,'pck_', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

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
                                mean_data=mean_data + (data / (length(inputyear2) * length(inputmonth)));
                            end
                        end
                        save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
                    else
                        load(matname);
                    end
                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    mean_data=mean_data.*mask_model;
                    
                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;

                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                    shading(gca,m_pcolor_shading_method);   

                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat(variable, ' mean, ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
                    caxis(colorbar_lev);
                    
                    disp(['M = ', num2str(mean(mean_data(:),'omitnan'))]);
                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    RemoveWhiteSpace([], 'file', pngname);
                    close all;
                    clear lon_rho mean_data
                end
            end
            fig_flag=0;
        end

% start-------------------- earlier YSBCW plot
        fig_flag=fig_flags{5,2};
        while (fig_flag)
            if (strcmp(regionname, 'YS')==1)
                pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', 'YSBCW','_',num2str(min(inputyear1),'%04i'), ...
                    '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                    variable='BT';
                    run(param_script);
                    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
                    if (exist(matname , 'file') ~= 2 || fig_flag==2)
                        for yearij=1:length(inputyear1)
                            tempyear=inputyear1(yearij);
                            yearstr=num2str(tempyear, '%04i');
                            tempmonth=8;
                            monthstr=num2str(tempmonth, '%02i');
                            filename=[filedir, yearstr, filesep,'pck_', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

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
                            mean_data=mean_data + (data / length(inputyear1));
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
                    titlename = strcat('YSBCW-Aug', ' mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
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
                    RemoveWhiteSpace([], 'file', pngname);
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
                pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', 'YSBCW','_',num2str(min(inputyear2),'%04i'), ...
                    '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2 || fig_flag==2)       
                    var='BT';
                    run(param_script);
                    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
                    if (exist(matname , 'file') ~= 2|| fig_flag==2)
                        for yearij=1:length(inputyear2)
                            tempyear=inputyear2(yearij);
                            yearstr=num2str(tempyear, '%04i');
                            tempmonth=8;
                            monthstr=num2str(tempmonth, '%02i');
                            filename=[filedir, yearstr, filesep,'pck_', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

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
                    saveas(gcf,pngname,'tif');
                    RemoveWhiteSpace([], 'file', pngname);
                    close all;
                    clear lon_rho mean_data
                end
            end
            fig_flag=0;
        end

% % start-------------------- earlier decadal current + SST plot
%         pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_uv_SST_',num2str(min(inputyear1),'%04i'), ...
%             '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
% %         if (exist(pngname , 'file') ~= 2)        
%             for yearij=1:length(inputyear1)
%                 tempyear=inputyear1(yearij);
%                 yearstr=num2str(tempyear, '%04i');
%                 for monthij=1:length(inputmonth)
%                     tempmonth=inputmonth(monthij);
%                     monthstr=num2str(tempmonth, '%02i');
%                     filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];
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
%                     mean_u=mean_u + (u / (length(inputyear1) * length(inputmonth)));
%                     mean_v=mean_v + (v / (length(inputyear1) * length(inputmonth)));
%                 end
%             end
%             u_rho = u2rho_2d(mean_u')';
%             v_rho = v2rho_2d(mean_v')';
%             if (exist('ref_vec_x_range' , 'var') ~= 1)
%                 ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
%                 ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
%                 ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
%                 ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
%             end
%             u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
%             v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
%             
%             var='SST';
%             run(param_script);
%             for yearij=1:length(inputyear1)
%                 tempyear=inputyear1(yearij);
%                 yearstr=num2str(tempyear, '%04i');
%                 for monthij=1:length(inputmonth)
%                     tempmonth=inputmonth(monthij);
%                     monthstr=num2str(tempmonth, '%02i');
%                     filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];
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
%                     mean_data=mean_data + (data / (length(inputyear1) * length(inputmonth)));
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
% %             titlename = strcat('UV mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
% %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
%                       
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%             saveas(gcf,pngname,'tif');
%             hold off
%             close all;
%             clear lon_rho mean_u ref_vec_x_range
% %         end
% % end-------------------- earlier decadal current + SST plot
    
    end
end


