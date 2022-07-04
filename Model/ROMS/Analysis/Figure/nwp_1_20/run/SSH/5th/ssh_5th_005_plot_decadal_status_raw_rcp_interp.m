close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'ens09', 'ens10'};
all_testname2 = {'ens08'};

% all_testname2 = {'NorESM1-M'};
% scenname ='rcp26';

% all_region2 ={'NWP', 'YS', 'AKP2'}
% all_region2 ={'YS'};
% all_region2 ={'NWP', 'AKP4'};
all_region2 ={'NWP'};

all_var2 = {'SST', 'SSH'};
% all_var2 = {'SSH'};

% all_var2 = {'SSH', 'SST', 'SSS'};
% all_var2 = {'BT'};

% all_region2 ={'NWP'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2 scenname
        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIAS
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\User\Dropbox';
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

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
%         inputyear1 = [2026:2035]; % % put year which you want to plot [year year ...]
%         inputyear2 = [2046:2055]; % % put year which you want to plot [year year ...]
        inputyear1 = [2006:2015]; % % put year which you want to plot [year year ...]
        inputyear2 = [2091:2100]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        
        switch(testname)
            case('ens09')
                scenname='rcp26';
            case('ens08')
                scenname='rcp45';
            case('ens10')
                scenname='rcp85';
        end
                
        
%         varname ='zeta'
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        for folding=1:1
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
        end
        
        % % %         flag configuration
        for folding=1:1
            fig_flags{1,1}='earlier decadal current plot';
            fig_flags{2,1}='later decadal current plot';
            fig_flags{3,1}='earlier decadal Surface variable plot';
            fig_flags{4,1}='later decadal Surface variable plot';
            fig_flags{5,1}='earlier YSBCW plot';
            fig_flags{6,1}='later YSBCW plot';
            fig_flags{7,1}='earlier decadal current + SST plot';
            fig_flags{8,1}='later decadal current + SST plot';

            for flagi=1:8
                fig_flags{flagi,2}=0;
            end

            fig_flags{1,2}=2;
            fig_flags{2,2}=2;
            fig_flags{3,2}=1;
            fig_flags{4,2}=1;
            fig_flags{5,2}=0;
            fig_flags{6,2}=0;
%             fig_flags{7,2}=1;
            fig_flags{8,2}=0;
        end
        
        valnum=0;
        run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\MEPL\project\SSH\figures\5th_year\figure\CMIP5\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_interped_', regionname, '.m'];
%             cmip5dir = strcat('G:\Data\Model\CMIP5\'); % % where data files are
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
            matdir = strcat('D:\Data\Model\CMIP5\', testname, '\', scenname, '\mean\');
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
% start-------------------- earlier decadal current plot
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            pngname=strcat(outfile, '_', testname,'_interp_',regionname, '_clim_uv_',num2str(min(inputyear1),'%04i'), ...
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

                            varname='uo';
                            filedir = strcat(cmip5dir, varname, '/', scenname, '/interp/', testname, '/'); % % where data files are
                            flag_file_in = false;
                            list = dir( [ filedir, '/', varname, '*' ]); 
                            for kk = 1 : length( list )
                                fname_in    = list(kk).name;
                                fname_split = strsplit( fname_in, {'_','.'} );
                                fyear_str   = strsplit( fname_split{end-1}, '-' );
                                fyear_start = str2num( fyear_str{1}(1:4) );
                                fyear_end   = str2num( fyear_str{1}(1:4) );
                                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                        strcmp( fname_split{2}, 'interp' ) &&         ...
                                        strcmp( fname_split{3}, scenname ) &&      ...                 
                                        strcmp( fname_split{4}, testname ) )
                                    flag_file_in = true;            break;
                                end         
                            end         
                            if( ~flag_file_in )
                                fprintf('Source File for %04i does not Exist. Continue...\n',year);   continue;
                            end
                            ufilename=[filedir, '\', fname_in];
                            lonfilename='D:\Data\Model\CMIP5\uo\rcp45\interp\IPSL-CM5A-LR\lon_interp_rcp45_IPSL-CM5A-LR.nc';
                            latfilename='D:\Data\Model\CMIP5\uo\rcp45\interp\IPSL-CM5A-LR\lat_interp_rcp45_IPSL-CM5A-LR.nc';
                            
                            tind_u=tempmonth;

                            varname='vo';
                            filedir = strcat(cmip5dir, varname, '/', scenname, '/interp/', testname, '/'); % % where data files are
                            flag_file_in = false;
                            list = dir( [ filedir, '/', varname, '*' ]); 
                            for kk = 1 : length( list )
                                fname_in    = list(kk).name;
                                fname_split = strsplit( fname_in, {'_','.'} );
                                fyear_str   = strsplit( fname_split{end-1}, '-' );
                                fyear_start = str2num( fyear_str{1}(1:4) );
                                fyear_end   = str2num( fyear_str{1}(1:4) );
                                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                        strcmp( fname_split{2}, 'interp' ) &&         ...
                                        strcmp( fname_split{3}, scenname ) &&      ...                 
                                        strcmp( fname_split{4}, testname ) )
                                    flag_file_in = true;            break;
                                end         
                            end         
                            if( ~flag_file_in )
                                fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                            end
                            vfilename=[filedir, '\', fname_in];

                            if (exist('lon_rho' , 'var') ~= 1)
                                lon_gcm=ncread(lonfilename, 'lon');
                                lat_gcm=ncread(latfilename, 'lat');
                                [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end

                            tind_v=tempmonth;

                            data_info = ncinfo(ufilename, 'uo'); 
                            u = ncread(ufilename,'uo',[lon_min(1) lat_min(1) 1 tind_u], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                            v = ncread(vfilename,'vo',[lon_min(1) lat_min(1) 1 tind_v], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)

                            if (exist('mean_u' , 'var') ~= 1)
                                mean_u=zeros(size(u));
                                mean_v=zeros(size(v));
                            end
                            mean_u=mean_u + (u / (length(inputyear1) * length(inputmonth)));
                            mean_v=mean_v + (v / (length(inputyear1) * length(inputmonth)));
                        end
                    end
        %             u_rho = u2rho_2d(mean_u')';
        %             v_rho = v2rho_2d(mean_v')';
                    u_rho = mean_u;
                    v_rho = mean_v;
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval;
        %                 ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
        %                 ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind+(m_quiver_x_interval/2));
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind+(m_quiver_y_interval/2));
                    end
                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho=u_rho.*mask_model;
                    v_rho=v_rho.*mask_model;
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
                    save(matname, 'u_rho','v_rho', 'cut_lon_rho', 'cut_lat_rho');
                else
                    load(matname);
                end

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
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_u ref_vec_x_range
            end
            fig_flag=0;
        end

% start-------------------- later decadal current plot
        fig_flag=fig_flags{2,2};
        while (fig_flag)
            pngname=strcat(outfile, '_', testname,'_interp_',regionname, '_clim_uv_',num2str(min(inputyear2),'%04i'), ...
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

                            varname='uo';
                            filedir = strcat(cmip5dir, varname, '/', scenname, '/interp/', testname, '/'); % % where data files are
                            flag_file_in = false;
                            list = dir( [ filedir, '/', varname, '*' ]); 
                            for kk = 1 : length( list )
                                fname_in    = list(kk).name;
                                fname_split = strsplit( fname_in, {'_','.'} );
                                fyear_str   = strsplit( fname_split{end-1}, '-' );
                                fyear_start = str2num( fyear_str{1}(1:4) );
                                fyear_end   = str2num( fyear_str{1}(1:4) );
                                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                        strcmp( fname_split{2}, 'interp' ) &&         ...
                                        strcmp( fname_split{3}, scenname ) &&      ...                 
                                        strcmp( fname_split{4}, testname ) )
                                    flag_file_in = true;            break;
                                end         
                            end         
                            if( ~flag_file_in )
                                fprintf('Source File for %04i does not Exist. Continue...\n',year);   continue;
                            end
                            ufilename=[filedir, '\', fname_in];
                            lonfilename='D:\Data\Model\CMIP5\uo\rcp45\interp\IPSL-CM5A-LR\lon_interp_rcp45_IPSL-CM5A-LR.nc';
                            latfilename='D:\Data\Model\CMIP5\uo\rcp45\interp\IPSL-CM5A-LR\lat_interp_rcp45_IPSL-CM5A-LR.nc';
                            tind_u=tempmonth;

                            varname='vo';
                            filedir = strcat(cmip5dir, varname, '/',scenname, '/interp/', testname, '/'); % % where data files are
                            flag_file_in = false;
                            list = dir( [ filedir, '/', varname, '*' ]); 
                            for kk = 1 : length( list )
                                fname_in    = list(kk).name;
                                fname_split = strsplit( fname_in, {'_','.'} );
                                fyear_str   = strsplit( fname_split{end-1}, '-' );
                                fyear_start = str2num( fyear_str{1}(1:4) );
                                fyear_end   = str2num( fyear_str{1}(1:4) );
                                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                        strcmp( fname_split{2}, 'interp' ) &&         ...
                                        strcmp( fname_split{3}, scenname ) &&      ...                 
                                        strcmp( fname_split{4}, testname ) )
                                    flag_file_in = true;            break;
                                end         
                            end         
                            if( ~flag_file_in )
                                fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                            end
                            vfilename=[filedir, '\', fname_in];

                            if (exist('lon_rho' , 'var') ~= 1)
                                lon_gcm=ncread(lonfilename, 'lon');
                                lat_gcm=ncread(latfilename, 'lat');
                                [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end

                            tind_v=tempmonth;

                            data_info = ncinfo(ufilename, 'uo'); 
                            u = ncread(ufilename,'uo',[lon_min(1) lat_min(1) 1 tind_u], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                            v = ncread(vfilename,'vo',[lon_min(1) lat_min(1) 1 tind_v], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)

                            if (exist('mean_u' , 'var') ~= 1)
                                mean_u=zeros(size(u));
                                mean_v=zeros(size(v));
                            end
                            mean_u=mean_u + (u / (length(inputyear2) * length(inputmonth)));
                            mean_v=mean_v + (v / (length(inputyear2) * length(inputmonth)));
                        end
                    end
        %             u_rho = u2rho_2d(mean_u')';
        %             v_rho = v2rho_2d(mean_v')';
                    u_rho = mean_u;
                    v_rho = mean_v;
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval;
        %                 ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
        %                 ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind+(m_quiver_x_interval/2));
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind+(m_quiver_y_interval/2));
                    end

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho=u_rho.*mask_model;
                    v_rho=v_rho.*mask_model;

                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
                    save(matname, 'u_rho','v_rho', 'cut_lon_rho', 'cut_lat_rho');
                else
                    load(matname);
                end
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
                pngname=strcat(outfile, '_', testname,'_interp_',regionname, '_clim_', variable,'_',num2str(min(inputyear1),'%04i'), ...
                    '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2 || fig_flag==2)        
                    run(param_script);
                    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
                    if (exist(matname , 'file') ~= 2 || fig_flag==2)
                        for yearij=1:length(inputyear1)
                            tempyear=inputyear1(yearij);
                            yearstr=num2str(tempyear, '%04i');
                            for monthij=1:length(inputmonth)
                                tempmonth=inputmonth(monthij);
                                monthstr=num2str(tempmonth, '%02i');

                                filedir = strcat(cmip5dir, varname, '\', scenname, '\interp\', testname, '\'); % % where data files are
                                flag_file_in = false;
                                list = dir( [ filedir, '\', varname, '*' ]); 
                                for kk = 1 : length( list )
                                    fname_in    = list(kk).name;
                                    fname_split = strsplit( fname_in, {'_','.'} );
                                    fyear_str   = strsplit( fname_split{end-1}, '-' );
                                    fyear_start = str2num( fyear_str{1}(1:4) );
                                    fyear_end   = str2num( fyear_str{1}(1:4) );
                                    if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                            strcmp( fname_split{2}, 'interp' ) &&         ...
                                            strcmp( fname_split{3}, scenname ) &&      ...                 
                                            strcmp( fname_split{4}, testname ) )
                                        flag_file_in = true;            break;
                                    end         
                                end         
                                if( ~flag_file_in )
                                    fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                                end
                                filename=[filedir, '\', fname_in];
                                tind=tempmonth;
                                lonfilename='D:\Data\Model\CMIP5\uo\rcp45\interp\IPSL-CM5A-LR\lon_interp_rcp45_IPSL-CM5A-LR.nc';
                                latfilename='D:\Data\Model\CMIP5\uo\rcp45\interp\IPSL-CM5A-LR\lat_interp_rcp45_IPSL-CM5A-LR.nc';
                                if (exist('lon_rho' , 'var') ~= 1)
                                    lon_gcm=ncread(lonfilename, 'lon');
                                    lat_gcm=ncread(latfilename, 'lat');
                                    [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                    cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                    cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                end
                                data_info = ncinfo(filename, varname); 

                                if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
    % %                                 zos = ncread(vfilename,'zos',[lon_min(1) lat_min(1) 1 tind_v], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
    % %                                 zosmask=zeros(size(zos));
    % %                                 zosmask(isfinite(zos))=1;
                                    data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 

                                elseif (strcmp(variable,'SSH')==1)
                                    data = ncread(filename,varname,[lon_min(1) lat_min(1) tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                                end

                                if (exist('mean_data' , 'var') ~= 1)
                                    mean_data=zeros(size(data));
                                end
                                mean_data=mean_data + (data / (length(inputyear1) * length(inputmonth)));
                            end
                        end
                        if(strcmp(variable, 'SST'))
                            mean_data=mean_data-273.15;
                        end
                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        mean_data=mean_data.*mask_model;
                        save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
                    else
                        load(matname);
                    end
                    
%                     mean_data= mean_data -0.2296 + 0.3740; %%% for AKP4 pcolor
                    
                    
                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;

                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                    shading(gca,m_pcolor_shading_method);   

                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat(variable, ' mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
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
            fig_flag=0;
        end

% start-------------------- later decadal SST, SSS plot
        fig_flag=fig_flags{4,2};
        while (fig_flag)
            for varind2=1:length(all_var2)
                variable=all_var2{varind2};
                pngname=strcat(outfile, '_', testname,'_interp_',regionname, '_clim_', variable,'_',num2str(min(inputyear2),'%04i'), ...
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

                                filedir = strcat(cmip5dir, varname, '\', scenname, '\interp\', testname, '\'); % % where data files are
                                flag_file_in = false;
                                list = dir( [ filedir, '\', varname, '*' ]); 
                                for kk = 1 : length( list )
                                    fname_in    = list(kk).name;
                                    fname_split = strsplit( fname_in, {'_','.'} );
                                    fyear_str   = strsplit( fname_split{end-1}, '-' );
                                    fyear_start = str2num( fyear_str{1}(1:4) );
                                    fyear_end   = str2num( fyear_str{1}(1:4) );
                                    if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                            strcmp( fname_split{2}, 'interp' ) &&         ...
                                            strcmp( fname_split{3}, scenname ) &&      ...                 
                                            strcmp( fname_split{4}, testname ) )
                                        flag_file_in = true;            break;
                                    end         
                                end         
                                if( ~flag_file_in )
                                    fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                                end
                                filename=[filedir, '\', fname_in];
                                tind=tempmonth;
                                lonfilename='D:\Data\Model\CMIP5\uo\rcp45\interp\IPSL-CM5A-LR\lon_interp_rcp45_IPSL-CM5A-LR.nc';
                                latfilename='D:\Data\Model\CMIP5\uo\rcp45\interp\IPSL-CM5A-LR\lat_interp_rcp45_IPSL-CM5A-LR.nc';
                            
                                if (exist('lon_rho' , 'var') ~= 1)
                                    lon_gcm=ncread(lonfilename, 'lon');
                                    lat_gcm=ncread(latfilename, 'lat');
                                    [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                    cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                    cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                end
                                data_info = ncinfo(filename, varname); 

                                if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
    %                                 zos = ncread(vfilename,'zos',[lon_min(1) lat_min(1) 1 tind_v], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
    %                                 zosmask=zeros(size(zos));
    %                                 zosmask(isfinite(zos))=1;
                                    data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 

                                elseif (strcmp(variable,'SSH')==1)
                                    data = ncread(filename,varname,[lon_min(1) lat_min(1) tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                                end

                                if (exist('mean_data' , 'var') ~= 1)
                                    mean_data=zeros(size(data));
                                end
                                mean_data=mean_data + (data / (length(inputyear2) * length(inputmonth)));
                            end
                        end
                        if(strcmp(variable, 'SST'))
                            mean_data=mean_data-273.15;
                        end
                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        mean_data=mean_data.*mask_model;
                        save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
                    else
                        load(matname);
                    end
                    
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

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
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
                pngname=strcat(outfile, '_', testname,'_interp_',regionname, '_clim_', 'YSBCW','_',num2str(min(inputyear1),'%04i'), ...
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

                            filedir = strcat(cmip5dir, varname, '\', scenname, '\interp\', testname, '\'); % % where data files are
                            flag_file_in = false;
                            list = dir( [ filedir, '\', varname, '*' ]); 
                            for kk = 1 : length( list )
                                fname_in    = list(kk).name;
                                fname_split = strsplit( fname_in, {'_','.'} );
                                fyear_str   = strsplit( fname_split{end-1}, '-' );
                                fyear_start = str2num( fyear_str{1}(1:4) );
                                fyear_end   = str2num( fyear_str{1}(1:4) );
                                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                        strcmp( fname_split{2}, 'interp' ) &&         ...
                                        strcmp( fname_split{3}, testname ) &&      ...                 
                                        strcmp( fname_split{4}, scenname ) )
                                    flag_file_in = true;            break;
                                end         
                            end         
                            if( ~flag_file_in )
                                fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                            end
                            filename=[filedir, '\', fname_in];
                            tind=tempmonth;

                            if (exist('lon_rho' , 'var') ~= 1)
                                lon_gcm=ncread(filename, 'lon');
                                lat_gcm=ncread(filename, 'lat');
                                [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
                            data_info = ncinfo(filename, varname); 

                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                            numx=size(data,1);
                            numy=size(data,2);
                            for lonind=1:numx
                                for latind=1:numy
                                    vert_temp=squeeze(data(lonind,latind,:,1));
                                    nanind=find(isnan(vert_temp));
                                    botind=min(nanind)-1;
                                    if (botind>0)
                                        bottomtemp(lonind,latind)=vert_temp(botind)-273.15;
                                    else
                                        bottomtemp(lonind,latind)=NaN;
                                    end
                                    clear nanind botind
                                end
                            end
                            if (exist('mean_data' , 'var') ~= 1)
                                mean_data=zeros(size(bottomtemp));
                            end
                            mean_data=mean_data + (bottomtemp / length(inputyear1));
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
                pngname=strcat(outfile, '_', testname,'_interp_',regionname, '_clim_', 'YSBCW','_',num2str(min(inputyear2),'%04i'), ...
                    '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2 || fig_flag==2)       
                    variable='BT';
                    run(param_script);
                    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
                    if (exist(matname , 'file') ~= 2 || fig_flag==2)
                        for yearij=1:length(inputyear2)
                            tempyear=inputyear2(yearij);
                            yearstr=num2str(tempyear, '%04i');
                            tempmonth=8;
                            monthstr=num2str(tempmonth, '%02i');

                            filedir = strcat(cmip5dir, varname, '\', scenname, '\interp\', testname, '\'); % % where data files are
                            flag_file_in = false;
                            list = dir( [ filedir, '\', varname, '*' ]); 
                            for kk = 1 : length( list )
                                fname_in    = list(kk).name;
                                fname_split = strsplit( fname_in, {'_','.'} );
                                fyear_str   = strsplit( fname_split{end-1}, '-' );
                                fyear_start = str2num( fyear_str{1}(1:4) );
                                fyear_end   = str2num( fyear_str{1}(1:4) );
                                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                        strcmp( fname_split{2}, 'interp' ) &&         ...
                                        strcmp( fname_split{3}, testname ) &&      ...                 
                                        strcmp( fname_split{4}, scenname ) )
                                    flag_file_in = true;            break;
                                end         
                            end         
                            if( ~flag_file_in )
                                fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                            end
                            filename=[filedir, '\', fname_in];
                            tind=tempmonth;

                            if (exist('lon_rho' , 'var') ~= 1)
                                lon_gcm=ncread(filename, 'lon');
                                lat_gcm=ncread(filename, 'lat');
                                [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
                            data_info = ncinfo(filename, varname); 

                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                            numx=size(data,1);
                            numy=size(data,2);
                            for lonind=1:numx
                                for latind=1:numy
                                    vert_temp=squeeze(data(lonind,latind,:,1));
                                    nanind=find(isnan(vert_temp));
                                    botind=min(nanind)-1;
                                    if (botind>0)
                                        bottomtemp(lonind,latind)=vert_temp(botind)-273.15;
                                    else
                                        bottomtemp(lonind,latind)=NaN;
                                    end
                                    clear nanind botind
                                end
                            end
                            if (exist('mean_data' , 'var') ~= 1)
                                mean_data=zeros(size(bottomtemp));
                            end
                            mean_data=mean_data + (bottomtemp / length(inputyear2));
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
                    close all;
                    clear lon_rho mean_data
                end
            end
            fig_flag=0;
        end

% start-------------------- earlier decadal current + SST plot
        fig_flag=fig_flags{7,2};
        while (fig_flag)
            pngname=strcat(outfile, '_', testname,'_interp_',regionname, '_clim_SST_uv_',num2str(min(inputyear1),'%04i'), ...
                '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)        
                for yearij=1:length(inputyear1)
                    tempyear=inputyear1(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');

                        varname='uo';
                        filedir = strcat(cmip5dir, varname, '/', scenname, '/interp/', testname, '/'); % % where data files are
                        flag_file_in = false;
                        list = dir( [ filedir, '/', varname, '*' ]); 
                        for kk = 1 : length( list )
                            fname_in    = list(kk).name;
                            fname_split = strsplit( fname_in, {'_','.'} );
                            fyear_str   = strsplit( fname_split{end-1}, '-' );
                            fyear_start = str2num( fyear_str{1}(1:4) );
                            fyear_end   = str2num( fyear_str{1}(1:4) );
                            if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                    strcmp( fname_split{2}, 'interp' ) &&         ...
                                    strcmp( fname_split{3}, testname ) &&      ...                 
                                    strcmp( fname_split{4}, scenname ) )
                                flag_file_in = true;            break;
                            end         
                        end         
                        if( ~flag_file_in )
                            fprintf('Source File for %04i does not Exist. Continue...\n',year);   continue;
                        end
                        ufilename=[filedir, '\', fname_in];
                        tind_u=tempmonth;

                        varname='vo';
                        filedir = strcat(cmip5dir, varname, '/', scenname, '/interp/', testname, '/'); % % where data files are
                        flag_file_in = false;
                        list = dir( [ filedir, '/', varname, '*' ]); 
                        for kk = 1 : length( list )
                            fname_in    = list(kk).name;
                            fname_split = strsplit( fname_in, {'_','.'} );
                            fyear_str   = strsplit( fname_split{end-1}, '-' );
                            fyear_start = str2num( fyear_str{1}(1:4) );
                            fyear_end   = str2num( fyear_str{1}(1:4) );
                            if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                    strcmp( fname_split{2}, 'interp' ) &&         ...
                                    strcmp( fname_split{3}, testname ) &&      ...                 
                                    strcmp( fname_split{4}, scenname ) )
                                flag_file_in = true;            break;
                            end         
                        end         
                        if( ~flag_file_in )
                            fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                        end
                        vfilename=[filedir, '\', fname_in];

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_gcm=ncread(ufilename, 'lon');
                            lat_gcm=ncread(ufilename, 'lat');
                            [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end

                        tind_v=tempmonth;

                        data_info = ncinfo(ufilename, 'uo'); 
                        u = ncread(ufilename,'uo',[lon_min(1) lat_min(1) 1 tind_u], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        v = ncread(vfilename,'vo',[lon_min(1) lat_min(1) 1 tind_v], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)

                        if (exist('mean_u' , 'var') ~= 1)
                            mean_u=zeros(size(u));
                            mean_v=zeros(size(v));
                        end
                        mean_u=mean_u + (u / (length(inputyear1) * length(inputmonth)));
                        mean_v=mean_v + (v / (length(inputyear1) * length(inputmonth)));
                    end
                end
    %             u_rho = u2rho_2d(mean_u')';
    %             v_rho = v2rho_2d(mean_v')';
                u_rho = mean_u;
                v_rho = mean_v;
                if (exist('ref_vec_x_range' , 'var') ~= 1)
                    ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                    ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval;
    %                 ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
    %                 ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2));
                    ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));
                end
                u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     

                for yearij=1:length(inputyear1)
                    tempyear=inputyear1(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');

                        filedir = strcat(cmip5dir, 'thetao', '\', scenname, '\interp\', testname, '\'); % % where data files are
                        flag_file_in = false;
                        list = dir( [ filedir, '\', 'thetao', '*' ]); 
                        for kk = 1 : length( list )
                            fname_in    = list(kk).name;
                            fname_split = strsplit( fname_in, {'_','.'} );
                            fyear_str   = strsplit( fname_split{end-1}, '-' );
                            fyear_start = str2num( fyear_str{1}(1:4) );
                            fyear_end   = str2num( fyear_str{1}(1:4) );
                            if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                    strcmp( fname_split{2}, 'interp' ) &&         ...
                                    strcmp( fname_split{3}, testname ) &&      ...                 
                                    strcmp( fname_split{4}, scenname ) )
                                flag_file_in = true;            break;
                            end         
                        end         
                        if( ~flag_file_in )
                            fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                        end
                        filename=[filedir, '\', fname_in];
                        tind=tempmonth;

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_gcm=ncread(filename, 'lon');
                            lat_gcm=ncread(filename, 'lat');
                            [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, 'thetao'); 

                        sst = ncread(filename,'thetao',[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 

                        if (exist('mean_sst' , 'var') ~= 1)
                            mean_sst=zeros(size(sst));
                        end
                        mean_sst=mean_sst + (sst / (length(inputyear1) * length(inputmonth)));
                    end
                end
                mean_sst=mean_sst-273.15;

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                mean_sst=mean_sst.*mask_model;
                u_rho=u_rho.*mask_model;
                v_rho=v_rho.*mask_model;

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_sst');
                shading(gca,m_pcolor_shading_method);   

    %             [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_sst', [6 8 10], m_contour_color, 'linewidth', m_contour_linewidth);
    %             clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
    %                 'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    %             titlename = strcat('YSBCW-Aug', ' mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
    %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                colorbar_title = '(^oC)';
                title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
                colorbar_lev = [0 30]
                caxis(colorbar_lev);

                uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    %             titlename = strcat('UV mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
    %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_u ref_vec_x_range mean_sst
            end
            fig_flag=0;
        end

    end
end