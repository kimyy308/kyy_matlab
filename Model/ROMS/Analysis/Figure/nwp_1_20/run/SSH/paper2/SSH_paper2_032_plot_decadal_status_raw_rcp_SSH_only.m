close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
all_testname2 = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
% all_testname2 = {'NorESM1-M'};
all_scenname2 ={'rcp26', 'rcp85'};

% all_region2 ={'NWP', 'YS', 'AKP2'}
all_region2 ={'NWP'};

% all_region2 ={'YS'};
% all_var2 = {'SST', 'SSS'};
all_var2 = {'SSH'};

% all_region2 ={'NWP'}
for scenind2=1:length(all_scenname2)
    for testnameind2=1:length(all_testname2)
        for regionind2=1:length(all_region2)
            close all;
            clearvars '*' -except regionind2 testnameind2 scenind2 all_region2 all_testname2 all_var2 all_scenname2

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
            scenname=all_scenname2{scenind2};

            inputyear1 = [2006:2015]; % % put year which you want to plot [year year ...]
            inputyear2 = [2051:2060]; % % put year which you want to plot [year year ...]
            inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
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

                fig_flags{3,2}=2;
                fig_flags{4,2}=2;
            end

            % % % for EKB
            % regionname='EKB';
            % lonlat = [127, 129.5, 38, 40.5];

    %         load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
    %             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
    %         filename = ['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
    %                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

            valnum=0;
            run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
            wrmap = bwrmap(51:100,:);

            if (strcmp(system_name,'PCWIN64'))
                % % for windows
                figrawdir =strcat('F:\OneDrive - 서울대학교\MEPL\project\SSH\5th_year\figure\CMIP5\',testname,'\',regionname,'\'); % % where figure files will be saved
                param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m']
    %             cmip5dir = strcat('G:\Data\Model\CMIP5\'); % % where data files are
                cmip5dir = strcat('E:\Data\Model\CMIP5\'); % % where data files are
            elseif (strcmp(system_name,'GLNXA64'))
            end

            run(param_script);

            figdir=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir) , 'dir') ~= 7)
                mkdir(strcat(figdir));
            end 
            outfile = strcat(figdir,regionname);


            data_path_ar5=[cmip5dir, '\WG1AR5_Ch13SM_datafiles\13.SM.1\']; 

            %---------- glacier AR5
            fname = [data_path_ar5,'/',scenname,'_glaciermid.txt'];
            fid = fopen( fname );
            tmp = textscan( fid, '%f %f' );
            fclose( fid );
            time_ar5    = tmp{1};
            glacier_ar5 = tmp{2};

            %---------- land water AR5
            fname = [data_path_ar5,'/',scenname,'_landwatermid.txt'];
            fid = fopen( fname );
            tmp = textscan( fid, '%f %f' );
            fclose( fid );
            landwater_ar5 = tmp{2};

            %---------- ice sheets AR5
            fname = [data_path_ar5,'/',scenname,'_antnetmid.txt'];
            fid = fopen( fname );
            tmp = textscan( fid, '%f %f' );
            fclose( fid );
            icesheets_ar5 = tmp{2};
            fname = [data_path_ar5,'/',scenname,'_greennetmid.txt'];
            fid = fopen( fname );
            tmp = textscan( fid, '%f %f' );
            fclose( fid );
            icesheets_ar5 = icesheets_ar5 + tmp{2};

            %-- add 2006
            glacier_time = [time_ar5(1)-1; time_ar5; time_ar5(end)+1];
            glacier_sle = [ 2*glacier_ar5(1)-glacier_ar5(2); glacier_ar5; 2*glacier_ar5(end)-glacier_ar5(end-1)];
            landwater_time = [time_ar5(1)-1; time_ar5; time_ar5(end)+1];
            landwater_sle = [ 2*landwater_ar5(1)-landwater_ar5(2); landwater_ar5; 2*landwater_ar5(end)-landwater_ar5(end-1)];
            icesheets_time = [time_ar5(1)-1; time_ar5; time_ar5(end)+1];
            icesheets_sle = [ 2*icesheets_ar5(1)-icesheets_ar5(2); icesheets_ar5; 2*icesheets_ar5(end)-icesheets_ar5(end-1)];
            % %     calculate correction values
    %         time_yr = (0.5:11.5)'/12 + year;
    %         time_yr = (0.5:11.5)'/12 + 2006:2100;
            time_yr = 2006+1/24:1/12:2101-1/24;

            glacier   = interp1( glacier_time,   glacier_sle,   time_yr );
            landwater = interp1( landwater_time, landwater_sle, time_yr );
            icesheets = interp1( icesheets_time, icesheets_sle, time_yr );

            %------ select zostoga file
    %         flag_file_ts = false;
    %         list = dir( [ path_ts, '/', TSname, '*' ] );
    
%             filedir = strcat(cmip5dir, 'zostoga', '\', scenname, '\Omon\', testname, '\'); % % where data files are
%             flag_file_ts = false;
%             list = dir( [ filedir, '\', 'zostoga', '*' ]); 
            fullyear=2006:2100;
            zostoga_ind=1;
            for yearij=1:length(fullyear)
                tempyear=fullyear(yearij);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');
                    
                    filedir = strcat(cmip5dir, 'zostoga', '\', scenname, '\Omon\', testname, '\'); % % where data files are
                    flag_file_ts = false;
                    list = dir( [ filedir, '\', 'zostoga', '*' ]); 
                    
                    for kk = 1 : length( list )
                        fname_ts    = list(kk).name;
                        fname_ts_split = strsplit( fname_ts, {'_','.'} );
                        fyear_ts_str   = strsplit( fname_ts_split{end-1}, '-' );
                        fyear_ts_start = str2num( fyear_ts_str{1}(1:4) );
                        fyear_ts_end   = str2num( fyear_ts_str{2}(1:4) );
                        if( tempyear >= fyear_ts_start && tempyear <= fyear_ts_end &&  ...
                                strcmp( fname_ts_split{2}, 'Omon' ) &&         ...
                                strcmp( fname_ts_split{3}, testname ) &&      ...
                                strcmp( fname_ts_split{4}, scenname ) )
                            flag_file_ts = true;
                            break;
                        end
                    end
                    if( ~flag_file_ts )
                        fprintf('%s File for %04i does not Exist. Continue...\n',TSname,year);  continue;
                    end
                    filename=[filedir, '\', fname_ts];
                    tind=(tempyear-fyear_ts_start)*12+tempmonth;
                    zostoga(zostoga_ind) = ncread(filename,'zostoga',[tind], [1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');
                    
                    filedir = strcat(cmip5dir, 'zos', '\', scenname, '\Omon\', testname, '\'); % % where data files are
                    flag_file_ts = false;
                    list = dir( [ filedir, '\', 'zos', '*' ]); 
                    
                    for kk = 1 : length( list )
                        fname_ts    = list(kk).name;
                        fname_ts_split = strsplit( fname_ts, {'_','.'} );
                        fyear_ts_str   = strsplit( fname_ts_split{end-1}, '-' );
                        fyear_ts_start = str2num( fyear_ts_str{1}(1:4) );
                        fyear_ts_end   = str2num( fyear_ts_str{2}(1:4) );
                        if( tempyear >= fyear_ts_start && tempyear <= fyear_ts_end &&  ...
                                strcmp( fname_ts_split{2}, 'Omon' ) &&         ...
                                strcmp( fname_ts_split{3}, testname ) &&      ...
                                strcmp( fname_ts_split{4}, scenname ) )
                            flag_file_ts = true;
                            break;
                        end
                    end
                    if( ~flag_file_ts )
                        fprintf('%s File for %04i does not Exist. Continue...\n',TSname,year);  continue;
                    end
                    filename=[filedir, '\', fname_ts];
                    tind=(tempyear-fyear_ts_start)*12+tempmonth;
                    
                    if (exist('lon_GCM' , 'var') ~= 1)
                        lon_GCM=ncread(filename, 'lon');
                        lat_GCM=ncread(filename, 'lat');
                    end     
                    
                    zos(:,:,zostoga_ind) = ncread(filename,'zos',[1 1 tind], [inf inf 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                    area_tgt = sum(reshape(~isnan(zos(:,:,zostoga_ind)).*cosd(lat_GCM),[],1), 'omitnan');
                    zos_m(tind) = sum(reshape(zos(:,:,zostoga_ind) .*cosd(lat_GCM),[],1),'omitnan') / area_tgt;
                    
                    zostoga_ind=zostoga_ind+1;
                end
            end
            %------ get global mean value
            

            
            
    % start-------------------- earlier decadal SSH plot
            fig_flag=fig_flags{3,2};
            while (fig_flag)
                for varind2=1:length(all_var2)
                    variable=all_var2{varind2};
                    pngname=strcat(outfile, '_', scenname, '_', testname,'_',regionname, '_clim_', variable,'_',num2str(min(inputyear1),'%04i'), ...
                        '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2 || fig_flag==2)        
                        run(param_script);
                        for yearij=1:length(inputyear1)
                            tempyear=inputyear1(yearij);
                            yearstr=num2str(tempyear, '%04i');
                            for monthij=1:length(inputmonth)
                                tempmonth=inputmonth(monthij);
                                monthstr=num2str(tempmonth, '%02i');

                                filedir = strcat(cmip5dir, varname, '\', scenname, '\Omon\', testname, '\'); % % where data files are
                                flag_file_in = false;
                                list = dir( [ filedir, '\', varname, '*' ]); 
                                for kk = 1 : length( list )
                                    fname_in    = list(kk).name;
                                    fname_split = strsplit( fname_in, {'_','.'} );
                                    fyear_str   = strsplit( fname_split{end-1}, '-' );
                                    fyear_start = str2num( fyear_str{1}(1:4) );
                                    fyear_end   = str2num( fyear_str{2}(1:4) );
                                    if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                            strcmp( fname_split{2}, 'Omon' ) &&         ...
                                            strcmp( fname_split{3}, testname ) &&      ...                 
                                            strcmp( fname_split{4}, scenname ) )
                                        flag_file_in = true;            break;
                                    end         
                                end         
                                if( ~flag_file_in )
                                    fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                                end
                                filename=[filedir, '\', fname_in];
                                tind=(tempyear-fyear_start)*12+tempmonth;

                                if (exist('lon_rho' , 'var') ~= 1)
                                    lon_rho=ncread(filename, 'lon');
                                    lat_rho=ncread(filename, 'lat');
                                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                    cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                    cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                end
                                data_info = ncinfo(filename, varname); 

                                data = ncread(filename,varname,[lon_min(1) lat_min(1) tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 

                                add_ind=(tempyear-2006)*12+tempmonth;
                                data = data - zos_m(add_ind) + glacier(add_ind) + landwater(add_ind) + icesheets(add_ind) + zostoga(add_ind);

                                if (exist('mean_data' , 'var') ~= 1)
                                    mean_data=zeros(size(data));
                                end
                                mean_data=mean_data + (data / (length(inputyear1) * length(inputmonth)));
                            end
                        end
                        if(strcmp(variable, 'SST'))
                            mean_data=mean_data-273.15;
                        end

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

    % start-------------------- later decadal SSH plot
            fig_flag=fig_flags{4,2};
            while (fig_flag)
                for varind2=1:length(all_var2)
                    variable=all_var2{varind2};
                    pngname=strcat(outfile, '_', scenname, '_', testname,'_',regionname, '_clim_', variable,'_',num2str(min(inputyear2),'%04i'), ...
                        '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2 || fig_flag==2)        
                        run(param_script);
                        for yearij=1:length(inputyear2)
                            tempyear=inputyear2(yearij);
                            yearstr=num2str(tempyear, '%04i');
                            for monthij=1:length(inputmonth)
                                tempmonth=inputmonth(monthij);
                                monthstr=num2str(tempmonth, '%02i');

                                filedir = strcat(cmip5dir, varname, '\', scenname, '\Omon\', testname, '\'); % % where data files are
                                flag_file_in = false;
                                list = dir( [ filedir, '\', varname, '*' ]); 
                                for kk = 1 : length( list )
                                    fname_in    = list(kk).name;
                                    fname_split = strsplit( fname_in, {'_','.'} );
                                    fyear_str   = strsplit( fname_split{end-1}, '-' );
                                    fyear_start = str2num( fyear_str{1}(1:4) );
                                    fyear_end   = str2num( fyear_str{2}(1:4) );
                                    if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                                            strcmp( fname_split{2}, 'Omon' ) &&         ...
                                            strcmp( fname_split{3}, testname ) &&      ...                 
                                            strcmp( fname_split{4}, scenname ) )
                                        flag_file_in = true;            break;
                                    end         
                                end         
                                if( ~flag_file_in )
                                    fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                                end
                                filename=[filedir, '\', fname_in];
                                tind=(tempyear-fyear_start)*12+tempmonth;

                                if (exist('lon_rho' , 'var') ~= 1)
                                    lon_rho=ncread(filename, 'lon');
                                    lat_rho=ncread(filename, 'lat');
                                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                                    cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                    cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                end
                                data_info = ncinfo(filename, varname); 

                                data = ncread(filename,varname,[lon_min(1) lat_min(1) tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 

                                add_ind=(tempyear-2006)*12+tempmonth;
                                data = data - zos_m(add_ind) + glacier(add_ind) + landwater(add_ind) + icesheets(add_ind) + zostoga(add_ind);

                                if (exist('mean_data' , 'var') ~= 1)
                                    mean_data=zeros(size(data));
                                end
                                mean_data=mean_data + (data / (length(inputyear2) * length(inputmonth)));
                            end
                        end
                        if(strcmp(variable, 'SST'))
                            mean_data=mean_data-273.15;
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

        end
    end
end