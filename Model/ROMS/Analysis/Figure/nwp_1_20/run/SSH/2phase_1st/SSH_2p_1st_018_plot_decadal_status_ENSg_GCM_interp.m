% %  Updated 09-Dec-2021 by Yong-Yub Kim,


close all; clear all;  clc;
warning off;


% % % % % % % % Ensemble of 3 members
all_testname2 = {'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
% scenname ='historical';
scenname ='ssp585';

% all_region2 ={'NWP', 'YS', 'AKP2'}
all_region2 ={'AKP4'}
% all_region2 ={'AKP4'};
% all_region2 ={'YS'};

all_var2 = {'SST', 'SSS', 'SSH'};
% all_var2 = {'SSH'};

% all_var2 = {'BT'};

% all_region2 ={'NWP'}
for regionind2=1:length(all_region2)
    close all;
%     clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2 scenname
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
    inputyear1 = [2050]; % % put year which you want to plot [year year ...]
%     inputyear1 = [1993:2014]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
%     inputmonth = [8]; % % put month which you want to plot [month month ...]

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
    
    valnum=0;
    run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
    wrmap = bwrmap(51:100,:);

    figrawdir =strcat('D:\MEPL\project\SSH\6th_year\figure\CMIP6\','ENS','\',scenname, filesep, regionname,'\'); % % where figure files will be saved
    param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP6_interped_', regionname, '.m'];
    cmip6dir = strcat('D:\Data\Model\CMIP6\'); % % where data files are
    dirs.figrawdir =strcat('D:\MEPL\project\SSH\6th_year\figure\nwp_1_20\'); % % where figure files will be saved
    dirs.enssavedir =strcat('D:\Data\Model\ROMS\nwp_1_20\2phase_1st\GCM_ENSg\mean\');
    tmp.fs=filesep;
    tmp.regionname=regionname;
%     tmp.testname=testname;
    RCM_info.years=inputyear1;
    
    run(param_script);

    figdir=[figrawdir,'CLIM\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);
% start-------------------- earlier decadal current plot
    pngname=strcat(outfile, '_ENSg_',scenname, '_interp_',regionname, '_clim_uv_',num2str(min(inputyear1),'%04i'), ...
        '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
    dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'vec', tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 
    tmp.tifname=strcat(dirs.figdir, 'GCM_ENSg', '_clim_uv_',num2str(min(RCM_info.years),'%04i'), ...
    '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg

    variable='UV';
    for testnameind2=1:length(all_testname2)
        testname=all_testname2{testnameind2}
        matdir = strcat('D:\Data\Model\CMIP6\',testname, filesep, scenname, filesep,'mean', filesep);
        if (exist(strcat(matdir) , 'dir') ~= 7)
            mkdir(strcat(matdir));
        end 
        matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
        if (exist(matname , 'file') ~= 2)
            for yearij=1:length(inputyear1)
                tempyear=inputyear1(yearij)
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');

                    varname='uo';
                    filedir = strcat(cmip6dir, varname, '/', scenname, '/interp/', testname, '/'); % % where data files are
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
                    filedir = strcat(cmip6dir, varname, '/', scenname, '/interp/', testname, '/'); % % where data files are
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
            
            u_rho = mean_u;
            v_rho = mean_v;
            save(matname, 'u_rho','v_rho', 'cut_lon_rho', 'cut_lat_rho');
            clear lon_rho mean_u ens_u_rho ens_v_rho
        else
            load(matname);
        end
        
        if (exist('ens_u_rho' , 'var') ~= 1)
            ens_u_rho=u_rho;
            ens_v_rho=v_rho;
        else
            ens_u_rho=ens_u_rho+u_rho;
            ens_v_rho=ens_v_rho+v_rho;
        end
    end
    ens_u_rho=ens_u_rho/length(all_testname2);
    ens_v_rho=ens_v_rho/length(all_testname2);
    u_rho=ens_u_rho;
    v_rho=ens_v_rho;
    
    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
    mask_model(mask_model==0)=NaN;
    u_rho=u_rho.*mask_model;
    v_rho=v_rho.*mask_model;
    ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
    ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval;
    ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2));
    ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));
    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;  
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
    titlename = strcat('UV mean, ','ENSg',',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
    saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
    enssavefilename = [dirs.enssavedir, 'ENSg_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
    save(enssavefilename, 'cut_lon_rho', 'cut_lat_rho', 'u_rho', 'v_rho');
    close all;
    clear lon_rho mean_u ens_u_rho ens_v_rho

% start-------------------- earlier decadal SST, SSS plot
    for varind2=1:length(all_var2)
        variable=all_var2{varind2};
        pngname=strcat(outfile, '_ENSg_',scenname, '_interp_',regionname, '_clim_', variable, '_',num2str(min(inputyear1),'%04i'), ...
        '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
        dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, variable, tmp.fs, ...
            num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
        if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.figdir));
        end 
        tmp.tifname=strcat(dirs.figdir, 'GCM_ENSg', '_clim_', variable, '_',num2str(min(RCM_info.years),'%04i'), ...
                     '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg
        run(param_script);
        for testnameind2=1:length(all_testname2)
            testname=all_testname2{testnameind2}
            matdir = strcat('D:\Data\Model\CMIP6\',testname, filesep, scenname, filesep,'mean', filesep);
            if (exist(strcat(matdir) , 'dir') ~= 7)
                mkdir(strcat(matdir));
            end 
            matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
            if (exist(matname , 'file') ~= 2)
                for yearij=1:length(inputyear1)
                    tempyear=inputyear1(yearij)
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');

                        filedir = strcat(cmip6dir, varname, '/', scenname, '/interp/', testname, '/'); % % where data files are
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
                        ufilename=[filedir, '\', fname_in];
                        tind_u=tempmonth;

                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_gcm=ncread(ufilename, 'lon');
                            lat_gcm=ncread(ufilename, 'lat');
                            [lat_rho, lon_rho]= meshgrid(lat_gcm, lon_gcm);                          
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end

                        data_info = ncinfo(ufilename, varname); 
%                         data = ncread(ufilename,varname,[lon_min(1) lat_min(1) 1 tind_u], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
                            data = ncread(ufilename,varname,[lon_min(1) lat_min(1) 1 tind_u], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        elseif (strcmp(variable,'SSH')==1)
                            data = ncread(ufilename,varname,[lon_min(1) lat_min(1) tind_u], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
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

            if (exist('ens_data' , 'var') ~= 1)
                ens_data=mean_data;
            else
                ens_data=ens_data+mean_data;
            end
            clear lon_rho mean_data
        end
        ens_data=ens_data/length(all_testname2);
        mean_data=ens_data;
        if (strcmp(variable,'SSH')==1)
            [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho)
            ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor('GCM_ENSg_historical');
            mean_data=mean_data-ssh_correction_for_fig;
        end
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mean_data=mean_data.*mask_model;

        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        hold on;

        m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
        shading(gca,m_pcolor_shading_method);   
        
        if(strcmp(variable, 'SST'))
            [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1.5, 'linestyle', '-');
            clabel(C,h2,'FontSize',13,'Color','k', ...
                 'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        end

                    
        m_gshhs_i('color',m_gshhs_line_color)  
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
        titlename = strcat(variable, ' mean, ','ENSg',',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
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
        saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        enssavefilename = [dirs.enssavedir, 'ENSg_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
        save(enssavefilename, 'cut_lon_rho', 'cut_lat_rho', 'mean_data');
        close all;
        clear lon_rho ens_data mean_data

    end



    % start-------------------- earlier YSBCW plot

    if (strcmp(regionname, 'YS')==1)
        pngname=strcat(outfile, '_ENS_',scenname, '_interp_',regionname, '_clim_', 'YSBCW', '_',num2str(min(inputyear1),'%04i'), ...
        '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
        variable='BT';
        run(param_script);
        for testnameind2=1:length(all_testname2)
            testname=all_testname2{testnameind2}
            matdir = strcat('D:\Data\Model\CMIP6\',testname, filesep, scenname, filesep,'mean', filesep);
            if (exist(strcat(matdir) , 'dir') ~= 7)
                mkdir(strcat(matdir));
            end 
            matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
            if (exist(matname , 'file') ~= 2)
                for yearij=1:length(inputyear1)
                    tempyear=inputyear1(yearij)
                    yearstr=num2str(tempyear, '%04i');
                    tempmonth=8;
                    monthstr=num2str(tempmonth, '%02i');

                    filedir = strcat(cmip6dir, varname, '\', scenname, '\interp\', testname, '\'); % % where data files are
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
                                bottomtemp(lonind,latind)=vert_temp(botind);
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
            if (exist('ens_data' , 'var') ~= 1)
                ens_data=mean_data;
            else
                ens_data=ens_data+mean_data;
            end
            clear lon_rho mean_data
        end
        ens_data=ens_data/length(all_testname2);
        mean_data=ens_data;

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
        titlename = strcat('YSBCW-Aug', ' mean, ','ENS',',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
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
        saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
        close all; 
        clear lon_rho ens_data mean_data
    end
end

