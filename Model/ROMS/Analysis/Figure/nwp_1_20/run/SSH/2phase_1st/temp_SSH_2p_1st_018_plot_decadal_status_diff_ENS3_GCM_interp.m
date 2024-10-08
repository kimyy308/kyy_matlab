close all; clear all;  clc;
warning off;


% % % % % % % % Ensemble of 3 members
all_testname2 = {'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2'};
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
    inputyear1 = [2015]; % % put year which you want to plot [year year ...]
    inputyear2 = [2050];
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

    figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\CMIP6\','ENS','\',scenname, filesep, regionname,'\'); % % where figure files will be saved
    param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP6_interped_', regionname, '.m'];
    cmip6dir = strcat('D:\Data\Model\CMIP6\'); % % where data files are
    dirs.figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\'); % % where figure files will be saved
    dirs.enssavedir =strcat('D:\Data\Model\ROMS\nwp_1_20\2phase_1st\GCM_ENS3\mean\');
    tmp.fs=filesep;
    tmp.regionname=regionname;
%     tmp.testname=testname;
    RCM_info.years=inputyear1:inputyear2;
    
    run(param_script);

    figdir=[figrawdir,'CLIM\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);
% start-------------------- earlier decadal current plot
    pngname=strcat(outfile, '_diff_ENS3_',scenname, '_interp_',regionname, '_clim_uv_',num2str(min(inputyear1),'%04i'), ...
        '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
    dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'vec', tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 
    tmp.tifname=strcat(dirs.figdir, 'diff_GCM_ENS3', '_clim_uv_',num2str(min(RCM_info.years),'%04i'), ...
    '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg

    variable='UV';
    
    early_enssavefilename = [dirs.enssavedir, 'ENS3_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
    load(early_enssavefilename)
    early_u_rho=u_rho;
    early_v_rho=v_rho;
    
    late_enssavefilename = [dirs.enssavedir, 'ENS3_', regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
    load(late_enssavefilename)
    late_u_rho=u_rho;
    late_v_rho=v_rho;
    
    u_rho = late_u_rho - early_u_rho;
    v_rho = late_v_rho - early_v_rho;
    
%     mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%     mask_model(mask_model==0)=NaN;
%     u_rho=u_rho.*mask_model;
%     v_rho=v_rho.*mask_model;

    ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
    ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval;
    ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2));
    ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));
    u_rho(ref_vec_x_range,ref_vec_y_range)=0.5;
    v_rho(ref_vec_x_range,ref_vec_y_range)=0.00001;  
    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    hold on;
    m_gshhs_i('color',m_gshhs_line_color)  
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                    v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, '0.5 m/s', 'FontSize', m_quiver_ref_text_fontsize); 
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    titlename = strcat('UV diff, ','ENS3',',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
    saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
    close all;
    clear lon_rho mean_u ens_u_rho ens_v_rho

% start-------------------- earlier decadal SST, SSS plot
    for varind2=1:length(all_var2)
        variable=all_var2{varind2};
        pngname=strcat(outfile, '_diff_ENS3_',scenname, '_interp_',regionname, '_clim_', variable, '_',num2str(min(inputyear1),'%04i'), ...
        '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
        dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, variable, tmp.fs, ...
            num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
        if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.figdir));
        end 
        tmp.tifname=strcat(dirs.figdir, 'diff_GCM_ENS3', '_clim_', variable, '_',num2str(min(RCM_info.years),'%04i'), ...
                     '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg
        run(param_script);
        
        early_enssavefilename = [dirs.enssavedir, 'ENS3_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
        load(early_enssavefilename)
        early_data=mean_data;

        late_enssavefilename = [dirs.enssavedir, 'ENS3_', regionname, '_', variable, '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), '.mat'];
        load(late_enssavefilename)
        late_data=mean_data;

        mean_data = late_data - early_data;
%         
%         mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%         mask_model(mask_model==0)=NaN;
%         mean_data=mean_data.*mask_model;

        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        hold on;

        m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
        shading(gca,m_pcolor_shading_method);   
        
        if(strcmp(variable, 'SST'))
            [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [0, 1, 2, 3], 'color','k', ...
                    'linewidth', 1.5, 'linestyle', '-');
            clabel(C,h2,'FontSize',13,'Color','k', ...
                 'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        end

                    
        m_gshhs_i('color',m_gshhs_line_color)  
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
        titlename = strcat(variable, ' diff, ','ENS3',',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
        set(h,'fontsize',colorbar_fontsize);
        title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
        if(strcmp(variable, 'SST'))
            caxis([-1 3]);
            colormap(jet);
        elseif(strcmp(variable, 'SSS'))
            caxis([-1 1]);
            [byrmap, error_status] = Func_0009_get_colormaps('byr3', dropboxpath)
            colormap(byrmap);
        elseif(strcmp(variable, 'SSH'))   
            caxis([0 0.4] + (inputyear2-2050)*0.006);
            colormap(jet);
        end

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
        saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
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

