close all; clear all;  clc;
warning off;

% if(isempty(gcp('nocreate')))
%     parpool(4);
% end

% % % % % % % % Ensemble of 3 members
% all_testname2 = {'test2102', 'test2103', 'test2104'};
all_testname2 = {'test2107', 'test2108', 'test2109'};
RCM_info.abbs = {'RCM-CNRM', 'RCM-EC-Veg', 'RCM-ACC'};

% RCM_info.abbs = {'RCM-CNRM', 'RCM-EC-Veg', 'RCM-ACC', 'RCM-CNRM-HR', 'RCM-CMCC'};

% scenname='historical';
scenname='ssp585';

% all_region2 ={'NWP', 'AKP4'};
all_region2 ={'AKP4'};

% all_region2 ={'YS'};

all_var2 = {'SST', 'SSH', 'SSS'};
% all_var2 = {'SSS'};
% all_var2 = {'SST'};
% all_var2 = {'SST', 'SSH', 'SSS'};
% all_var2 = {'BT'};

% all_region2 ={'NWP'}
for testnameind2=1:length(all_testname2)

for regionind2=1:length(all_region2)
    close all;
    clearvars '*' -except regionind2 testnameind2 RCM_info all_region2 all_testname2 scenname all_var2
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
%     inputyear1 = [1993:2014]; % % put year which you want to plot [year year ...]
%     inputyear1 = [1985:2014]; % % put year which you want to plot [year year ...]
    inputyear1 = [2015]; % % put year which you want to plot [year year ...]
    inputyear2 = [2050]; % % put year which you want to plot [year year ...]

%     inputmonth = [1:12]; % % put month which you want to plot [month month ...]
    inputmonth = [8]; % % put month which you want to plot [month month ...]

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

    
    figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\ENS\',scenname, filesep, regionname,filesep); % % where figure files will be saved
    param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
    dirs.figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\'); % % where figure files will be saved
    dirs.enssavedir =strcat('D:\Data\Model\ROMS\nwp_1_20\2phase_1st\RCM_ENS3\mean\');

    tmp.fs=filesep;
    tmp.regionname=regionname;
%     tmp.testname=testname;
    RCM_info.years=inputyear1:inputyear2;
    RCM_info.months=inputmonth;
    RCM_info.name=all_testname2;
    tmp.testname=RCM_info.name{testnameind2};   % % need to change
    tmp.abb=RCM_info.abbs{testnameind2};    % % need to change
    
    dirs.matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\mean\');

    run(param_script);
    
    figdir=[figrawdir,'CLIM\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);

% start-------------------- decadal current plot

    pngname=strcat(outfile, '_diff_', tmp.testname, '_',regionname, '_clim_uv_',num2str(min(inputyear1),'%04i'), ...
        '_',num2str(max(inputyear1),'%04i'), ...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'),'.tif'); %% ~_year_month.jpg
    variable='UV';
    dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'vec', tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 
    tmp.tifname=strcat(dirs.figdir, 'diff_RCM_', tmp.testname, '_clim_uv_',num2str(min(RCM_info.years),'%04i'), ...
    '_',num2str(max(RCM_info.years),'%04i'), ...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'),'.tif'); %% ~_year_month.jpg
    

                
    early_enssavefilename = [dirs.matdir, tmp.testname, '_', regionname, '_', variable, ...
        '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'),...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'), '.mat'];
    load(early_enssavefilename)
    early_u_rho=u_rho;
    early_v_rho=v_rho;
    
    late_enssavefilename = [dirs.matdir, tmp.testname, '_', regionname, '_', variable, ...
        '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), ...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'), '.mat'];
    load(late_enssavefilename)
    late_u_rho=u_rho;
    late_v_rho=v_rho;
    
    u_rho = late_u_rho - early_u_rho;
    v_rho = late_v_rho - early_v_rho;
    
%     RCM_grid.mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,RCM_grid.domain(:,1),RCM_grid.domain(:,2)));
%     RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
%     u_rho=u_rho(1:size(RCM_grid.mask_model,1),:).*RCM_grid.mask_model;
%     v_rho=v_rho(:,1:size(RCM_grid.mask_model,2)).*RCM_grid.mask_model;  
%     
    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
    mask_model(mask_model==0)=NaN;
    u_rho=u_rho(1:size(mask_model,1),:).*mask_model;
    v_rho=v_rho(:,1:size(mask_model,2)).*mask_model;  
    
    if (exist('ref_vec_x_range' , 'var') ~= 1)
        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
    end
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
    titlename = strcat('UV diff, ',tmp.testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
    saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
%                 system(['magick ', pngname, ' -trim ', pngname]);
    close all;
    clear lon_rho ens_u_rho ref_vec_x_range
    
    
% start-------------------- decadal SST, SSH, SSS plot
    for varind2=1:length(all_var2)
        variable=all_var2{varind2};
        pngname=strcat(outfile, '_diff_', tmp.testname, '_',regionname, '_clim_',variable, '_', num2str(min(inputyear1),'%04i'), ...
            '_',num2str(max(inputyear1),'%04i'), ...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'),'.tif'); %% ~_year_month.jpg
        dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, variable, tmp.fs, ...
            num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
        if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.figdir));
        end 
        tmp.tifname=strcat(dirs.figdir, 'diff_RCM_', tmp.testname, '_clim_', variable, '_',num2str(min(RCM_info.years),'%04i'), ...
                     '_',num2str(max(RCM_info.years),'%04i'), ...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'),'.tif'); %% ~_year_month.jpg
        run(param_script);
        
        early_enssavefilename = [dirs.matdir, tmp.testname, '_', regionname, '_', variable, ...
            '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'),...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'), '.mat'];
        load(early_enssavefilename)
        early_data=mean_data;

        late_enssavefilename = [dirs.matdir, tmp.testname, '_', regionname, '_', variable, ...
            '_mean_', num2str(min(inputyear2),'%04i'), '-', num2str(max(inputyear2),'%04i'), ...
                    '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'),'.mat'];
        load(late_enssavefilename)
        late_data=mean_data;

        mean_data = late_data - early_data;
       
        
%         if (strcmp(variable,'SSH')==1)
%             ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor('RCM_ENS3_historical');
%             mean_data=mean_data-ssh_correction_for_fig;
%         end

        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mean_data=mean_data.*mask_model;

        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        hold on;

        [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho)
        m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
        shading(gca,m_pcolor_shading_method);   
        if(strcmp(variable, 'SST'))
            [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [1, 1.5, 2], 'color','k', ...
                    'linewidth', 1.5, 'linestyle', '-');
                clabel(C,h2,'FontSize',13,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        end

        m_gshhs_i('color',m_gshhs_line_color)  
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
        titlename = strcat(variable, ' diff, ',tmp.testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
        colormap(jet);

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

        disp(['M = ', num2str(mean(mean_data(:),'omitnan'))]);
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
        saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        
        clear lon_rho mean_data ens_data
    end
    
    
% % % % start-------------------- decadal YSBCW plot
% % %     if (strcmp(regionname, 'YS')==1)
% % %         pngname=strcat(outfile, '_ENS_',regionname, '_clim_YSBCW', '_', num2str(min(inputyear1),'%04i'), ...
% % %             '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
% % %         variable='BT';
% % %         run(param_script);
% % %         for testnameind2=1:length(all_testname2)
% % %             testname=all_testname2{testnameind2};
% % %             filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\output\', testname, '\run\packed_monthly\'); % % where data files are          
% % %             matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\mean\');
% % %             griddir = strcat('E:\Data\Model\ROMS\nwp_1_20\output\', testname, '\run\packed_monthly\'); % % where data files are            
% % % 
% % %             matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'), '.mat'];
% % %             load(matname);
% % %             if (exist('ens_data' , 'var') ~= 1)
% % %                 ens_data=mean_data;
% % %             else
% % %                 ens_data=ens_data+mean_data;
% % %             end
% % %         end
% % %         ens_data=ens_data/length(all_testname2);
% % %         mean_data=ens_data;
% % %         
% % %         m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % %         hold on;
% % % 
% % %         m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
% % %         shading(gca,m_pcolor_shading_method);   
% % % 
% % %         [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [6 8 10], m_contour_color, 'linewidth', m_contour_linewidth);
% % %         clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
% % %             'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
% % % 
% % %         m_gshhs_i('color',m_gshhs_line_color)  
% % %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % 
% % %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % %         titlename = strcat('YSBCW-Aug', ' mean, ','ENS',',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
% % %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % 
% % %         % set colorbar 
% % %         h = colorbar;
% % %         colormap(jet);
% % %         set(h,'fontsize',colorbar_fontsize);
% % %         title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
% % %         caxis(colorbar_lev);
% % % 
% % %         set(gcf, 'PaperUnits', 'points');
% % %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % %         saveas(gcf,pngname,'tif');
% % %         RemoveWhiteSpace([], 'file', pngname);
% % %         close all;
% % %         clear lon_rho mean_data
% % %     end
        
        
end
end


