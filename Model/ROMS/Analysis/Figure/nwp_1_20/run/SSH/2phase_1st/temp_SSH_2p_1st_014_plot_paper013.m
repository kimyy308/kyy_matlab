close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','AKP2'}
all_region2 ={'NWP', 'AKP4'}

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

%         shadlev = [0 35];
%         rms_shadlev = [0 4];
%     %     trendlev = [-3 3];  %% trend lev
%         trendlev = [-10 10];  %% trend lev
%         abstrendlev =[4 7];
%         reltrendlev =[-5 5];
%         conlev  = 0:5:35;
%         meanplotlev =[-0.3 0.3];
%         trendplotlev = [3 7];
%         sshlev =[-0.7 1.3];
%         sshdifflev = [40 70];

    % for snu_desktopd
    testname='cmems';    % % need to change
    inputyear1 = [1993:2014]; % % put year which you want to plot [year year ...]
%     inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

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
        case('AKP3') %% Around Korea Peninsula
            refpolygon=akp3polygon;
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

        
%         filename = ['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

    valnum=0;
    run('C:\Users\user\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
    wrmap = bwrmap(51:100,:);

    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\MEPL\project\SSH\6th_year\figure\CMEMS\',regionname,'\'); % % where figure files will be saved
        param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_cmems_', regionname, '.m'];
%         filedir = strcat('H:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
    elseif (strcmp(system_name,'GLNXA64'))
    end


    figdir=[figrawdir,'CLIM\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);

    seasons_group={'all', 'spring', 'summer', 'fall', 'winter'};
          
    for seasons_groupi=1:length(seasons_group)
        season=seasons_group{seasons_groupi};
        switch(season)
            case 'all'
                inputmonth =1:12;  
            case 'spring'
                inputmonth =[3,4,5];  
            case 'summer'
                inputmonth =[6,7,8];  
            case 'fall'
                inputmonth =[9,10,11];  
            case 'winter'
                inputmonth =[12,1,2];  
        end  
    
% start-------------------- earlier decadal current plot
%         pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_uv_',num2str(min(inputyear1),'%04i'), ...
%             '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
    load(['D:\Data\Observation\CMEMS\',regionname, ...
            'cmems_gos_',num2str(min(inputyear1),'%04i'),'_',num2str(max(inputyear1),'%04i'), ...
            '_', num2str(inputmonth(1),'%02i'), '-', num2str(inputmonth(end),'%02i'), '.mat']);
        
        run(param_script);

        
        pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_uv_',num2str(min(inputyear1),'%04i'), ...
            '_',num2str(max(inputyear1),'%04i'), ...
            '_', num2str(inputmonth(1),'%02i'), '-', num2str(inputmonth(end),'%02i'), '.tif'); %% ~_year_month.jpg
%         tmp.tifname=strcat(dirs.figdir, tmp.testname, '_clim_uv_',num2str(min(RCM_info.years),'%04i'), ...
%             '_',num2str(max(RCM_info.years),'%04i'), ...
%             '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.tif'); %% ~_year_month.jpg

%     if (exist(pngname , 'file') ~= 2)        
        
        cut_lon_rho = cmems_lon2;
        cut_lat_rho = cmems_lat2;
        u_rho=mean_u;
        v_rho=mean_v;
        if (exist('ref_vec_x_range' , 'var') ~= 1)
            ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
            ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
            ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+1;
%             ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
            ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));        
        end
        u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
        v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     

        yearstr_min=num2str(inputyear1(1));
        yearstr_max=num2str(inputyear1(end));
        save(['D:\Data\Observation\CMEMS\',regionname, '_cmems_', 'vec','_',yearstr_min, '_', yearstr_max, ...
            num2str(inputmonth(1), '%02i'), '_', num2str(inputmonth(end), '%02i'), '.mat'], 'cut_lon_rho', 'cut_lat_rho', 'u_rho','v_rho');

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
%         titlename = strcat('UV mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
        if min(inputyear1) == max(inputyear1)
            titlename = strcat('UV, ', season(1:2), ', ', testname,',(',num2str(min(inputyear1),'%04i'),') ');  %% + glacier contribution
        else
            titlename = strcat('UV, ', season(1:2), ', ', testname, ',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
        end
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
        saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
        close all;
        clear lon_rho mean_u ref_vec_x_range
    end
%     end
% end-------------------- earlier decadal current plot
end