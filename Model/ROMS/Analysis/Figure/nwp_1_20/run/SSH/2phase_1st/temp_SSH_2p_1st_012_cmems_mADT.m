close all; clear all;  clc;
% % horizontal SST trend plot (connected with SSH_mid_report1)
warning off;
% all_region2 ={'AKP','NWP','ES', 'SS', 'YS'}
% all_region2 ={'NWP', 'AKP4'};
all_region2 ={'AKP4'};

all_testname2 = {'test2102'};

% all_region2 ={'NWP'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2
    % % % 
    % % % Read Model SST
    % % % interp
    % % % get RMS
    % % % get BIAS
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        dropboxpath='C:\Users\user\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
        addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end

    % for snu_desktop
    testname=all_testname2{testnameind2};
    inputyear = [1993:2014]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        varname ='temp';
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
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
        case('AKP') %% Around Korea Peninsula
            refpolygon=akppolygon;
        case('AKP2') %% Around Korea Peninsula
            refpolygon=akp2polygon;
        case('AKP3') %% Around Korea Peninsula
            refpolygon=akp3polygon;
        case('AKP4') %% Around Korea Peninsula
            refpolygon=akp4polygon;
        case('EKB')
            refpolygon=ekbpolygon;    
        otherwise
            ('?')
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));
        cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

  % start-------------------- earlier decadal SST, SSS plot

%         for varind2=1:length(all_var2)
            variable='SSH';
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_cmems_', regionname, '.m'];
            valnum=0;
            run('C:\Users\user\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
            wrmap = bwrmap(51:100,:);
        
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\CMEMS\',regionname,'\'); % % where figure files will be saved
            figdir=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir) , 'dir') ~= 7)
                mkdir(strcat(figdir));
            end 
            outfile = strcat(figdir,regionname);
            pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', variable,'_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    %                 if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
            run(param_script);
            matname = [filedir, 'CMEMS_', regionname, '_', variable, '_mean_', num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), '.mat'];
                
            cut_lon_rho=ncread(cmemsfilename, 'lon_cmems');
            cut_lat_rho=ncread(cmemsfilename, 'lat_cmems');
            cmems_adt=ncread(cmemsfilename, 'cmems_adt');

            mean_data=mean(cmems_adt,3);
            [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho)
            ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor('CMEMS');
            mean_data=mean_data-ssh_correction_for_fig;
            save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');

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
            titlename = strcat(variable, ' mean, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
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
            saveas(gcf,pngname,'tif'); RemoveWhiteSpace([], 'file', pngname);
            close all;
            clear lon_rho mean_data
%         end
    end
end

% SSH_4th_mid_report6