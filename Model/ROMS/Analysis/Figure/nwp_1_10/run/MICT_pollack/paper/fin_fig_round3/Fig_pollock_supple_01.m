close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test53', 'test54','test55','ens06','test56'};
all_testname2 = {'test06'};

% all_region2 ={'NWP', 'YS', 'AKP2'}
% all_region2 ={'NWP','AKP2'}

all_region2 ={'NWP'};

% all_var2 = {'SST', 'SSS', 'SSH', 'BT'};
all_var2 = {'H'};

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
            dropboxpath='C:\Users\user\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_user']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_user']));
        end
        var =all_var2{1};
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
        testname=all_testname2{testnameind2}    % % need to change
        inputyear2 = [2017:2017]; % % put year which you want to plot [year year ...]
        inputmonth = [1]; % % put month which you want to plot [month month ...]

%         varname ='zeta'
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        switch(regionname)
            case('NWP') %% North western Pacific
                lonlat = [115, 162, 15, 52];  %% whole data area
                refpolygon(1,1)=lonlat(1);
                refpolygon(2,1)=lonlat(2);
                refpolygon(1,2)=lonlat(3);
                refpolygon(2,2)=lonlat(4);
%             case('NWP_10') %% North western Pacific
%                 lonlat = [115, 162, 15, 52];  %% whole data area
%                 refpolygon(1,1)=lonlat(1);
%                 refpolygon(2,1)=lonlat(2);
%                 refpolygon(1,2)=lonlat(3);
%                 refpolygon(2,2)=lonlat(4);
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

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
%         filename = ['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\user\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        load('C:\Users\user\Dropbox\source\matlab\Common\Figure\gmt_ocean_mod2.mat')  % % set colormap (gmt_ocean, nonwhite)

        wrmap = bwrmap(51:100,:);

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\Research\Ph_D_course\2021_Particle tracking of the walleye pollock\02_review\Review_round_02\figure_review_round_02\'); % % where figure files will be saved
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_', regionname, '.m'];
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\input\', testname, '\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        run(param_script);

        figdir=[figrawdir,'\'];
%         if (exist(strcat(figdir) , 'dir') ~= 7)
%             mkdir(strcat(figdir));
%         end 
        outfile = strcat(figdir,regionname);




% start-------------------- later decadal SST, SSS plot
        for varind2=1:length(all_var2)
            var=all_var2{varind2};
            pngname=strcat(outfile, '_', testname,'_',regionname, '_topog','.tif'); %% ~_year_month.jpg
%             if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear2)
                    tempyear=inputyear2(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        filename=[filedir, 'roms_grid_nwp_1_10_', testname, '.nc'];
     
                        if (exist('lon_rho' , 'var') ~= 1)
                            lon_rho=ncread(filename, 'lon_rho');
                            lat_rho=ncread(filename, 'lat_rho');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, 'h'); 
                        
                        data = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);  %% cut horizontal area [x,y,z] (wider than target area)
%                         if (exist('mean_data' , 'var') ~= 1)
%                             mean_data=zeros(size(data));
%                         end
%                         mean_data=mean_data + (data / (length(inputyear2) * length(inputmonth)));
                    end
                end

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
                mean_data= data .* mask_model;
                
                switch(regionname)
                    case('NWP')
                        m_pcolor(cut_lon_rho',cut_lat_rho', -mean_data');
                    case('AKP2') %% North western Pacific
%                         [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [200, 2000], 'color','k', ...
%                             'linewidth', 1.5, 'linestyle', '-');
%                         clabel(C,h2,'FontSize',13,'Color','k', ...
%                             'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                        m_pcolor(cut_lon_rho',cut_lat_rho', -mean_data');
                    otherwise
                        m_pcolor(cut_lon_rho',cut_lat_rho', -mean_data');
                end

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

%                 m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'YAxisLocation', 'right');
                switch(regionname)
                    case('NWP') %% North western Pacific
                        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type);
                    case('AKP2') %% North western Pacific
%                         m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'YAxisLocation', 'right');                        
                        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type);  
                    otherwise
                        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type);  
                end

%                 titlename = strcat(var, ' mean, ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
%                 title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
%                 % set colorbar 
                h = colorbar;
                colormap(gmt_ocean_mod2);
%                 cptcmap('GMT_ocean','mapping','direct');
%                 set(h, 'Limits', [-5000, 0])
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(m)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
%             end
        end
% end-------------------- later decadal SST plot

    end
end