close all; clear all;  clc;
warning off;

% all_testname1  = {'test57', 'test58', 'test59', 'test60'};
% all_testname2  = {'test65', 'test66', 'test67', 'test68'};

testnames = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR', ...
                'test61', 'test62', 'test63', 'test64'};
for i=1:8
    scennames{i} = 'rcp26';
end

inputyears = [2006, 2015; ...
              2026, 2035; ...
              2046, 2055; ...
              2091, 2100; ...
              2006, 2015; ...
              2026, 2035; ...
              2046, 2055; ...
              2091, 2100];

all_region2 ={'AKP4'};

 inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
        yearstr_min=num2str(inputyear(1));
        yearstr_max=num2str(inputyear(end));
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        system_name=computer;
        for folding=1:1
            if (strcmp(system_name,'PCWIN64'))
                % % for windows
                dropboxpath='C:\Users\User\Dropbox';
                addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
                addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
                addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
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
            
            run('nwp_polygon_point.m');
            regionname=all_region2{regionind2};
            switch(regionname)
            case('NWP')
                refpolygon=nwppolygon;
            case('AKP4') %% Around Korea Peninsula
                refpolygon=akp4polygon;
            end
            lonlat(1)=min(refpolygon(:,1));
            lonlat(2)=max(refpolygon(:,1));
            lonlat(3)=min(refpolygon(:,2));
            lonlat(4)=max(refpolygon(:,2));

            load('C:\Users\User\Dropbox\source\matlab\Common\Figure\gmt_ocean_mod2.mat')  % % set colormap (gmt_ocean, nonwhite)
        end
       
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig\'); % % where figure files will be saved
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[-0.1000,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        byrmap = customcolormap_preset('red-yellow-blue');
        yrmap = byrmap(129:256,:);
        tifname=strcat(figrawdir, 'figS06','.tif'); %% ~_year_month.jpg

%         test_text={'(A) RCM-IPSL-L', '(B) RCM-IPSL-M', '(C) RCM-Nor', '(D) RCM-MPI','(E) RCM-IPSL-L', '(F) RCM-IPSL-M', '(G) RCM-Nor', '(H) RCM-MPI'};
%         test_text2={'(RCP 4.5)', '(RCP 4.5)', '(RCP 4.5)', '(RCP 4.5)','(RCP 8.5)', '(RCP 8.5)', '(RCP 8.5)', '(RCP 8.5)'};


variable = 'SSH';
% scenname='rcp85';
for regionind2=1:length(all_region2)
        close all;

        for testnameind=1:4
            testname=all_testname1{testnameind};    
            temp_testname ='test53';
            tempyear=2005;
            varname='zeta';
            filedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            landdir = strcat('J', ':\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            tidedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            interpedfilename = strcat(tidedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            lon_cmems=ncread(interpedfilename, 'lon_cmems');
            lat_cmems=ncread(interpedfilename, 'lat_cmems');
            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
            lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            diffname=[tidedir, testname, '_slr_ind_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
            clear slr_ind
            load(diffname)
            
            mean_data=slr_ind;
            mask_model = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;
            mean_data=mean_data .*mask_model;
            load([tidedir,regionname, '_', testname, '_model_interped_land','.mat']);

            amplev=([0 70]);       
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(8,4,[1+(testnameind-1)*8 2+(testnameind-1)*8 5+(testnameind-1)*8 6+(testnameind-1)*8]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', model_land','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [], 'xtick',[120, 130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', mean_data','parent',ax{testnameind,2});
            colormap(ax{testnameind,2},jet);
            caxis(amplev);
            shading(gca,m_pcolor_shading_method);   
            hold on
            
           if testnameind==4
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            else
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [], 'xtick',[120,  130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            end
            txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 
            txt{testnameind,3}=m_text(119, 48, test_text2{testnameind}, 'FontSize', m_grid_fontsize+4); 
        end

        for testnameind=5:8
            testname=all_testname2{testnameind-4};    
            temp_testname ='test53';
            tempyear=2005;
            varname='zeta';
            filedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            landdir = strcat('J', ':\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            tidedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            interpedfilename = strcat(tidedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            lon_cmems=ncread(interpedfilename, 'lon_cmems');
            lat_cmems=ncread(interpedfilename, 'lat_cmems');
            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
            lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            diffname=[tidedir, testname, '_slr_ind_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
            clear slr_ind
            load(diffname)
            
            mean_data=slr_ind;
            mask_model = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;
            mean_data=mean_data .*mask_model;
            load([tidedir,regionname, '_', testname, '_model_interped_land','.mat']);
            
            amplev=([0 70]);       
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(8,4,[3+(testnameind-5)*8 4+(testnameind-5)*8 7+(testnameind-5)*8 8+(testnameind-5)*8]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', model_land','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [], 'xtick',[120, 130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', mean_data','parent',ax{testnameind,2});
            colormap(ax{testnameind,2},jet);
            caxis(amplev);
            shading(gca,m_pcolor_shading_method);   
            hold on
            
           if testnameind==8
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            else
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            end
            txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 
            txt{testnameind,3}=m_text(119, 48, test_text2{testnameind}, 'FontSize', m_grid_fontsize+4); 
        end

        
        h = colorbar('eastoutside','AxisLocation','out');
        caxis(amplev);
        set(h,'fontsize',m_grid_fontsize+7);
%         h_title=title(h,'(cm)','fontsize',m_grid_fontsize+3);
        set(h, 'Position', [pos_sb{end}(1)+0.352, pos_sb{end}(2)+0, 0.0231, pos_sb{end}(3)*2+0.09])  % right, up, width, height

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        
end



function [cut_lon_rho, cut_lat_rho, mean_data, land] = prepare_data_for_fig (testname, scenname, varname, inputyear);

    










end