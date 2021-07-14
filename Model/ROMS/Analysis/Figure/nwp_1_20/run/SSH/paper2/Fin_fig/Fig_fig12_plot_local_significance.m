close all; clear all;  clc;
warning off;

all_testname2  = {'test57', 'test58', 'test59', 'test60'};
all_testname3  = {'test65', 'test66', 'test67', 'test68'};
all_region2 ={'AKP4'};

variable = 'SSH';
% scenname='rcp45';
for regionind2=1:length(all_region2)
        close all;

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
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[0,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        tifname=strcat(figrawdir, 'fig12','.tif'); %% ~_year_month.jpg

%         f1=figure(1);
%         jet = customcolormap_preset('red-yellow-blue');
%         yrmap = jet(129:256,:);
        
        test_text={'(A) RCM, RCP 4.5', '(B) RCM, RCP 8.5'};
% start-------------------- validiation point dot
        testname ='test65';
        tempyear=2005;
        varname='zeta';
        filedir = strcat('G', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        tidedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        lon_cmems=ncread(interpedfilename, 'lon_cmems');
        lat_cmems=ncread(interpedfilename, 'lat_cmems');
        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
        cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
        cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
        
        testnameind = 1;        
        param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
        run(param_script);
        
        for temp_testind=1:length(all_testname2)
            temp_testname=all_testname2{temp_testind};
            temp_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
%                         temp_filename=[temp_filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
%                         temp_tcon=ncread(temp_filename, 'tcon');
%                         all_m2_amp(1:size(temp_tcon,1), 1:size(temp_tcon,2), temp_testind)=temp_tcon(:,:,tide_info.index(1),1).*100;
            temp_diffname=[temp_filedir, temp_testname, '_slr_ind_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
            clear slr_ind
            load(temp_diffname)
            all_slr_ind(:,:,temp_testind)=slr_ind;
        end
        mean_data=mean(all_slr_ind,3);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mean_data=mean_data .*mask_model;
        load([filedir,regionname, '_', testname, '_model_interped_land','.mat']);

        amplev=([0 70]);
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(2,4,[1 2 5 6]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on
        
        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},jet);
        caxis(amplev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 

            
            
            
        testnameind = 2;        
        param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
        run(param_script);
        
        for temp_testind=1:length(all_testname3)
            temp_testname=all_testname3{temp_testind};
            temp_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
%                         temp_filename=[temp_filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
%                         temp_tcon=ncread(temp_filename, 'tcon');
%                         all_m2_amp(1:size(temp_tcon,1), 1:size(temp_tcon,2), temp_testind)=temp_tcon(:,:,tide_info.index(1),1).*100;
            temp_diffname=[temp_filedir, temp_testname, '_slr_ind_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
            clear slr_ind
            load(temp_diffname)
            all_slr_ind(:,:,temp_testind)=slr_ind;
        end
        mean_data=mean(all_slr_ind,3);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mean_data=mean_data .*mask_model;
        load([filedir,regionname, '_', testname, '_model_interped_land','.mat']);

        amplev=([0 70]);
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(2,4,[3 4 7 8]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on
        
        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},jet);
        caxis(amplev);
        shading(gca,m_pcolor_shading_method);   
        hold on
       
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 

        
        h{1} = colorbar(ax{testnameind,2}, 'eastoutside','AxisLocation','out');
%         title(h{1},'()','fontsize',colorbar_title_fontsize);
        set(h{1},'fontsize',colorbar_fontsize);
        set(h{1}, 'Position', [pos_sb{2}(1)+0.382, pos_sb{2}(2)+0.015, 0.0231, pos_sb{2}(3)*2+0.09])  % right, up, width, height
        

            
            
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width*2 paper_position_height*1.2]) 
        print('-dtiff','-r500',tifname)

%         saveas(gcf,tifname,'tif');
        hold off;
        close all;
        
end