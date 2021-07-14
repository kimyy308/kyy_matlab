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
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\'); % % where figure files will be saved
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[0,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        tifname=strcat(figrawdir, 'fig11','.tif'); %% ~_year_month.jpg

%         f1=figure(1);
        byrmap = customcolormap_preset('red-yellow-blue');
        yrmap = byrmap(129:256,:);
        
        test_text={'(A) RCM, RCP 4.5', '(B) RCM, RCP 8.5'};
% start-------------------- validiation point dot
        testname ='test53';
        tempyear=2005;
        varname='zeta';
        filedir = strcat('J', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        tidedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        filename=[tidedir, num2str(tempyear,'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
        lon_rho=ncread(filename, 'lon_rho');
        lat_rho=ncread(filename, 'lat_rho');
        
        testnameind = 1;        
        param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
        run(param_script);
        
        for temp_testind=1:length(all_testname2)
            temp_testname=all_testname2{temp_testind};
            temp_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            temp_diffname=[temp_filedir, temp_testname, '_diff_4con_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
            clear con4_amp_diff
            load(temp_diffname)
            all_con4_amp_diff(:,:,temp_testind)=con4_amp_diff;
        end
        testname1='ens08';
        mean_data.(testname1)=mean(all_con4_amp_diff,3);
        mask_model = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mean_data.(testname1)=mean_data.(testname1) .*mask_model;
        load([filedir,regionname, '_', testname, '_model_land','.mat']);

        amplev=([-15 15]);
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(2,4,[1 2 5 6]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on
        
        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', mean_data.(testname1)','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},byrmap);
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
            temp_diffname=[temp_filedir, temp_testname, '_diff_4con_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
            clear con4_amp_diff
            load(temp_diffname)
            all_con4_amp_diff(:,:,temp_testind)=con4_amp_diff;
        end
        testname2='ens10';
        mean_data.(testname2)=mean(all_con4_amp_diff,3);
        mask_model = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mean_data.(testname2)=mean_data.(testname2) .*mask_model;
        load([filedir,regionname, '_', testname, '_model_land','.mat']);

        amplev=([-15 15]);
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(2,4,[3 4 7 8]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on
        
        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', mean_data.(testname2)','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},byrmap);
        caxis(amplev);
        shading(gca,m_pcolor_shading_method);   
        hold on
       
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 

        
        h{1} = colorbar(ax{testnameind,2}, 'eastoutside','AxisLocation','out');
        title(h{1},'(cm)','fontsize',colorbar_title_fontsize);
        set(h{1},'fontsize',colorbar_fontsize);
        set(h{1}, 'Position', [pos_sb{2}(1)+0.382, pos_sb{2}(2)+0.015, 0.0231, pos_sb{2}(3)*2+0.09])  % right, up, width, height
        

            
            
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width*2 paper_position_height*1.2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        
        RCP_corrcoef=corrcoef(mean_data.(testname1)(logical(isfinite(mean_data.(testname1)).*isfinite(mean_data.(testname2)))), mean_data.(testname2)(logical(isfinite(mean_data.(testname1)).*isfinite(mean_data.(testname2)))), 'Alpha',0.05);
        
        RCP_correcoef_r2=RCP_corrcoef(1,2)^2;
        RCP_ratio = mean_data.(testname2) ./ mean_data.(testname1);
        RCP_diff = mean_data.(testname2) - mean_data.(testname1);
        median(abs(mean_data.(testname1)(:)), 'omitnan')
        pcolor(lon_rho', lat_rho',RCP_ratio'); shading flat; colorbar; caxis([-2 2]); set(gca, 'fontsize', 30);
        pcolor(lon_rho', lat_rho',RCP_diff'); shading flat; colorbar; set(gca, 'fontsize', 30);
%         historical std mean value = 2.0024, 2.2464, 3.1554, 3.3565
        criteria_val= (2.0024 + 2.2464 + 3.1554 + 3.3565)/4.0;
        mask_crit1=NaN(size(mean_data.(testname1)));
        mask_crit1(mean_data.(testname1)>criteria_val)=1;
        mask_crit2=NaN(size(mean_data.(testname2)));
        mask_crit2(mean_data.(testname2)>criteria_val)=1;
        mask_crit3=mask_crit1 .* mask_crit2;
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        RCP_ratio_crit=RCP_ratio .* mask_crit3;
        m_pcolor(lon_rho',lat_rho',RCP_ratio_crit'); shading flat; colorbar;
        m_gshhs_i('patch', [0.8 0.8 0.8]);
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
            'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
            'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        set(gca, 'fontsize', 30);
        mean(RCP_ratio_crit(:), 'omitnan')
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        m_pcolor(lon_rho',lat_rho',RCP_ratio'); shading flat; colorbar;
        m_gshhs_i('patch', [0.8 0.8 0.8]);
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
            'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
            'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        set(gca, 'fontsize', 30);
        mean(RCP_ratio(:), 'omitnan')
        caxis([-2 2])
        tifname=strcat(figrawdir, 'fig11','.tif'); %% ~_year_month.jpg
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
        tifname=strcat(figrawdir, 'R3Q5C15','.tif'); %% ~_year_month.jpg
        print('-dtiff','-r500',tifname)
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        RCP_diff_crit=RCP_diff .* mask_crit3;
        m_pcolor(lon_rho',lat_rho',RCP_diff_crit'); shading flat; colorbar;
        m_gshhs_i('patch', [0.8 0.8 0.8]);
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
            'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
            'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        set(gca, 'fontsize', 30);
        mean(RCP_diff_crit(:), 'omitnan')
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        RCP_sum = mean_data.(testname1) + mean_data.(testname2);
        RCP_ratio_rel = RCP_diff./ (RCP_sum./2.0);
        RCP_ratio_rel = RCP_diff./ mean_data.(testname1);
        m_pcolor(lon_rho',lat_rho', RCP_ratio_rel'); shading flat; colorbar;
        m_gshhs_i('patch', [0.8 0.8 0.8]);
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
            'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
            'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        set(gca, 'fontsize', 30);
        mean(RCP_ratio_rel(:), 'omitnan');
        caxis([-10 10])
end