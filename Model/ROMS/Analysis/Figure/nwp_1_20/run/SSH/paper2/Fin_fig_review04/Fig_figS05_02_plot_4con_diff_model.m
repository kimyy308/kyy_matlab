close all; clear all;  clc;
warning off;

all_testname1  = {'test57', 'test58', 'test59', 'test60'};
all_testname2  = {'test65', 'test66', 'test67', 'test68'};
RCM_testnames = {'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI', 'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'};
scennames={'rcp45', 'rcp85'};
all_region2 ={'AKP4'};

variable = 'SSH';
% scenname='rcp85';
for regionind2=1:length(all_region2)
        close all;
        
        addpath(genpath('C:\Users\User\Dropbox\source\matlab\function\'))
        [dropboxpath, erorr_status] = Func_0008_set_dropbox_path(computer);
        
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
            figrawdir =strcat('D:\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r4\'); % % where figure files will be saved
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[-0.1000,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        byrmap = customcolormap_preset('red-yellow-blue');
        yrmap = byrmap(129:256,:);
        tifname=strcat(figrawdir, 'figS05','.tif'); %% ~_year_month.jpg
        test_text={'(A) RCM-IPSL-L', '(B) RCM-IPSL-M', '(C) RCM-Nor', '(D) RCM-MPI','(E) RCM-IPSL-L', '(F) RCM-IPSL-M', '(G) RCM-Nor', '(H) RCM-MPI'};
        test_text2={'(RCP 4.5)', '(RCP 4.5)', '(RCP 4.5)', '(RCP 4.5)','(RCP 8.5)', '(RCP 8.5)', '(RCP 8.5)', '(RCP 8.5)'};

        
%      ---------- RCP 4.5
        for testnameind=1:4
            testname=all_testname1{testnameind};    
            temp_testname ='test53';
            tempyear=2005;
            varname='zeta';
            filedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            landdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            tidedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            filename=[filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
            lon_rho=ncread(filename, 'lon_rho');
            lat_rho=ncread(filename, 'lat_rho');
            
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            diffname=[tidedir, testname, '_diff_4con_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
            clear con4_amp_diff
            load(diffname)
            
            mean_data=con4_amp_diff;
            mask_model = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;
            mean_data=mean_data .*mask_model;
            load([landdir,regionname, '_', temp_testname, '_model_land','.mat']);

            amplev=([-15 15]);       
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
            colormap(ax{testnameind,2},byrmap);
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
            RCM_lon{testnameind}=lon_rho;
            RCM_lat{testnameind}=lat_rho;
            RCM_tidal_amp_change{testnameind}=mean_data;
            RCM_tidal_amp_mean(testnameind)=mean(mean_data(:), 'omitnan');
            ysecsmask=Func_0005_get_mask_from_polygon(lon_rho, lat_rho, ysecspolygon);
            tidal_mean_data_ysecs = mean_data.* ysecsmask;
            RCM_tidal_amp_mean_YSECS(testnameind)=mean(tidal_mean_data_ysecs(:), 'omitnan');
            esmask=Func_0005_get_mask_from_polygon(lon_rho, lat_rho, espolygon);
            tidal_mean_data_ES = mean_data.* esmask;
            RCM_tidal_amp_mean_ES(testnameind)=mean(tidal_mean_data_ES(:), 'omitnan');
            [indw, inde, inds, indn]=Func_0012_findind_Y(1/20, [121.2 121.2 32.8 32.8],lon_rho,lat_rho); % Jiangsu intertidal flat
            mean_data_Jiangsu=mean_data(indw-10:inde+10, inds-10:indn+10);
            RCM_tidal_amp_mean_Jiangsu(testnameind)=mean(mean_data_Jiangsu(:), 'omitnan');
             slrfilename=[tidedir, filesep, testname, '_AKP4_ssh_trend_', ...
                num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.nc'];
            minssh=mean(ncread(slrfilename, 'raw_ssh', [indw+10, inds+10, 1], [20, 20, 12]), 3);
            maxssh=mean(ncread(slrfilename, 'raw_ssh', [indw+10, inds+10, 1140-11], [20, 20, 12]), 3);
            slr=maxssh-minssh;
            RCM_slr_mean_Jiangsu(testnameind)=mean(slr(:), 'omitnan');
            
            [indw, inde, inds, indn]=Func_0012_findind_Y(1/20, [126.5 126.5 37.3 37.3],lon_rho,lat_rho); % Gyeonggi intertidal flat
            mean_data_GB=mean_data(indw-10:inde+10, inds-10:indn+10);  
            RCM_tidal_amp_mean_GB(testnameind)=mean(mean_data_GB(:), 'omitnan');
            minssh=mean(ncread(slrfilename, 'raw_ssh', [indw+10, inds+10, 1], [20, 20, 12]), 3);
            maxssh=mean(ncread(slrfilename, 'raw_ssh', [indw+10, inds+10, 1140-11], [20, 20, 12]), 3);
            slr=maxssh-minssh;
            RCM_slr_mean_GB(testnameind)=mean(slr(:), 'omitnan');
            
%             ncinfo('D:\Data\Model\ROMS\nwp_1_20\test57\run\test57_AKP4_ssh_trend_2006_2100.nc')
%             disp(' ')
            
        end
%      ---------- RCP 4.5

%      ---------- RCP8.5
        for testnameind=5:8
            testname=all_testname2{testnameind-4};    
            temp_testname ='test53';
            tempyear=2005;
            varname='zeta';
            filedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            landdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
            tidedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            filename=[filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
            lon_rho=ncread(filename, 'lon_rho');
            lat_rho=ncread(filename, 'lat_rho');
            
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            diffname=[tidedir, testname, '_diff_4con_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
            clear con4_amp_diff
            load(diffname)
            
            mean_data=con4_amp_diff;
            mask_model = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;
            mean_data=mean_data .*mask_model;
            load([landdir,regionname, '_', temp_testname, '_model_land','.mat']);

            amplev=([-15 15]);       
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
            colormap(ax{testnameind,2},byrmap);
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
            RCM_lon{testnameind}=lon_rho;
            RCM_lat{testnameind}=lat_rho;
            RCM_tidal_amp_change{testnameind}=mean_data;
            RCM_tidal_amp_mean(testnameind)=mean(mean_data(:), 'omitnan');
            ysecsmask=Func_0005_get_mask_from_polygon(lon_rho, lat_rho, ysecspolygon);
            tidal_mean_data_ysecs = mean_data.* ysecsmask;
            RCM_tidal_amp_mean_YSECS(testnameind)=mean(tidal_mean_data_ysecs(:), 'omitnan');
            esmask=Func_0005_get_mask_from_polygon(lon_rho, lat_rho, espolygon);
            tidal_mean_data_ES = mean_data.* esmask;
            RCM_tidal_amp_mean_ES(testnameind)=mean(tidal_mean_data_ES(:), 'omitnan');
            [indw, inde, inds, indn]=Func_0012_findind_Y(1/20, [121.2 121.2 32.8 32.8],lon_rho,lat_rho); % Jiangsu intertidal flat
            mean_data_Jiangsu=mean_data(indw-10:inde+10, inds-10:indn+10);
            RCM_tidal_amp_mean_Jiangsu(testnameind)=mean(mean_data_Jiangsu(:), 'omitnan');
             slrfilename=[tidedir, filesep, testname, '_AKP4_ssh_trend_', ...
                num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.nc'];
            minssh=mean(ncread(slrfilename, 'raw_ssh', [indw+10, inds+10, 1], [20, 20, 12]), 3);
            maxssh=mean(ncread(slrfilename, 'raw_ssh', [indw+10, inds+10, 1140-11], [20, 20, 12]), 3);
            slr=maxssh-minssh;
            RCM_slr_mean_Jiangsu(testnameind)=mean(slr(:), 'omitnan');
            
            [indw, inde, inds, indn]=Func_0012_findind_Y(1/20, [126.5 126.5 37.3 37.3],lon_rho,lat_rho); % Gyeonggi intertidal flat
            mean_data_GB=mean_data(indw-10:inde+10, inds-10:indn+10);  
            RCM_tidal_amp_mean_GB(testnameind)=mean(mean_data_GB(:), 'omitnan');
            minssh=mean(ncread(slrfilename, 'raw_ssh', [indw+10, inds+10, 1], [20, 20, 12]), 3);
            maxssh=mean(ncread(slrfilename, 'raw_ssh', [indw+10, inds+10, 1140-11], [20, 20, 12]), 3);
            slr=maxssh-minssh;
            RCM_slr_mean_GB(testnameind)=mean(slr(:), 'omitnan');
        end
%      ---------- RCP 8.5
%         lon_test=ncread(slrfilename, 'lon_rho', [indw+10, inds+10], [20, 20]);
%         lat_test=ncread(slrfilename, 'lat_rho', [indw+10, inds+10], [20, 20]);

        h = colorbar('eastoutside','AxisLocation','out');
        caxis(amplev);
        set(h,'fontsize',m_grid_fontsize+7);
        h_title=title(h,'(cm)','fontsize',m_grid_fontsize+3);
        set(h, 'Position', [pos_sb{end}(1)+0.352, pos_sb{end}(2)+0, 0.0231, pos_sb{end}(3)*2+0.09])  % right, up, width, height

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)
        
        

        hold off;
        close all;
        save ('Z:\내 드라이브\research\Ph_D_course\2020_SSH_CMIP5_decadal_variation_around_korea_using_downscaling_ocean_model\data_dryad\DataS05_Tidal_amp_changes_4con.mat', ...
            'RCM_lon', 'RCM_lat', 'RCM_tidal_amp_change', ...
            'RCM_testnames', 'scennames')
end




% figH = figure;
% axLH = gca;
% axRH = axes('color','none');
% mslplot{1}=plot(inputyear,Model.amp_M2(3,:), 'b','parent',axLH);
% mslplot{2}=plot(inputyear,Obs.amp_M2(3,:), 'k','parent',axRH);
% ylabel(axLH,'Model M2 Tidal amplitude (cm)')
% ylabel(axRH,'Obs M2 Tidal amplitude (cm)')
% ax_pos = get(axLH,'position');
% set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
% set(axRH,'color','none','yaxislocation','right','xtick', inputyear, 'position', ax_pos+[0 0.02 -0.01 -0.02]);
% %             set(axRH,'color','none','yaxislocation','right');
% 
% %             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
% set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
% set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
% set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
% xlabel(axRH, 'Year');
% 
% title([num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
% % datetick('x','yyyy','keepticks')
% axis tight;
% % ylim(meanplotlev2)
% set(mslplot{1},'LineWidth',2);
% set(mslplot{2},'LineWidth',2);
% grid on
% 
% lgd=legend([mslplot{1} mslplot{2}], 'Model','TG-UST');
% 
% %             lgd=legend('Model','TG-UST');
% set(lgd,'FontSize',15);
% set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
% set(lgd,'Orientation','horizontal');
% 
% set(gcf,'PaperPosition', [0 0 36 12]) 
% saveas(gcf,pngname,'tif');
% grid off
% close all;