close all; clear all;  clc;
warning off;

all_testname1  = {'test57', 'test58', 'test59', 'test60'};
all_testname2  = {'test65', 'test66', 'test67', 'test68'};

all_region2 ={'AKP4'};

variable = 'SSH';
% scenname='rcp85';
for regionind2=1:length(all_region2)
        close all;

        inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
        yearstr_min=num2str(inputyear(1));
        yearstr_max=num2str(inputyear(end));
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        system_name=computer;
        for folding=1:1
            dropboxpath='C:\Users\user\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\SSH\function']));
            error_status=Func_0008_set_dropbox_path(system_name);
            run('nwp_polygon_point.m');
            regionname=all_region2{regionind2};
            [refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);

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
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[-0.1000,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        byrmap = customcolormap_preset('red-yellow-blue');
        yrmap = byrmap(129:256,:);
        
        tifname=strcat(figrawdir, 'figS10','.tif'); %% ~_year_month.jpg
        
        test_text={'(A) RCM-IPSL-L', '(B) RCM-IPSL-M', '(C) RCM-Nor', '(D) RCM-MPI','(E) RCM-IPSL-L', '(F) RCM-IPSL-M', '(G) RCM-Nor', '(H) RCM-MPI'};
        test_text2={'(RCP 4.5)', '(RCP 4.5)', '(RCP 4.5)', '(RCP 4.5)','(RCP 8.5)', '(RCP 8.5)', '(RCP 8.5)', '(RCP 8.5)'};

        
%      ---------- RCP 4.5
        for testnameind=1:4
            testname=all_testname1{testnameind};    
            varname='zeta';
            
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;
            
            rcmdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

%             load([rcmdir, rcmtestname,'_',regionname, '_rcm_gcm_glo_interped_slc_hist_rcp.mat']);
            load([rcmdir,testname,'_local_vuln_sl_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'])

            all_rcm_slc(:,:,testnameind) = rcm_slc;
            all_diff_sum4_interped(:,:,testnameind) = diff_sum4_interped;
            all_rcp_interped_std(:,:,testnameind) = rcp_interped_std;
            all_rcp_interped_det_std(:,:,testnameind) = rcp_interped_det_std;
            all_local_vuln_slr(:,:,testnameind) = local_vuln_slr;
            all_local_slr_ratio(:,:,testnameind) = local_slr_ratio;
            
%             all_mean_slc(testnameind,1) = mean_glo_slc;
            
                    
            load([rcmdir,testname,'_',regionname,'land_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
            land_model_glo2 = land_model_glo;
            
            filename = strcat('D:\Data\Model\ROMS\nwp_1_20\test53\run\2005\ocean_rst2.nc');
            land_model_interped=get_land_grid(filename,lonlat,lon_glo,lat_glo);
            
            mean_data=all_rcp_interped_det_std(:,:,testnameind);
            
            ysmask=Func_0005_get_mask_from_polygon(lon_glo, lat_glo, ys_khoapolygon);
            ys_std=mean_data.*ysmask;
            [mean_ys_std(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(ys_std, lon_glo, lat_glo);

            
            
%             mask_model = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
%             mask_model(mask_model==0)=NaN;
%             mean_data=mean_data .*mask_model;
            
% load([tidedir,regionname, '_', testname, '_model_interped_land','.mat']);

            amplev=([0 15]);       
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(8,4,[1+(testnameind-1)*8 2+(testnameind-1)*8 5+(testnameind-1)*8 6+(testnameind-1)*8]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(lon_glo',lat_glo', land_model_interped','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [], 'xtick',[120, 130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(lon_glo',lat_glo', mean_data','parent',ax{testnameind,2});
            colormap(ax{testnameind,2},parula);
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
%      ---------- RCP 4.5

%      ---------- RCP8.5
        for testnameind=5:8
            testname=all_testname2{testnameind-4};    
            varname='zeta';
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;
            
            rcmdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

            load([rcmdir,testname,'_local_vuln_sl_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'])

            all_rcm_slc(:,:,testnameind) = rcm_slc;
            all_diff_sum4_interped(:,:,testnameind) = diff_sum4_interped;
            all_rcp_interped_std(:,:,testnameind) = rcp_interped_std;
            all_rcp_interped_det_std(:,:,testnameind) = rcp_interped_det_std;
            all_local_vuln_slr(:,:,testnameind) = local_vuln_slr;
            all_local_slr_ratio(:,:,testnameind) = local_slr_ratio;         
                    
            load([rcmdir,testname,'_',regionname,'land_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
            land_model_glo2 = land_model_glo;
            
            filename = strcat('D:\Data\Model\ROMS\nwp_1_20\test53\run\2005\ocean_rst2.nc');
            land_model_interped=get_land_grid(filename,lonlat,lon_glo,lat_glo);
            
            mean_data=all_rcp_interped_det_std(:,:,testnameind);
            ysmask=Func_0005_get_mask_from_polygon(lon_glo, lat_glo, ys_khoapolygon);
            ys_std=mean_data.*ysmask;
            [mean_ys_std(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(ys_std, lon_glo, lat_glo);
            
            amplev=([0 15]);       
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(8,4,[3+(testnameind-5)*8 4+(testnameind-5)*8 7+(testnameind-5)*8 8+(testnameind-5)*8]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(lon_glo',lat_glo', land_model_interped','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [], 'xtick',[120, 130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(lon_glo',lat_glo', mean_data','parent',ax{testnameind,2});
            colormap(ax{testnameind,2},parula);
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
%      ---------- RCP 8.5

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

function land_model_interped = get_land_grid(filename,lonlat,lon_glo,lat_glo)
modelinfo=ncinfo(filename);
lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);

lon_west = abs(lon - (lonlat(1)-1));
[min_lon_west, lon_min]=min(lon_west);
lon_east = abs(lon - (lonlat(2)+1));
[min_lon_east, lon_max]=min(lon_east);
lat_south = abs(lat - (lonlat(3)-1));
[min_lat_south, lat_min]=min(lat_south);
lat_north = abs(lat - (lonlat(4)+1));
[min_lat_north, lat_max]=min(lat_north);

lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
mask_rho = ncread(filename,'mask_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
interped_mask_rho = griddata(double(lon), double(lat), mask_rho,double(lon_glo),double(lat_glo));   
land_model_interped=interped_mask_rho;
land_model_interped(interped_mask_rho==0)=1;
interped_mask_rho(interped_mask_rho==0)=NaN;
land_model_interped(isfinite(interped_mask_rho))=NaN;
end