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
        er_status=Func_0008_set_dropbox_path(system_name);
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        [refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        run(param_script);
        
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\'); % % where figure files will be saved
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[-0.200,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        tifname=strcat(figrawdir, 'R3Q5C13','.tif'); %% ~_year_month.jpg

%         f1=figure(1);
%         byrmap = customcolormap_preset('red-yellow-blue');
%         yrmap = byrmap(129:256,:);
        
        test_text={'(A) GCM, RCP 4.5', '(B) GCM, RCP 8.5', '(C) RCM, RCP 4.5', '(D) RCM, RCP 8.5'};
        
        
        
        
        all_testname2  = {'test57', 'test58', 'test59', 'test60'};
        for testnameind=1:4
            rcmtestname = all_testname2{testnameind};  
            rcmdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', rcmtestname, '\run\'); % % where data files are
            load([rcmdir, rcmtestname,'_',regionname, '_rcm_gcm_glo_interped_slc_hist_rcp.mat']);
            all_gcm_slc(:,:,testnameind) = glo_slc.*100;
            all_gcm_steric_slc(:,:,testnameind) = glo_steric_slc.*100;
            all_gcm_mano_slc(:,:,testnameind) = glo_mano_slc.*100;
            all_gcm_mean_slc(testnameind,1) = mean_glo_slc.*100;
            all_gcm_mean_steric_slc(testnameind,1) = mean_glo_steric_slc.*100;
            all_gcm_mean_mano_slc(testnameind,1) = mean_glo_mano_slc.*100;
            load([rcmdir,rcmtestname,'_',regionname,'land_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
            land_model_glo2 = land_model_glo;
            
            all_rcm_slc(:,:,testnameind) = rcm_slc.*100;
            all_rcm_steric_slc(:,:,testnameind) = rcm_steric_slc.*100;
            all_rcm_mano_slc(:,:,testnameind) = rcm_mano_slc.*100;
            all_rcm_mean_slc(testnameind,1) = mean_rcm_slc.*100;
            all_rcm_mean_steric_slc(testnameind,1) = mean_rcm_steric_slc.*100;
            all_rcm_mean_mano_slc(testnameind,1) = mean_rcm_mano_slc.*100;
                
        end
        for loni=1:size(all_gcm_slc,1)
            for lati=1:size(all_gcm_slc,2)
                ens_spread(loni,lati,1)=std(all_gcm_slc(loni,lati,:));
                ens_spread(loni,lati,3)=std(all_rcm_slc(loni,lati,:));
            end
        end
        

        all_testname2  = {'test65', 'test66', 'test67', 'test68'};
        for testnameind=1:4
            rcmtestname = all_testname2{testnameind};  
            rcmdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', rcmtestname, '\run\'); % % where data files are
            load([rcmdir, rcmtestname,'_',regionname, '_rcm_gcm_glo_interped_slc_hist_rcp.mat']);
            all_gcm_slc(:,:,testnameind) = glo_slc.*100;
            all_gcm_steric_slc(:,:,testnameind) = glo_steric_slc.*100;
            all_gcm_mano_slc(:,:,testnameind) = glo_mano_slc.*100;
            all_gcm_mean_slc(testnameind,1) = mean_glo_slc.*100;
            all_gcm_mean_steric_slc(testnameind,1) = mean_glo_steric_slc.*100;
            all_gcm_mean_mano_slc(testnameind,1) = mean_glo_mano_slc.*100;
            load([rcmdir,rcmtestname,'_',regionname,'land_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
            land_model_glo2 = land_model_glo;
            
            all_rcm_slc(:,:,testnameind) = rcm_slc.*100;
            all_rcm_steric_slc(:,:,testnameind) = rcm_steric_slc.*100;
            all_rcm_mano_slc(:,:,testnameind) = rcm_mano_slc.*100;
            all_rcm_mean_slc(testnameind,1) = mean_rcm_slc.*100;
            all_rcm_mean_steric_slc(testnameind,1) = mean_rcm_steric_slc.*100;
            all_rcm_mean_mano_slc(testnameind,1) = mean_rcm_mano_slc.*100;
                
        end
        for loni=1:size(all_gcm_slc,1)
            for lati=1:size(all_gcm_slc,2)
                ens_spread(loni,lati,2)=std(all_gcm_slc(loni,lati,:));
                ens_spread(loni,lati,4)=std(all_rcm_slc(loni,lati,:));
            end
        end
%         pcolor(ens_spread(:,:,2)'); shading flat; colorbar;
        
% start-------------------- validiation point dot
        
    for testnameind=1:4
        mean_data=ens_spread(:,:,testnameind);
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        
        if testnameind<3
            sb{testnameind}=subplot(4,4,[1+2*(testnameind-1) 2+2*(testnameind-1) 5+2*(testnameind-1) 6+2*(testnameind-1)]);  % Auto-fitted to the figure.
        else
            sb{testnameind}=subplot(4,4,[9+2*(testnameind-3) 10+2*(testnameind-3) 13+2*(testnameind-3) 14+2*(testnameind-3)]);  % Auto-fitted to the figure.
        end

        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        if testnameind==2 || testnameind ==4
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
        else
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig ;
        end
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(lon_glo',lat_glo', land_model_glo2','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        if testnameind==1
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [], 'xtick',[120, 130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        elseif testnameind ==2 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [], 'xtick',[120, 130, 140], ...
                'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1}); 
        elseif testnameind ==3
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});
        elseif testnameind ==4
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1}); 
        end
        hold on
        
        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(lon_glo',lat_glo', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},parula);
        caxis([0 20]);
        shading(gca,m_pcolor_shading_method);   
        hold on
        if testnameind==1
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [], 'xtick',[120,  130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        elseif testnameind ==2 || testnameind==4
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        elseif testnameind ==3 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2}); 
        elseif testnameind ==4
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        end
        txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize-2); 
    end   
        
        h{1} = colorbar(ax{testnameind,2}, 'eastoutside','AxisLocation','out');
        title(h{1},'(cm)','fontsize',colorbar_title_fontsize);
        set(h{1},'fontsize',colorbar_fontsize);
        set(h{1}, 'Position', [pos_sb{4}(1)+0.292, pos_sb{4}(2)+0.015, 0.0231, pos_sb{2}(3)*2+0.09])  % right, up, width, height
        

            
            
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width*2 paper_position_height*1.2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        
end