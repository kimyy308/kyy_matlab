% %  Updated 26-Apr-2021 by Yong-Yub Kim, 

% start-------------------- earlier decadal SST, SSS plot
for varind2=1:length(RCM_info.vars)
%% set variable name & figure directory
    tmp.variable=RCM_info.vars{varind2};
    dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 

%% set figure file name
    tmp.tifname=strcat(dirs.figdir, tmp.testname, '_surf_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
    if (exist(tmp.tifname , 'file') ~= 2 || fig_flag==2)      
        run(tmp.param_script);
%% set data file name (mat) for reading
        tmp.testname_reana = 'test06';
        tmp.matname_reana = [dirs.matdir_reana, tmp.testname_reana, '_', tmp.regionname, '_', tmp.variable,...
            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
        '_', RCM_info.season, '.mat'];
        if (exist(tmp.matname_reana , 'file') ~= 2)
            disp('please get data from subroutine_004_001 first')
        else
            load(tmp.matname_reana);
            RCM_data_reana=RCM_data;
            RCM_grid_reana=RCM_grid;
        end 

        tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
        '_', RCM_info.season, '.mat'];
        if (exist(tmp.matname , 'file') ~= 2)
            disp('please get data from subroutine_004_001 first')
        else
            load(tmp.matname);
        end 
        
%% interpolation
        RCM_data_interped.all=NaN(size(RCM_data_reana.all));
        RCM_data_interped.se=NaN(size(RCM_data_reana.all));
        RCM_data_interped.bias=NaN(size(RCM_data_reana.all));

        for yearij=1:size(RCM_data.all,3)
            for monthij = 1:size(RCM_data.all,4)
                tmp.data=squeeze(RCM_data.all(:,:,yearij,monthij));
                tmp.data_interped=griddata(RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho, tmp.data, RCM_grid_reana.cut_lon_rho, RCM_grid_reana.cut_lat_rho);
                RCM_data_interped.all(:,:,yearij,monthij)=tmp.data_interped;
%                 for lonij=1:size(RCM_data_interped.all,1)
%                     for latij=1:size(RCM_data_interped.all,2)
                        RCM_data_interped.se(:, :, yearij, monthij) =  (RCM_data_interped.all(:,:,yearij,monthij) - RCM_data_reana.all(:,:,yearij,monthij)).^2;
                        RCM_data_interped.bias(:, :, yearij, monthij) =  RCM_data_interped.all(:,:,yearij,monthij) - RCM_data_reana.all(:,:,yearij,monthij);
%                     end
%                 end
            end
        end
        
        tmp.data_interped_se=reshape(RCM_data_interped.se, [size(RCM_data_interped.all,1), size(RCM_data_interped.all,2), size(RCM_data_interped.all,3)*size(RCM_data_interped.all,4)]);
        RCM_data_interped.mse=mean(tmp.data_interped_se,3);
        RCM_data_interped.rmse=sqrt(RCM_data_interped.mse);
        tmp.data_interped_bias=reshape(RCM_data_interped.bias, [size(RCM_data_interped.all,1), size(RCM_data_interped.all,2), size(RCM_data_interped.all,3)*size(RCM_data_interped.all,4)]);
        RCM_data_interped.mbias=mean(tmp.data_interped_bias,3);
%         pcolor(RCM_data_interped.rmse'); shading flat; colorbar;
%         pcolor(RCM_data_interped.mbias'); shading flat; colorbar;
        [RCM_data_interped.rmse_m, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_interped.rmse, RCM_grid_reana.cut_lon_rho, RCM_grid_reana.cut_lat_rho);
        [RCM_data_interped.mbias_m, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_interped.mbias, RCM_grid_reana.cut_lon_rho, RCM_grid_reana.cut_lat_rho);
        disp([tmp.testname, ' M RMSE = ', num2str(RCM_data_interped.rmse_m)]);
        disp([tmp.testname, ' M bias = ', num2str(RCM_data_interped.mbias_m)]);
        
        tmp.data_interped_all=reshape(RCM_data_interped.all, [size(RCM_data_interped.all,1), size(RCM_data_interped.all,2), size(RCM_data_interped.all,3)*size(RCM_data_interped.all,4)]);
        tmp.data_interped_m=mean(tmp.data_interped_all,3);
        tmp.data_reana_all=reshape(RCM_data_reana.all, [size(RCM_data_interped.all,1), size(RCM_data_interped.all,2), size(RCM_data_interped.all,3)*size(RCM_data_interped.all,4)]);
        tmp.data_reana_m=mean(tmp.data_reana_all,3);
        
        tmp.aaa= tmp.data_interped_m(logical(~isnan(tmp.data_interped_m).*~isnan(tmp.data_reana_m)));
        tmp.bbb= tmp.data_reana_m(logical(~isnan(tmp.data_interped_m).*~isnan(tmp.data_reana_m)));
        [RCM_data_interped.pat_corr, RCM_data_interped.p, RCM_data_interped.corr_l]=corrcoef(tmp.aaa,tmp.bbb, 'alpha', 0.05);
        
        disp([tmp.testname, ' pattern corr = ', num2str(RCM_data_interped.pat_corr(1,2))]);
        
% % %         aaa= movmean(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))),1);
% % %          bbb= movmean(m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))),1);
% % %          corrcoef(aaa,bbb);
% % %          ccc= movmean(m_gcm_interped_ssh(logical(~isnan(m_gcm_interped_ssh).*~isnan(m_cmems_adt))),1);
% % %          ddd= movmean(m_cmems_adt(logical(~isnan(m_gcm_interped_ssh).*~isnan(m_cmems_adt))),1);
% % %          corrcoef(ccc,ddd);
% % % 
% % %          eee= ysecs_m_interped_ssh(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
% % %          fff= ysecs_m_cmems_adt(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
% % %          [ysecscorr, ysecsp, ysecscorr_l, ysecescorr_u]=corrcoef(eee,fff, 'alpha', 0.05)
% % %          ysecscorr_rank = corr(eee,fff, 'Type', 'Spearman')
% % %          ysecscorr_tau = corr(eee,fff, 'Type', 'Kendall')
% % % 
% % %          eeeg= ysecs_m_interped_ssh(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt).*~isnan(ysecs_m_gcm_interped_ssh)));
% % %          fffg= ysecs_m_cmems_adt(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt).*~isnan(ysecs_m_gcm_interped_ssh)));
% % %          [rcmg_ysecscorr, rcmg_ysecsp, rcmg_ysecscorr_l, rcmg_ysecescorr_u]=corrcoef(eeeg,fffg, 'alpha', 0.05)
% % %          rcmg_ysecscorr_rank = corr(eeeg,fffg, 'Type', 'Spearman')
% % %          rcmg_ysecscorr_tau = corr(eeeg,fffg, 'Type', 'Kendall')
        
        
% % %         RCM_grid.mask_model = double(inpolygon(RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
% % %         RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
% % %         RCM_data.mean=RCM_data.mean.*RCM_grid.mask_model;

% % %         m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
% % %         hold on;
% % % 
% % %         m_pcolor(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho',RCM_data.mean');
% % %         shading(gca,param.m_pcolor_shading_method);   
% % %         if(strcmp(tmp.variable, 'SST'))
% % %             [C,h2]=m_contour(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho', RCM_data.mean', [2, 5], 'color','k', ...
% % %                     'linewidth', 1.5, 'linestyle', '-');
% % %                 clabel(C,h2,'FontSize',13,'Color','k', ...
% % %                     'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
% % %         end
% % %  
% % %         m_gshhs_i('color',param.m_gshhs_line_color)  
% % %         m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land
% % % 
% % %         m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
% % %         if min(RCM_info.years) == max(RCM_info.years)
% % %             tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
% % %         else
% % %             tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
% % %         end
% % %         title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title
% % % 
% % %         % set colorbar 
% % %         h = colorbar;
% % %         colormap(jet);
% % % 
% % %         set(h,'fontsize',param.colorbar_fontsize);
% % %         title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
% % %         switch tmp.variable
% % %             case {'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'u', 'v', }
% % %                 caxis(param.colorbar_lev);
% % %             case {'wstrcurl', 'wcurl'}
% % % %                 caxis([-max(abs(RCM_data.mean(:))), max(abs(RCM_data.mean(:)))]);
% % %                 caxis(param.colorbar_lev);
% % %                 colormap(cmaps.byrmap);
% % %         end
% % % 
% % %         disp(['M = ', num2str(tmp.m_value)]);
% % % %         m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 
% % % 
% % %         set(gcf, 'PaperUnits', 'points');
% % %         set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
% % %         set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
% % %         saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
% % %         close all;
% % %         clear RCM_data.mean
% % %         RCM_grid=rmfield(RCM_grid, 'lon_rho');
    end
end