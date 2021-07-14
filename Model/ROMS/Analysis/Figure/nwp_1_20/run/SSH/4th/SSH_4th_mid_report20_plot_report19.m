close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
all_testname2 = {'test52'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
all_region2 ={'NWP'}

% all_region2 ={'AKP2'};

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
            dropboxpath='C:\Users\KYY\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end

        shadlev = [0 35];
        rms_shadlev = [0 4];
    %     trendlev = [-3 3];  %% trend lev
        trendlev = [-10 10];  %% trend lev
        msllev= [-300 300];
        msl_abs_lev=[-500 1500];
%         msl_abs_lev=[-100 500];
        msl_diff_lev=[-300 300];
        conlev  = 0:5:35;
        meanplotlev =[-0.3 0.3];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1980:2008]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

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
            case('EKB') %% East Korea Bay
                refpolygon=ekbpolygon;
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

        load(['E:\Data\Model\ROMS\nwp_1_20\input\',testname,'\',testname,'_',regionname, ...
            'wind_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

        valnum=0;
%         plot(squeeze(mean(mean(clim_interped_trend_divided,1,'omitnan'),2,'omitnan')))
%         plot(squeeze(mean(mean(clim_recon_trend_divided,1,'omitnan'),2,'omitnan')))     

    % %     valid cell number
    %      for vi=1:size(comb_spatial_meanressh,1)
    %          for vj=1:size(comb_spatial_meanressh,2)
    %              if (isnan(comb_spatial_meanressh(vi,vj,1))~=1)
    %                 valnum=valnum+1;
    %              end
    %          end
    %      end

    
    %     plot(squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan')))
    %     isize = size(comb_interped_data_filtered,1)
    %     jsize = size(comb_interped_data_filtered,2)
    %     lsize = size(comb_interped_data_filtered,3)
    %     comb_yearly_interped_data_filtered=reshape(comb_interped_data_filtered,[isize, jsize, 12, lsize/12]);
    %     mean_yearly_interped_data_filtered=squeeze(mean(mean(mean(comb_yearly_interped_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
    %     trendtime=14:29;
    %     p=polyfit(trendtime,mean_yearly_interped_data_filtered(14:29)',1);
    %     yearly_interped_trend=p(1);
    %     yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_wind_',regionname,'.m']
            filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
%             recondir='E:\Data\Observation\OISST\monthly\';
            run(param_script)
        elseif (strcmp(system_name,'GLNXA64'))
        end
%         model_std = std(mean(mean(comb_interped_data_filtered,'omitnan'),'omitnan'))
%         recon_std = std(mean(mean(comb_recon_data_filtered,'omitnan'),'omitnan'))

        figdir=[figrawdir,'Trend\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        
%         % relative trend plot
%             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%             hold on;
%             m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'-mean_trend_filtered));
%     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
%     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
%             shading(gca,m_pcolor_shading_method);
%             m_gshhs_i('color',m_gshhs_line_color);
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%     %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
%             titlename = strcat('SSH trend(rel), ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)));  %% + glacier contribution
% 
%             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%             % set colorbar 
%             h = colorbar;
%             colormap(bwrmap);
%             set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis(trendlev);
% 
%             % set grid
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%             jpgname=strcat(outfile, '_', testname,'_',regionname, '_relative_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%             saveas(gcf,jpgname,'tif');
% 
%             disp(' ')
%             disp(['clim_', num2str(tempmonth), '_', 'relative_ssh_trend', ' plot is created.'])
%             disp(' ')
%             disp([' File path is : ',jpgname])
%             disp(' ')
% 
%             hold off
%             close all;
 
%             isize = size(comb_interped_data_filtered,1)
%             jsize = size(comb_interped_data_filtered,2)
%             lsize = size(comb_interped_data_filtered,3)
%             comb_monthly_interped_data=reshape(comb_interped_data,[isize, jsize, 12, lsize/12]);
%             comb_yearly_interped_data=squeeze(mean(comb_monthly_interped_data,3,'omitnan'));
            figdir2=[figrawdir,'WIND\'];
            outfile2 = strcat(figdir2,regionname);
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end

        
          for yearij = 1:length(inputyear)
          % % % %             % msl pcolor (absolute)
            
            tempyear= inputyear(yearij);
         
            figdir3=[figrawdir,'WIND\monthly\'];
            outfile3 = strcat(figdir3,regionname);
            if (exist(strcat(figdir3) , 'dir') ~= 7)
                mkdir(strcat(figdir3));
            end
            
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                jpgname=strcat(outfile3, '_', testname,'_',regionname, '_wind_',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'), '.tif'); %% ~_year_month.jpg
                tind=12*(yearij-1)+monthij;
                if (exist(jpgname , 'file') ~= 2)
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    interval=1;
%                     tempU=comb_interped_Udata(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,tind) .*mask_recon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,tind);
%                     tempV=comb_interped_Vdata(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,tind) .*mask_recon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,tind);
                    tempU=comb_interped_Udata(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,tind);
                    tempV=comb_interped_Vdata(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end,tind);
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        [recon_lat2 recon_lon2]=meshgrid(recon_lat, recon_lon);
                        ref_vec_x_ind = find(abs(recon_lon2(1:m_quiver_x_interval:end,1)-m_quiver_ref_text_x_location) == min(abs(recon_lon2(1:m_quiver_x_interval:end,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(recon_lat2(1,1:m_quiver_y_interval:end)-m_quiver_ref_text_y_location) == min(abs(recon_lat2(1,1:m_quiver_y_interval:end)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*1;
%                         ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
%                         ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2));
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));
                    end
                    
                    
                    if (strcmp(regionname,'NWP')==1)
                        emptyind=3;
                        emptyxrange= ref_vec_x_range-emptyind :ref_vec_x_range+emptyind;
                        emptyyrange= ref_vec_y_range-emptyind :ref_vec_y_range+emptyind;
                        tempU(emptyxrange ,emptyyrange)=NaN;
                        tempV(emptyxrange ,emptyyrange)=NaN;   
                    end
                    tempU(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    tempV(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;   
                    uvplot=m_quiver(recon_lon2(1:m_quiver_x_interval:end,1:m_quiver_x_interval:end),recon_lat2(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end),...
                        tempU* m_quiver_vector_size,tempV* m_quiver_vector_size, ...
                        'AutoScale','off','color',m_quiver_vector_color,'LineWidth', m_quiver_LineWidth); 
                  

    %         m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'));
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('Wind',',(',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'),') ');
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                    
                    if (strcmp(regionname,'NWP')==1)
                        m_text(m_quiver_ref_text_x_location-0.2, m_quiver_ref_text_y_location-0.7, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    else
                        m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    end
                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['clim_', num2str(tempmonth), '_', 'absolute_recon_ssh_trend', ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all; 
                end
            end
          end
            
            
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
    end
end