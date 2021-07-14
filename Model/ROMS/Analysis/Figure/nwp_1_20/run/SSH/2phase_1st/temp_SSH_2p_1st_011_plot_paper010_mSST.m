close all; clear all;  clc;
% % horizontal SST trend plot (connected with SSH_mid_report1)
warning off;
% all_region2 ={'AKP','NWP','ES', 'SS', 'YS'}
all_region2 ={'AKP4'};

all_testname2 = {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};

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
    testname=all_testname2{testnameind2}
    inputyear = [1985:2014]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        
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

        load(['D:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);

            valnum=0;
    % %     valid cell number
         for vi=1:size(comb_spatial_meanavhrr,1)
             for vj=1:size(comb_spatial_meanavhrr,2)
                 if (isnan(comb_spatial_meanavhrr(vi,vj,1))~=1)
                    valnum=valnum+1;
                 end
             end
         end

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\','avhrr','\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            avhrrdir='Z:\내 드라이브\Data\Observation\OISST\monthly\';
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        variable='SST';
        run(param_script);

        figdir=[figrawdir,'CLIM\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
            % CLIM AVHRR temp plot
            mean_avhrr=mean(comb_spatial_meanavhrr,3,'omitnan');
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(avhrr_lon),avhrr_lat,mean_avhrr);
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname, ', SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)));  %% + glacier contribution
            titlename = strcat('AVHRR SST ', '(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
            
            [C,h2]=m_contour(double(avhrr_lon),avhrr_lat, mean_avhrr, [10, 15, 20], 'color','k', ...
                    'linewidth', 1.5, 'linestyle', '-');
                clabel(C,h2,'FontSize',13,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'(^oC)','fontsize',colorbar_title_fontsize);
            colorbar_lev = [-2 33];
            caxis(colorbar_lev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            jpgname=strcat(outfile, '_', testname, '_',regionname, '_msst_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            saveas(gcf,jpgname,'jpg'); RemoveWhiteSpace([], 'file', jpgname);

            disp(' ')
            disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
    end
end

% SSH_4th_mid_report6