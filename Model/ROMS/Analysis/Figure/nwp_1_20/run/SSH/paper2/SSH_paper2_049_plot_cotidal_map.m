close all; clear all; clc;

dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\t_tide_v1.3beta']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

all_region2 ={'AKP4'};

% all_testname2 = {'test57', 'test58', 'test59', 'test60', 'test65', 'test66', 'test67', 'test68'};
% all_testname2 = {'test57', 'test58', 'test59', 'test60'};
all_testname2 = {'test65', 'test66', 'test67', 'test68'};
% all_testname2 = {'test53', 'test54', 'test55', 'test56'};

% inputyear=[2100];
inputyear=[2006, 2100];

for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        for yearind2=1:length(inputyear)
            close all;
            clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2 inputyear yearind2
            varname='zeta';
            % % %         flag configuration
            for folding=1:1
                fig_flags{1,1}='shade amplitude(M2)';
                fig_flags{2,1}='shade amplitude sum(M2, S2, K1, O1)';
                fig_flags{3,1}='con phase(M2)';
                fig_flags{4,1}='con amplitude(M2)';
                fig_flags{5,1}='SSH amplitude difference(M2) (end-start)';
                fig_flags{6,1}='phase difference(M2) (end-start)';
                fig_flags{7,1}='SSH amplitude difference(S2) (end-start)';
                fig_flags{8,1}='current amplitude difference(M2) (end-start)';
                fig_flags{9,1}='current amplitude (M2) ';
                fig_flags{10,1}='current amplitude (S2) ';
                fig_flags{11,1}='current amplitude difference(S2) (end-start)';
                fig_flags{12,1}='current amplitude (K1) ';
                fig_flags{13,1}='current amplitude (O1) ';
                fig_flags{14,1}='amplitude difference(4con) (end-start)';
                fig_flags{15,1}='M2 tidal validation (line)';
                fig_flags{16,1}='M2 tidal validation (scatter)';
                fig_flags{17,1}='4con tidal validation (line)';
                fig_flags{18,1}='4con tidal validation (scatter)';
                fig_flags{19,1}='tidal validation point';
                fig_flags{20,1}='4con tidal validation (line) - refined';
                fig_flags{21,1}='4con tidal validation (scatter) - refined';
                fig_flags{22,1}='4con tidal validation (scatter) - refined, all test';
                fig_flags{23,1}='amplitude difference(4con) (end-start), all test';
                fig_flags{24,1}='SSH amplitude difference(K1) (end-start)';
                fig_flags{25,1}='SSH amplitude difference(O1) (end-start)';
            end
            for flagi=1:25
                fig_flags{flagi,2}=0;
            end
% 
%             fig_flags{1,2}=0;
%             fig_flags{2,2}=0;
%             fig_flags{3,2}=0;
%             fig_flags{4,2}=0;
%             fig_flags{5,2}=2;
%             fig_flags{6,2}=0;
%             fig_flags{7,2}=2;
%             fig_flags{8,2}=0;
%             fig_flags{9,2}=0;
%             fig_flags{10,2}=0;
%             fig_flags{11,2}=0;
%             fig_flags{12,2}=0;
%             fig_flags{13,2}=0;
%             fig_flags{14,2}=0;
%             fig_flags{15,2}=0;
%             fig_flags{16,2}=0;
%             fig_flags{17,2}=0;
%             fig_flags{18,2}=0;
%             fig_flags{19,2}=0;
%             fig_flags{20,2}=0;
%             fig_flags{21,2}=0;
%             fig_flags{22,2}=0;
%             fig_flags{23,2}=0;
%             fig_flags{24,2}=2;
%             fig_flags{25,2}=2;
            
           fig_flags{1,2}=0;
            fig_flags{2,2}=0;
            fig_flags{3,2}=0;
            fig_flags{4,2}=2;
            fig_flags{5,2}=0;
            fig_flags{6,2}=0;
            fig_flags{7,2}=0;
            fig_flags{8,2}=0;
            fig_flags{9,2}=0;
            fig_flags{10,2}=0;
            fig_flags{11,2}=0;
            fig_flags{12,2}=0;
            fig_flags{13,2}=0;
            fig_flags{14,2}=0;
            fig_flags{15,2}=0;
            fig_flags{16,2}=0;
            fig_flags{17,2}=0;
            fig_flags{18,2}=0;
            fig_flags{19,2}=0;
            fig_flags{20,2}=0;
            fig_flags{21,2}=0;
            fig_flags{22,2}=0;
            fig_flags{23,2}=0;
            fig_flags{24,2}=0;
            fig_flags{25,2}=0;
            
            testname=all_testname2{testnameind2};
            tempyear=inputyear(yearind2);
            outputdir=['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear)];
            run('nwp_polygon_point.m');
            regionname=all_region2{regionind2};
            for folding=1:1
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
                    case('YSECS') %% East China Sea
                        refpolygon=ysecspolygon;
                    case('AKP') %% Around Korea Peninsula
                        refpolygon=akppolygon;
                    case('AKP2') %% Around Korea Peninsula
                        refpolygon=akp2polygon;
                    case('AKP4') %% Around Korea Peninsula
                        refpolygon=akp4polygon;
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

%                 lonfilename=['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear),'\ocean_his_lon_rho_AKP4', '.nc'];
%                 latfilename=['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear),'\ocean_his_lat_rho_AKP4', '.nc'];
                lonfilename=['D:\Data\Model\ROMS\nwp_1_20\ocean_his_lon_rho_AKP4', '.nc'];
                latfilename=['D:\Data\Model\ROMS\nwp_1_20\ocean_his_lat_rho_AKP4', '.nc'];

                if (exist('lon_rho' , 'var') ~= 1)
                    lon_rho=ncread(lonfilename, 'lon_rho');
                    lat_rho=ncread(latfilename, 'lat_rho');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                    cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
            end

            % % for windows
            figrawdir =strcat('D:\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            
            run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
            wrmap = bwrmap(51:100,:);
            byrmap = customcolormap_preset('red-yellow-blue');
            yrmap = byrmap(129:256,:);
            run(param_script);
            
            filename=[filedir, num2str(tempyear,'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
            lon_rho=ncread(filename, 'lon_rho');
            lat_rho=ncread(filename, 'lat_rho');
            
            figdir=[figrawdir,'Tide\', varname, '\'];
            if (exist(strcat(figdir) , 'dir') ~= 7)
                mkdir(strcat(figdir));
            end 
            outfile = strcat(figdir, regionname);
            
            figdir2=[figrawdir,'Tide\', 'current', '\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile2 = strcat(figdir2, regionname);
            
            tide_info.name{1}='M2  ';
            tide_info.name{2}='S2  ';
            tide_info.name{3}='K1  ';
            tide_info.name{4}='O1  ';
            
            tname=ncread(filename, 'tname')';
            num_tide_all=size(tname,1);
            num_tide_tgt=length(tide_info.name);
            for coni=1:num_tide_all
                for tide_namei=1:num_tide_tgt
                    if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                        tide_info.index(tide_namei)=coni;
                    end
                end
            end
            
            tcon=ncread(filename, 'tcon');
            
% start-------------------- m2 amp
            fig_flag=fig_flags{1,2};
            figname=fig_flags{1,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_amp_m2_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)

                    m2_amp=ncread(filename, 'm2_amp');
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(m2_amp'.*100));

                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('M2 amp, ',testname, ',(',num2str(tempyear,'%04i'),')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(wrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
    %                     caxis(abstrendlev);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end

    % start-------------------- 4con amp sum 
            fig_flag=fig_flags{2,2};
            figname=fig_flags{2,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_sum_4con_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_amp=tcon(:,:,tide_info.index(1),1).*100;
                    s2_amp=tcon(:,:,tide_info.index(2),1).*100;
                    k1_amp=tcon(:,:,tide_info.index(3),1).*100;
                    o1_amp=tcon(:,:,tide_info.index(4),1).*100;
                    sum4=m2_amp+s2_amp+k1_amp+o1_amp;
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(sum4'));

                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('4con amp sum, ',testname, ',(',num2str(tempyear,'%04i'),')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(wrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
    %                     caxis(abstrendlev);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 

    % start-------------------- m2 phase con
            fig_flag=fig_flags{3,2};
            figname=fig_flags{3,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_m2_phase_con_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_phase=tcon(:,:,tide_info.index(1),3);

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(m2_phase'));
                    [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(m2_phase'), [0:30:360], 'color','k', ...
                        'linewidth', 1.5, 'linestyle', '-');
    %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    %                         'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                    clabel(C,h2,'FontSize',13,'Color','k', ...
                        'Rotation', 0,'fontweight', 'bold');


                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('m2 phase, ',testname, ',(',num2str(tempyear,'%04i'),')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'^o','fontsize',colorbar_title_fontsize);
    %                     caxis(abstrendlev);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 

    % start-------------------- m2 amp con
            fig_flag=fig_flags{4,2};
            figname=fig_flags{4,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_m2_amp_con_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_amp=tcon(:,:,tide_info.index(1),1).*100;

                    amp1_value= min(m2_amp(80:100,100:140), [], 'all', 'omitnan');
                    amp1_lon=lon_rho(m2_amp==amp1_value);
                    amp1_lat=lat_rho(m2_amp==amp1_value);
                    [amp1_lonind,amp1_latind]=ind2sub([size(m2_amp)],find(m2_amp(:,:)==amp1_value));
                    m2_amp(amp1_lonind-2:amp1_lonind+2, amp1_latind-2:amp1_latind+2)
                    disp(['1st amphidoromic points, ', num2str(amp1_lon), ', ', num2str(amp1_lat)]);
                    amp2_value= min(m2_amp(110:130,170:200), [], 'all', 'omitnan');
                    amp2_lon=lon_rho(m2_amp==amp2_value);
                    amp2_lat=lat_rho(m2_amp==amp2_value);
                    disp(['2nd amphidoromic points, ', num2str(amp2_lon), ', ', num2str(amp2_lat)]);
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(m2_amp'));
                    [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(m2_amp'), [0:20:240], 'color','k', ...
                        'linewidth', 1.5, 'linestyle', '-');
    %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    %                         'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                    clabel(C,h2,'FontSize',13,'Color','k', ...
                        'Rotation', 0,'fontweight', 'bold');


                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('M2 amp, ',testname, ',(',num2str(tempyear,'%04i'),')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(wrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([0 250]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 

    % start-------------------- m2 amp diff
            fig_flag=fig_flags{5,2};
            figname=fig_flags{5,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_m2_amp_diff_', ...
                    num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    tcon_end=ncread(endfilename, 'tcon');
                    startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    tcon_start=ncread(startfilename, 'tcon');

                    m2_amp_end=tcon_end(:,:,tide_info.index(1),1).*100;
                    m2_amp_start=tcon_start(:,:,tide_info.index(1),1).*100;
                    m2_amp_diff=m2_amp_end-m2_amp_start;

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(m2_amp_diff'));
    % %                     [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(m2_amp_diff'), [0:20:240], 'color','k', ...
    % %                         'linewidth', 1.5, 'linestyle', '-');
    % %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    % %                         'Rotation', 0,'fontweight', 'bold');

                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('M2 amp diff, ',testname, ',(',num2str(inputyear(end),'%04i'),'-',num2str(inputyear(1),'%04i'), ')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([-10 10]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');
                    save([filedir,testname,'_diff_m2_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'], 'm2_amp_diff');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 

    % start-------------------- s2 amp con
            fig_flag=fig_flags{6,2};
            figname=fig_flags{6,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_s2_amp_con_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    s2_amp=tcon(:,:,tide_info.index(2),1).*100;

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(s2_amp'));
                    [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(s2_amp'), [0:20:240], 'color','k', ...
                        'linewidth', 1.5, 'linestyle', '-');
    %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    %                         'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                    clabel(C,h2,'FontSize',13,'Color','k', ...
                        'Rotation', 0,'fontweight', 'bold');


                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('M2 amp, ',testname, ',(',num2str(tempyear,'%04i'),')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(wrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([0 250]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 

    % start-------------------- s2 amp diff
            fig_flag=fig_flags{7,2};
            figname=fig_flags{7,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_s2_amp_diff_', ...
                    num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    tcon_end=ncread(endfilename, 'tcon');
                    startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    tcon_start=ncread(startfilename, 'tcon');

                    s2_amp_end=tcon_end(:,:,tide_info.index(2),1).*100;
                    s2_amp_start=tcon_start(:,:,tide_info.index(2),1).*100;
                    s2_amp_diff=s2_amp_end-s2_amp_start;

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(s2_amp_diff'));
    % %                     [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(s2_amp_diff'), [0:20:240], 'color','k', ...
    % %                         'linewidth', 1.5, 'linestyle', '-');
    % %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    % %                         'Rotation', 0,'fontweight', 'bold');

                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('S2 amp diff, ',testname, ',(',num2str(inputyear(end),'%04i'),'-',num2str(inputyear(1),'%04i'), ')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([-10 10]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');
                    save([filedir,testname,'_diff_s2_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'], 's2_amp_diff');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 

    % start-------------------- m2 current diff plot
            fig_flag=fig_flags{8,2};
            while (fig_flag)
                jpgname=strcat(outfile2, '_', testname,'_',regionname, '_', '_m2_current_diff_', ...
                        num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)    

                    ubar_endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', 'harmonic_analysis_', 'ubar_eastward', '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    ubar_tcon_end=ncread(ubar_endfilename, 'tcon');
                    ubar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'ubar_eastward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    ubar_tcon_start=ncread(ubar_startfilename, 'tcon');

                    vbar_endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', 'harmonic_analysis_', 'vbar_northward', '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    vbar_tcon_end=ncread(vbar_endfilename, 'tcon');
                    vbar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'vbar_northward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    vbar_tcon_start=ncread(vbar_startfilename, 'tcon');

                    ubar_m2_amp_end=ubar_tcon_end(:,:,tide_info.index(1),1);
                    ubar_m2_amp_start=ubar_tcon_start(:,:,tide_info.index(1),1);
                    u_rho=ubar_m2_amp_end-ubar_m2_amp_start;

                    vbar_m2_amp_end=vbar_tcon_end(:,:,tide_info.index(1),1);
                    vbar_m2_amp_start=vbar_tcon_start(:,:,tide_info.index(1),1);
                    v_rho=vbar_m2_amp_end-vbar_m2_amp_start;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho = u_rho.*mask_model;
                    v_rho = v_rho.*mask_model;
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value /10.0;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value /10.0;     

                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size *10, ...
                                    v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size *10, ...
                                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);
                    titlename = strcat('M2 current diff, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,jpgname,'tif');
                    close all;
                    clear lon_rho mean_u ref_vec_x_range
                end
                fig_flag=0;
            end

    % start-------------------- m2 current plot
            fig_flag=fig_flags{9,2};
            while (fig_flag)
                jpgname=strcat(outfile2, '_', testname,'_',regionname, '_', '_m2_current_', ...
                        num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)    

                    ubar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'ubar_eastward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    ubar_tcon_start=ncread(ubar_startfilename, 'tcon');

                    vbar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'vbar_northward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    vbar_tcon_start=ncread(vbar_startfilename, 'tcon');

                    ubar_m2_amp_start=ubar_tcon_start(:,:,tide_info.index(1),1);
                    u_rho=ubar_m2_amp_start;

                    vbar_m2_amp_start=vbar_tcon_start(:,:,tide_info.index(1),1);
                    v_rho=vbar_m2_amp_start;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho = u_rho.*mask_model;
                    v_rho = v_rho.*mask_model;
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     

                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);
                    titlename = strcat('M2 current, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,jpgname,'tif');
                    close all;
                    clear lon_rho mean_u ref_vec_x_range
                end
                fig_flag=0;
            end

    % start-------------------- S2 current plot
            fig_flag=fig_flags{10,2};
            while (fig_flag)
                jpgname=strcat(outfile2, '_', testname,'_',regionname, '_', '_s2_current_', ...
                        num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)    

                    ubar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'ubar_eastward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    ubar_tcon_start=ncread(ubar_startfilename, 'tcon');

                    vbar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'vbar_northward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    vbar_tcon_start=ncread(vbar_startfilename, 'tcon');

                    ubar_m2_amp_start=ubar_tcon_start(:,:,tide_info.index(2),1);
                    u_rho=ubar_m2_amp_start;

                    vbar_m2_amp_start=vbar_tcon_start(:,:,tide_info.index(2),1);
                    v_rho=vbar_m2_amp_start;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho = u_rho.*mask_model;
                    v_rho = v_rho.*mask_model;
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     

                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);
                    titlename = strcat('S2 current, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,jpgname,'tif');
                    close all;
                    clear lon_rho mean_u ref_vec_x_range
                end
                fig_flag=0;
            end

    % start-------------------- S2 current diff plot
            fig_flag=fig_flags{11,2};
            while (fig_flag)
                jpgname=strcat(outfile2, '_', testname,'_',regionname, '_', '_s2_current_diff_', ...
                        num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)    

                    ubar_endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', 'harmonic_analysis_', 'ubar_eastward', '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    ubar_tcon_end=ncread(ubar_endfilename, 'tcon');
                    ubar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'ubar_eastward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    ubar_tcon_start=ncread(ubar_startfilename, 'tcon');

                    vbar_endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', 'harmonic_analysis_', 'vbar_northward', '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    vbar_tcon_end=ncread(vbar_endfilename, 'tcon');
                    vbar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'vbar_northward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    vbar_tcon_start=ncread(vbar_startfilename, 'tcon');

                    ubar_m2_amp_end=ubar_tcon_end(:,:,tide_info.index(2),1);
                    ubar_m2_amp_start=ubar_tcon_start(:,:,tide_info.index(2),1);
                    u_rho=ubar_m2_amp_end-ubar_m2_amp_start;

                    vbar_m2_amp_end=vbar_tcon_end(:,:,tide_info.index(2),1);
                    vbar_m2_amp_start=vbar_tcon_start(:,:,tide_info.index(2),1);
                    v_rho=vbar_m2_amp_end-vbar_m2_amp_start;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho = u_rho.*mask_model;
                    v_rho = v_rho.*mask_model;
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value/10.0;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value/10.0;     

                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size *10, ...
                                    v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size *10, ...
                                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);
                    titlename = strcat('S2 current diff, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,jpgname,'tif');
                    close all;
                    clear lon_rho mean_u ref_vec_x_range
                end
                fig_flag=0;
            end

    % start-------------------- K1 current plot
            fig_flag=fig_flags{12,2};
            while (fig_flag)
                jpgname=strcat(outfile2, '_', testname,'_',regionname, '_', '_k1_current_', ...
                        num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)    

                    ubar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'ubar_eastward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    ubar_tcon_start=ncread(ubar_startfilename, 'tcon');

                    vbar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'vbar_northward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    vbar_tcon_start=ncread(vbar_startfilename, 'tcon');

                    ubar_m2_amp_start=ubar_tcon_start(:,:,tide_info.index(3),1);
                    u_rho=ubar_m2_amp_start;

                    vbar_m2_amp_start=vbar_tcon_start(:,:,tide_info.index(3),1);
                    v_rho=vbar_m2_amp_start;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho = u_rho.*mask_model;
                    v_rho = v_rho.*mask_model;
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     

                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);
                    titlename = strcat('K1 current, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,jpgname,'tif');
                    close all;
                    clear lon_rho mean_u ref_vec_x_range
                end
                fig_flag=0;
            end

    % start-------------------- O1 current plot
            fig_flag=fig_flags{13,2};
            while (fig_flag)
                jpgname=strcat(outfile2, '_', testname,'_',regionname, '_', '_o1_current_', ...
                        num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)    

                    ubar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'ubar_eastward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    ubar_tcon_start=ncread(ubar_startfilename, 'tcon');

                    vbar_startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', 'vbar_northward', '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    vbar_tcon_start=ncread(vbar_startfilename, 'tcon');

                    ubar_m2_amp_start=ubar_tcon_start(:,:,tide_info.index(4),1);
                    u_rho=ubar_m2_amp_start;

                    vbar_m2_amp_start=vbar_tcon_start(:,:,tide_info.index(4),1);
                    v_rho=vbar_m2_amp_start;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
                    u_rho = u_rho.*mask_model;
                    v_rho = v_rho.*mask_model;
                    if (exist('ref_vec_x_range' , 'var') ~= 1)
                        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
                        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
                    end
                    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
                    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     

                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;
                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                    u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

                    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);
                    titlename = strcat('O1 current, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,jpgname,'tif');
                    close all;
                    clear lon_rho mean_u ref_vec_x_range
                end
                fig_flag=0;
            end

% start-------------------- 4con amp diff
            fig_flag=fig_flags{14,2};
            figname=fig_flags{14,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_4con_amp_diff_', ...
                    num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    tcon_end=ncread(endfilename, 'tcon');
                    startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    tcon_start=ncread(startfilename, 'tcon');

                    m2_amp_end=tcon_end(:,:,tide_info.index(1),1).*100;
                    s2_amp_end=tcon_end(:,:,tide_info.index(2),1).*100;
                    k1_amp_end=tcon_end(:,:,tide_info.index(3),1).*100;
                    o1_amp_end=tcon_end(:,:,tide_info.index(4),1).*100;
                    con4_amp_end= m2_amp_end + s2_amp_end + k1_amp_end + o1_amp_end;

                    m2_amp_start=tcon_start(:,:,tide_info.index(1),1).*100;
                    s2_amp_start=tcon_end(:,:,tide_info.index(2),1).*100;
                    k1_amp_start=tcon_end(:,:,tide_info.index(3),1).*100;
                    o1_amp_start=tcon_end(:,:,tide_info.index(4),1).*100;
                    con4_amp_start= m2_amp_start + s2_amp_start + k1_amp_start + o1_amp_start;

                    con4_amp_diff=con4_amp_end-con4_amp_start;
                    
%                     mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%                     mask_model(mask_model==0)=NaN;    
%                     con4_amp_diff=con4_amp_diff.*mask_model;
                    
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
%                     ax1=axes;
                    plot1=m_pcolor(double(lon_rho)',lat_rho',squeeze(con4_amp_diff'));
%                     plot2=m_pcolor(double(lon_rho)',lat_rho',squeeze(con4_amp_diff'));

%                     m_surface(double(lon_rho)',lat_rho',squeeze(con4_amp_diff'),[0.8 0.8 0.8])

    % %                     [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(m2_amp_diff'), [0:20:240], 'color','k', ...
    % %                         'linewidth', 1.5, 'linestyle', '-');
    % %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    % %                         'Rotation', 0,'fontweight', 'bold');
                    shading(gca,m_pcolor_shading_method);
%                     set(plot1, 'facecolor', 'r');
%                     m_gshhs_i('color',m_gshhs_line_color);
%                     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%                     set(gca1, 'color', [0.8,0.8,0.8]);

%                     surface('Tag','m_grid_color','Parent',Parent1,'ZData',ZData1,'YData',YData1,...
%                         'XData',XData1,...
%                         'LineStyle','none',...
%                         'FaceColor',[0.8 0.8 0.8]);
%                     m_patch([130 140  140 130], [30 30 35 35], 'w', 'LineStyle', 'none')

                    titlename = strcat('4con amp diff, ',testname, ',(',num2str(inputyear(end),'%04i'),'-',num2str(inputyear(1),'%04i'), ')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([-15 15]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);
%                     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');
                    save([filedir,testname,'_diff_4con_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'], 'con4_amp_diff');
                    
                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 

% start-------------------- m2 amp validiation plot
            fig_flag=fig_flags{15,2};
            figname=fig_flags{15,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_m2_amp_val_line_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_amp=tcon(:,:,tide_info.index(1),1).*100;

                    opts = spreadsheetImportOptions("NumVariables", 12);
                    opts.Sheet = "Sheet1";
                    opts.DataRange = "B2:M50";
                    opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
                    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
                    tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data.xlsx", opts, "UseExcel", false);
                    clear opts

    %                 if (exist('tcon_model') ~= 1)
                    tconname=[filedir, num2str(tempyear,'%04i'), '\', testname, '_tcon_val_zeta_', num2str(tempyear, '%04i'), '.mat'];
                    if (exist(tconname , 'file') ~= 2)
                        tcon_model=cell(size(tcon_obs));
                        tcon_model(:,1)=table2cell(tcon_obs(:,1));
                        tcon_model(:,2)=table2cell(tcon_obs(:,2));
                        for ind_sta=1:size(tcon_obs,1)
                            tgtlat=table2array(tcon_obs(ind_sta,3));
                            tgtlon=table2array(tcon_obs(ind_sta,4));

                            for i=1:size(lon_rho,1)
                                for j=1:size(lat_rho,2)
                                    if isfinite(tcon(i,j,tide_info.index(1),1))
                                        dist(i,j)=m_lldist([lon_rho(i,j), tgtlon], [lat_rho(i,j), tgtlat]);
                                    else
                                        dist(i,j)=NaN;
                                    end
                                end
                            end

                            ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
                            ind_lon=mod(ind_sta_model(ind_sta),size(lon_rho,1));
                            ind_lat=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
                            tcon_model{ind_sta,4}=ind_lon;
                            tcon_model{ind_sta,3}=ind_lat;
                            tcon_model{ind_sta,5}=tcon(ind_lon,ind_lat,tide_info.index(1),1).*100;
                            tcon_model{ind_sta,6}=tcon(ind_lon,ind_lat,tide_info.index(2),1).*100;
                            tcon_model{ind_sta,7}=tcon(ind_lon,ind_lat,tide_info.index(3),1).*100;
                            tcon_model{ind_sta,8}=tcon(ind_lon,ind_lat,tide_info.index(4),1).*100;
                            tcon_model{ind_sta,9}=tcon(ind_lon,ind_lat,tide_info.index(1),3).*100;
                            tcon_model{ind_sta,10}=tcon(ind_lon,ind_lat,tide_info.index(2),3).*100;
                            tcon_model{ind_sta,11}=tcon(ind_lon,ind_lat,tide_info.index(3),3).*100;
                            tcon_model{ind_sta,12}=tcon(ind_lon,ind_lat,tide_info.index(4),3).*100;
                        end
                        save(tconname, 'tcon_model')
                    else
                        load(tconname)
                    end

                    num_sta=size(tcon_obs,1);
                    mslplot{1}=plot(cell2mat(tcon_model(:,5)), 'b');
                    hold on
                    mslplot{2}=plot(table2array(tcon_obs(:,5)), 'k');
                    hold off
                    xticks(1:num_sta)
                    xticklabels(tcon_model(:,2))
                    xtickangle(45)

                    xlabel('Station')
                    ylabel('M2 Tidal amplitude (cm)')
                    title([testname, ', M2 amp, ', num2str(tempyear,'%04i')])
                    % datetick('x','yyyy','keepticks')
                    axis tight;
                    % ylim(meanplotlev2)
                    set(mslplot{1},'LineWidth',2);
                    set(mslplot{2},'LineWidth',2);

                    set(gca,'FontSize',15);
                    grid on
                    lgd=legend('Model','TG-Choi');
                    set(lgd,'FontSize',15);
                    set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                    set(lgd,'Orientation','horizontal');
                    % for tind=1:size(RCM_interped_sla_yearly_mean,2)
                    %     model_std(tind)=std(RCM_interped_sla_yearly_mean(:,tind));
                    % end
                    % meanstd=mean(model_std);
                    % txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 

                    set(gcf,'PaperPosition', [0 0 36 12]) 
                    % hold off
                    saveas(gcf,jpgname,'tif');
                    grid off
                    close all;
                end
                fig_flag=0;
            end 

% start-------------------- m2 amp validiation scatter
            fig_flag=fig_flags{16,2};
            figname=fig_flags{16,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_m2_amp_val_scatter_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_amp=tcon(:,:,tide_info.index(1),1).*100;

                    opts = spreadsheetImportOptions("NumVariables", 12);
                    opts.Sheet = "Sheet1";
                    opts.DataRange = "B2:M50";
                    opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
                    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
                    tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data.xlsx", opts, "UseExcel", false);
                    clear opts

    %                 if (exist('tcon_model') ~= 1)
                    tconname=[filedir, num2str(tempyear,'%04i'), '\', testname, '_tcon_val_zeta_', num2str(tempyear, '%04i'), '.mat'];
                    if (exist(tconname , 'file') ~= 2)
                        tcon_model=cell(size(tcon_obs));
                        tcon_model(:,1)=table2cell(tcon_obs(:,1));
                        tcon_model(:,2)=table2cell(tcon_obs(:,2));
                        for ind_sta=1:size(tcon_obs,1)
                            tgtlat=table2array(tcon_obs(ind_sta,3));
                            tgtlon=table2array(tcon_obs(ind_sta,4));

                            for i=1:size(lon_rho,1)
                                for j=1:size(lat_rho,2)
                                    if isfinite(tcon(i,j,tide_info.index(1),1))
                                        dist(i,j)=m_lldist([lon_rho(i,j), tgtlon], [lat_rho(i,j), tgtlat]);
                                    else
                                        dist(i,j)=NaN;
                                    end
                                end
                            end

                            ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
                            ind_lon=mod(ind_sta_model(ind_sta),size(lon_rho,1));
                            ind_lat=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
                            tcon_model{ind_sta,4}=ind_lon;
                            tcon_model{ind_sta,3}=ind_lat;
                            tcon_model{ind_sta,5}=tcon(ind_lon,ind_lat,tide_info.index(1),1).*100;
                            tcon_model{ind_sta,6}=tcon(ind_lon,ind_lat,tide_info.index(2),1).*100;
                            tcon_model{ind_sta,7}=tcon(ind_lon,ind_lat,tide_info.index(3),1).*100;
                            tcon_model{ind_sta,8}=tcon(ind_lon,ind_lat,tide_info.index(4),1).*100;
                            tcon_model{ind_sta,9}=tcon(ind_lon,ind_lat,tide_info.index(1),3).*100;
                            tcon_model{ind_sta,10}=tcon(ind_lon,ind_lat,tide_info.index(2),3).*100;
                            tcon_model{ind_sta,11}=tcon(ind_lon,ind_lat,tide_info.index(3),3).*100;
                            tcon_model{ind_sta,12}=tcon(ind_lon,ind_lat,tide_info.index(4),3).*100;
                        end
                        save(tconname, 'tcon_model')
                    else
                        load(tconname)
                    end



    %                 num_sta=size(tcon_obs,1);
                    mslplot{1}=scatter(table2array(tcon_obs(:,5)), cell2mat(tcon_model(:,5)), 'x');
                    hold on
                    maxval=max([table2array(tcon_obs(:,5)), cell2mat(tcon_model(:,5))]);
                    mslplot{2}=plot(0:maxval, 0:maxval);
                    hold off

                    xlabel('M2 Observed Tidal amplitude (cm)')
                    ylabel('M2 Model Tidal amplitude (cm)')
                    title([testname, ', M2 amp, ', num2str(tempyear,'%04i')])
                    % datetick('x','yyyy','keepticks')
                    axis tight;
                    % ylim(meanplotlev2)
                    set(mslplot{1},'LineWidth',2);
                    set(mslplot{2},'LineWidth',2);

                    set(gca,'FontSize',15);
                    grid on
    %                 lgd=legend('Model','TG-Choi');
    %                 set(lgd,'FontSize',15);
    %                 set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
    %                 set(lgd,'Orientation','horizontal');

                    set(gcf,'PaperPosition', [0 0 36 12]) 
                    % hold off
                    saveas(gcf,jpgname,'tif');
                    grid off
                    close all;
                end
                fig_flag=0;
            end 

% start-------------------- 4con amp validiation plot
            fig_flag=fig_flags{17,2};
            figname=fig_flags{17,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_4con_amp_val_line_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_amp=tcon(:,:,tide_info.index(1),1).*100;

                    opts = spreadsheetImportOptions("NumVariables", 12);
                    opts.Sheet = "Sheet1";
                    opts.DataRange = "B2:M50";
                    opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
                    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
                    tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data.xlsx", opts, "UseExcel", false);
                    clear opts

    %                 if (exist('tcon_model') ~= 1)
                    tconname=[filedir, num2str(tempyear,'%04i'), '\', testname, '_tcon_val_zeta_', num2str(tempyear, '%04i'), '.mat'];
                    if (exist(tconname , 'file') ~= 2)
                        tcon_model=cell(size(tcon_obs));
                        tcon_model(:,1)=table2cell(tcon_obs(:,1));
                        tcon_model(:,2)=table2cell(tcon_obs(:,2));
                        for ind_sta=1:size(tcon_obs,1)
                            tgtlat=table2array(tcon_obs(ind_sta,3));
                            tgtlon=table2array(tcon_obs(ind_sta,4));

                            for i=1:size(lon_rho,1)
                                for j=1:size(lat_rho,2)
                                    if isfinite(tcon(i,j,tide_info.index(1),1))
                                        dist(i,j)=m_lldist([lon_rho(i,j), tgtlon], [lat_rho(i,j), tgtlat]);
                                    else
                                        dist(i,j)=NaN;
                                    end
                                end
                            end

                            ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
                            ind_lon=mod(ind_sta_model(ind_sta),size(lon_rho,1));
                            ind_lat=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
                            tcon_model{ind_sta,4}=ind_lon;
                            tcon_model{ind_sta,3}=ind_lat;
                            tcon_model{ind_sta,5}=tcon(ind_lon,ind_lat,tide_info.index(1),1).*100;
                            tcon_model{ind_sta,6}=tcon(ind_lon,ind_lat,tide_info.index(2),1).*100;
                            tcon_model{ind_sta,7}=tcon(ind_lon,ind_lat,tide_info.index(3),1).*100;
                            tcon_model{ind_sta,8}=tcon(ind_lon,ind_lat,tide_info.index(4),1).*100;
                            tcon_model{ind_sta,9}=tcon(ind_lon,ind_lat,tide_info.index(1),3).*100;
                            tcon_model{ind_sta,10}=tcon(ind_lon,ind_lat,tide_info.index(2),3).*100;
                            tcon_model{ind_sta,11}=tcon(ind_lon,ind_lat,tide_info.index(3),3).*100;
                            tcon_model{ind_sta,12}=tcon(ind_lon,ind_lat,tide_info.index(4),3).*100;
                        end
                        save(tconname, 'tcon_model')
                    else
                        load(tconname)
                    end

                    model_4con= cell2mat(tcon_model(:,5)) + cell2mat(tcon_model(:,6)) + cell2mat(tcon_model(:,7)) + cell2mat(tcon_model(:,8));
                    obs_4con =  table2array(tcon_obs(:,5)) + table2array(tcon_obs(:,6)) + table2array(tcon_obs(:,7)) + table2array(tcon_obs(:,8));
                    num_sta=size(tcon_obs,1);
                    mslplot{1}=plot(model_4con, 'b');
                    hold on
                    mslplot{2}=plot(obs_4con, 'k');
                    hold off
                    xticks(1:num_sta)
                    xticklabels(tcon_model(:,2))
                    xtickangle(45)

                    xlabel('Station')
                    ylabel('4con Tidal amplitude (cm)')
                    title([testname, ', 4con amp, ', num2str(tempyear,'%04i')])
                    % datetick('x','yyyy','keepticks')
                    axis tight;
                    % ylim(meanplotlev2)
                    set(mslplot{1},'LineWidth',2);
                    set(mslplot{2},'LineWidth',2);

                    set(gca,'FontSize',15);
                    grid on
                    lgd=legend('Model','TG-Choi');
                    set(lgd,'FontSize',15);
                    set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                    set(lgd,'Orientation','horizontal');
                    % for tind=1:size(RCM_interped_sla_yearly_mean,2)
                    %     model_std(tind)=std(RCM_interped_sla_yearly_mean(:,tind));
                    % end
                    % meanstd=mean(model_std);
                    % txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 

                    set(gcf,'PaperPosition', [0 0 36 12]) 
                    % hold off
                    saveas(gcf,jpgname,'tif');
                    grid off
                    close all;
                end
                fig_flag=0;
            end 

% start-------------------- 4con amp validiation scatter
            fig_flag=fig_flags{18,2};
            figname=fig_flags{18,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_4con_amp_val_scatter_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_amp=tcon(:,:,tide_info.index(1),1).*100;

                    opts = spreadsheetImportOptions("NumVariables", 12);
                    opts.Sheet = "Sheet1";
                    opts.DataRange = "B2:M50";
                    opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
                    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
                    tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data.xlsx", opts, "UseExcel", false);
                    clear opts

    %                 if (exist('tcon_model') ~= 1)
                    tconname=[filedir, num2str(tempyear,'%04i'), '\', testname, '_tcon_val_zeta_', num2str(tempyear, '%04i'), '.mat'];
                    if (exist(tconname , 'file') ~= 2)
                        tcon_model=cell(size(tcon_obs));
                        tcon_model(:,1)=table2cell(tcon_obs(:,1));
                        tcon_model(:,2)=table2cell(tcon_obs(:,2));
                        for ind_sta=1:size(tcon_obs,1)
                            tgtlat=table2array(tcon_obs(ind_sta,3));
                            tgtlon=table2array(tcon_obs(ind_sta,4));

                            for i=1:size(lon_rho,1)
                                for j=1:size(lat_rho,2)
                                    if isfinite(tcon(i,j,tide_info.index(1),1))
                                        dist(i,j)=m_lldist([lon_rho(i,j), tgtlon], [lat_rho(i,j), tgtlat]);
                                    else
                                        dist(i,j)=NaN;
                                    end
                                end
                            end

                            ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
                            ind_lon=mod(ind_sta_model(ind_sta),size(lon_rho,1));
                            ind_lat=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
                            tcon_model{ind_sta,4}=ind_lon;
                            tcon_model{ind_sta,3}=ind_lat;
                            tcon_model{ind_sta,5}=tcon(ind_lon,ind_lat,tide_info.index(1),1).*100;
                            tcon_model{ind_sta,6}=tcon(ind_lon,ind_lat,tide_info.index(2),1).*100;
                            tcon_model{ind_sta,7}=tcon(ind_lon,ind_lat,tide_info.index(3),1).*100;
                            tcon_model{ind_sta,8}=tcon(ind_lon,ind_lat,tide_info.index(4),1).*100;
                            tcon_model{ind_sta,9}=tcon(ind_lon,ind_lat,tide_info.index(1),3).*100;
                            tcon_model{ind_sta,10}=tcon(ind_lon,ind_lat,tide_info.index(2),3).*100;
                            tcon_model{ind_sta,11}=tcon(ind_lon,ind_lat,tide_info.index(3),3).*100;
                            tcon_model{ind_sta,12}=tcon(ind_lon,ind_lat,tide_info.index(4),3).*100;
                        end
                        save(tconname, 'tcon_model')
                    else
                        load(tconname)
                    end


                    model_4con= cell2mat(tcon_model(:,5)) + cell2mat(tcon_model(:,6)) + cell2mat(tcon_model(:,7)) + cell2mat(tcon_model(:,8));
                    obs_4con =  table2array(tcon_obs(:,5)) + table2array(tcon_obs(:,6)) + table2array(tcon_obs(:,7)) + table2array(tcon_obs(:,8));

    %                 num_sta=size(tcon_obs,1);
                    mslplot{1}=scatter(obs_4con, model_4con, 'x')
                    hold on
                    maxval=max([obs_4con,  model_4con]);
                    mslplot{2}=plot(0:maxval, 0:maxval);
                    hold off

                    xlabel('4con Observed Tidal amplitude (cm)')
                    ylabel('4con Model Tidal amplitude (cm)')
                    title([testname, ', 4con amp, ', num2str(tempyear,'%04i')])
                    % datetick('x','yyyy','keepticks')
                    axis tight;
                    % ylim(meanplotlev2)
                    set(mslplot{1},'LineWidth',2);
                    set(mslplot{2},'LineWidth',2);

                    set(gca,'FontSize',15);
                    grid on
    %                 lgd=legend('Model','TG-Choi');
    %                 set(lgd,'FontSize',15);
    %                 set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
    %                 set(lgd,'Orientation','horizontal');

                    set(gcf,'PaperPosition', [0 0 36 12]) 
                    % hold off
                    saveas(gcf,jpgname,'tif');
                    grid off
                    close all;
                end
                fig_flag=0;
            end 

% start-------------------- validiation point dot
            fig_flag=fig_flags{19,2};
            figname=fig_flags{19,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_amp_val_point_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    opts = spreadsheetImportOptions("NumVariables", 12);
                    opts.Sheet = "Sheet1";
                    opts.DataRange = "B2:M20";
                    opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
                    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
                    tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data4_obs_and_khoa.xlsx", opts, "UseExcel", false);
                    clear opts

                    tconname=[filedir, num2str(tempyear,'%04i'), '\', testname, '_tcon_val_zeta_', num2str(tempyear, '%04i'), '.mat'];
                    load(tconname)
                    
                    amplev=([0 420]);

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    hold on;
                    amp_4con=sum(table2array(tcon_obs(:,5:8)),2);
                    for linei=1:size(tcon_obs,1)
                        if (amp_4con(linei)<=min(amplev))
                            colind=1;
                        elseif (amp_4con(linei)>=max(amplev))
                            colind=128;
                        else
                            colind=round((amp_4con(linei)-min(amplev))/diff(amplev) *128.0);
                        end
                        m_line(table2array(tcon_obs(linei,4)),table2array(tcon_obs(linei,3)),'marker','o','color',yrmap(colind,:),'linewi',2,...
                          'linest','none','markersize',8,'markerfacecolor',yrmap(colind,:));
                    end
%                     m_plot(table2array(tcon_obs(:,4)),table2array(tcon_obs(:,3)),'marker','o','color','k', ...
%                     'markerfacecolor',bwrmap(colind,:)); 
                    h = colorbar;
                    colormap(yrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
                    caxis(amplev);

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat('Obs station');  %% + glacier contribution
%                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;

                end
            fig_flag=0;
            end 

% start-------------------- 4con amp validiation plot -refined
            fig_flag=fig_flags{20,2};
            figname=fig_flags{20,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_4con_amp_val_ref_line_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_amp=tcon(:,:,tide_info.index(1),1).*100;

                    opts = spreadsheetImportOptions("NumVariables", 12);
                    opts.Sheet = "Sheet1";
                    opts.DataRange = "B2:M20";
                    opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
                    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
                    tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data4_obs_and_khoa.xlsx", opts, "UseExcel", false);
                    clear opts

    %                 if (exist('tcon_model') ~= 1)
                    tconname=[filedir, num2str(tempyear,'%04i'), '\', testname, '_tcon_val_zeta_ref_', num2str(tempyear, '%04i'), '.mat'];
%                     if (exist(tconname , 'file') ~= 2)
                        tcon_model=cell(size(tcon_obs));
                        tcon_model(:,1)=table2cell(tcon_obs(:,1));
                        tcon_model(:,2)=table2cell(tcon_obs(:,2));
                        for ind_sta=1:size(tcon_obs,1)
                            tgtlat=table2array(tcon_obs(ind_sta,3));
                            tgtlon=table2array(tcon_obs(ind_sta,4));

                            for i=1:size(lon_rho,1)
                                for j=1:size(lat_rho,2)
                                    if isfinite(tcon(i,j,tide_info.index(1),1))
                                        dist(i,j)=m_lldist([lon_rho(i,j), tgtlon], [lat_rho(i,j), tgtlat]);
                                    else
                                        dist(i,j)=NaN;
                                    end
                                end
                            end

                            ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
                            ind_lon=mod(ind_sta_model(ind_sta),size(lon_rho,1));
                            ind_lat=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
                            tcon_model{ind_sta,4}=ind_lon;
                            tcon_model{ind_sta,3}=ind_lat;
                            tcon_model{ind_sta,5}=tcon(ind_lon,ind_lat,tide_info.index(1),1).*100;
                            tcon_model{ind_sta,6}=tcon(ind_lon,ind_lat,tide_info.index(2),1).*100;
                            tcon_model{ind_sta,7}=tcon(ind_lon,ind_lat,tide_info.index(3),1).*100;
                            tcon_model{ind_sta,8}=tcon(ind_lon,ind_lat,tide_info.index(4),1).*100;
                            tcon_model{ind_sta,9}=tcon(ind_lon,ind_lat,tide_info.index(1),3).*100;
                            tcon_model{ind_sta,10}=tcon(ind_lon,ind_lat,tide_info.index(2),3).*100;
                            tcon_model{ind_sta,11}=tcon(ind_lon,ind_lat,tide_info.index(3),3).*100;
                            tcon_model{ind_sta,12}=tcon(ind_lon,ind_lat,tide_info.index(4),3).*100;
                        end
                        save(tconname, 'tcon_model')
%                     else
                        load(tconname)
%                     end

                    model_4con= cell2mat(tcon_model(:,5)) + cell2mat(tcon_model(:,6)) + cell2mat(tcon_model(:,7)) + cell2mat(tcon_model(:,8));
                    obs_4con =  table2array(tcon_obs(:,5)) + table2array(tcon_obs(:,6)) + table2array(tcon_obs(:,7)) + table2array(tcon_obs(:,8));
                    num_sta=size(tcon_obs,1);
                    mslplot{1}=plot(model_4con, 'b');
                    hold on
                    mslplot{2}=plot(obs_4con, 'k');
                    hold off
                    xticks(1:num_sta)
                    xticklabels(tcon_model(:,2))
                    xtickangle(45)

                    xlabel('Station')
                    ylabel('4con Tidal amplitude (cm)')
                    title([testname, ', 4con amp, ', num2str(tempyear,'%04i')])
                    % datetick('x','yyyy','keepticks')
                    axis tight;
                    % ylim(meanplotlev2)
                    set(mslplot{1},'LineWidth',2);
                    set(mslplot{2},'LineWidth',2);

                    set(gca,'FontSize',15);
                    grid on
                    lgd=legend('Model','TG-Choi');
                    set(lgd,'FontSize',15);
                    set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                    set(lgd,'Orientation','horizontal');
                    % for tind=1:size(RCM_interped_sla_yearly_mean,2)
                    %     model_std(tind)=std(RCM_interped_sla_yearly_mean(:,tind));
                    % end
                    % meanstd=mean(model_std);
                    % txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 

                    set(gcf,'PaperPosition', [0 0 36 12]) 
                    % hold off
                    saveas(gcf,jpgname,'tif');
                    grid off
                    
                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')
                    
                    close all;
                end
                fig_flag=0;
            end 

% start-------------------- 4con amp validiation scatter -refined
            fig_flag=fig_flags{21,2};
            figname=fig_flags{21,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_4con_amp_val_ref_scatter_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m2_amp=tcon(:,:,tide_info.index(1),1).*100;

                    opts = spreadsheetImportOptions("NumVariables", 12);
                    opts.Sheet = "Sheet1";
                    opts.DataRange = "B2:M20";
                    opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
                    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
                    tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data4_obs_and_khoa.xlsx", opts, "UseExcel", false);
                    clear opts

    %                 if (exist('tcon_model') ~= 1)
                    tconname=[filedir, num2str(tempyear,'%04i'), '\', testname, '_tcon_val_zeta_ref_', num2str(tempyear, '%04i'), '.mat'];
                    if (exist(tconname , 'file') ~= 2)
                        tcon_model=cell(size(tcon_obs));
                        tcon_model(:,1)=table2cell(tcon_obs(:,1));
                        tcon_model(:,2)=table2cell(tcon_obs(:,2));
                        for ind_sta=1:size(tcon_obs,1)
                            tgtlat=table2array(tcon_obs(ind_sta,3));
                            tgtlon=table2array(tcon_obs(ind_sta,4));

                            for i=1:size(lon_rho,1)
                                for j=1:size(lat_rho,2)
                                    if isfinite(tcon(i,j,tide_info.index(1),1))
                                        dist(i,j)=m_lldist([lon_rho(i,j), tgtlon], [lat_rho(i,j), tgtlat]);
                                    else
                                        dist(i,j)=NaN;
                                    end
                                end
                            end

                            ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
                            ind_lon=mod(ind_sta_model(ind_sta),size(lon_rho,1));
                            ind_lat=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
                            tcon_model{ind_sta,4}=ind_lon;
                            tcon_model{ind_sta,3}=ind_lat;
                            tcon_model{ind_sta,5}=tcon(ind_lon,ind_lat,tide_info.index(1),1).*100;
                            tcon_model{ind_sta,6}=tcon(ind_lon,ind_lat,tide_info.index(2),1).*100;
                            tcon_model{ind_sta,7}=tcon(ind_lon,ind_lat,tide_info.index(3),1).*100;
                            tcon_model{ind_sta,8}=tcon(ind_lon,ind_lat,tide_info.index(4),1).*100;
                            tcon_model{ind_sta,9}=tcon(ind_lon,ind_lat,tide_info.index(1),3).*100;
                            tcon_model{ind_sta,10}=tcon(ind_lon,ind_lat,tide_info.index(2),3).*100;
                            tcon_model{ind_sta,11}=tcon(ind_lon,ind_lat,tide_info.index(3),3).*100;
                            tcon_model{ind_sta,12}=tcon(ind_lon,ind_lat,tide_info.index(4),3).*100;
                        end
                        save(tconname, 'tcon_model')
                    else
                        load(tconname)
                    end


                    model_4con= cell2mat(tcon_model(:,5)) + cell2mat(tcon_model(:,6)) + cell2mat(tcon_model(:,7)) + cell2mat(tcon_model(:,8));
                    obs_4con =  table2array(tcon_obs(:,5)) + table2array(tcon_obs(:,6)) + table2array(tcon_obs(:,7)) + table2array(tcon_obs(:,8));

    %                 num_sta=size(tcon_obs,1);
                    mslplot{1}=scatter(obs_4con, model_4con, 'x')
                    hold on
                    maxval=max([obs_4con,  model_4con]);
                    mslplot{2}=plot(0:maxval, 0:maxval);
                    hold off

                    xlabel('4con Observed Tidal amplitude (cm)')
                    ylabel('4con Model Tidal amplitude (cm)')
                    title([testname, ', 4con amp, ', num2str(tempyear,'%04i')])
                    % datetick('x','yyyy','keepticks')
                    axis tight;
                    % ylim(meanplotlev2)
                    set(mslplot{1},'LineWidth',2);
                    set(mslplot{2},'LineWidth',2);

                    set(gca,'FontSize',15);
                    grid on
    %                 lgd=legend('Model','TG-Choi');
    %                 set(lgd,'FontSize',15);
    %                 set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
    %                 set(lgd,'Orientation','horizontal');

                    set(gcf,'PaperPosition', [0 0 36 12]) 
                    % hold off
                    saveas(gcf,jpgname,'tif');
                    grid off
                    
                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')
                    close all;
                end
                fig_flag=0;
            end 
            
 % start-------------------- 4con amp validiation scatter (all test)
            fig_flag=fig_flags{22,2};
            figname=fig_flags{22,1};
            figalldir =strcat('D:\MEPL\project\SSH\5th_year\figure\nwp_1_20\','all','\',regionname,'\Tide\'); % % where figure files will be saved
            if (exist(strcat(figalldir) , 'dir') ~= 7)
                mkdir(strcat(figalldir));
            end 
            while (fig_flag)
                jpgname=strcat(figalldir, 'all','_',regionname, '_', varname, '_4con_amp_val_ref_scatter_', ...
                    num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    
                    
                    opts = spreadsheetImportOptions("NumVariables", 12);
                    opts.Sheet = "Sheet1";
                    opts.DataRange = "B2:M20";
                    opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
                    opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
                    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
                    tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data4_obs_and_khoa.xlsx", opts, "UseExcel", false);
                    clear opts
                    
                    for temp_testind=1:length(all_testname2)
                        temp_testname=all_testname2{temp_testind};
                        temp_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
%                         temp_filename=[temp_filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
%                         temp_tcon=ncread(temp_filename, 'tcon');
%                         all_m2_amp(1:size(temp_tcon,1), 1:size(temp_tcon,2), temp_testind)=temp_tcon(:,:,tide_info.index(1),1).*100;
                        temp_tconname=[temp_filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_tcon_val_zeta_ref_', num2str(tempyear, '%04i'), '.mat'];
                        load(temp_tconname)
                        size1=size(tcon_model,1);
                        size2=size(tcon_model,2);
                        for cellx=1:size1
                            for celly=1:4
                                all_tcon_model(size1*(temp_testind-1)+cellx,celly)=tcon_model{cellx,celly+4};
                                all_tcon_obs(size1*(temp_testind-1)+cellx,celly)=table2array(tcon_obs(cellx,celly+4));
                            end
                        end
                    end
                    
    %                 if (exist('tcon_model') ~= 1)

                    model_4con= sum(all_tcon_model,2);
                    obs_4con =  sum(all_tcon_obs,2);

    %                 num_sta=size(tcon_obs,1);
                    mslplot{1}=scatter(obs_4con(1:cellx), model_4con(1:cellx), 'o');
                    hold on
                    mslplot{2}=scatter(obs_4con(cellx+1:cellx*2), model_4con(cellx+1:cellx*2), '+');
                    mslplot{3}=scatter(obs_4con(cellx*2+1:cellx*3), model_4con(cellx*2+1:cellx*3), '*');
                    mslplot{4}=scatter(obs_4con(cellx*3+1:cellx*4), model_4con(cellx*3+1:cellx*4), 'x');
                    maxval=max([obs_4con,  model_4con]);
                    mslplot{5}=plot(0:maxval, 0:maxval,'r');
                    hold off
                    
                    xlabel('Observed Tidal amplitude  (cm)')
                    ylabel('RCM Tidal amplitude (cm)')
%                     title(['all test', ', 4con amp, ', num2str(tempyear,'%04i')])
                  
                    % datetick('x','yyyy','keepticks')
                    axis tight;
                    % ylim(meanplotlev2)
                    set(mslplot{1},'LineWidth',2);
                    set(mslplot{2},'LineWidth',2);
                    set(mslplot{3},'LineWidth',2);
                    set(mslplot{4},'LineWidth',2);
                    set(mslplot{5},'LineWidth',5);
                    set(mslplot{1},'SizeData',100);
                    set(mslplot{2},'SizeData',100);
                    set(mslplot{3},'SizeData',100);
                    set(mslplot{4},'SizeData',100);
                    
                    set(gca,'FontSize',25);
                    grid on
                    
%                     lgd=legend(all_testname2);
                    lgd=legend({'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'});

                    set(lgd,'FontSize',25);
                    set(lgd,'Position',[0.70 0.20, 0.20, 0.03]);
                    set(lgd,'Orientation','vertical');

                    set(gcf,'PaperPosition', [0 0 24 24]) 
                    % hold off
                    saveas(gcf,jpgname,'tif');
                    grid off
                    
                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')
                    close all;
                end
                fig_flag=0;
            end 
            
% start-------------------- 4con amp diff all test
            fig_flag=fig_flags{23,2};
            figname=fig_flags{23,1};
            while (fig_flag)
                jpgname=strcat(figalldir, all_testname2{1},'_',all_testname2{end},'_',regionname, '_', varname, '_all_4con_amp_diff_', ...
                    num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    
                    for temp_testind=1:length(all_testname2)
                        temp_testname=all_testname2{temp_testind};
                        temp_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
%                         temp_filename=[temp_filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
%                         temp_tcon=ncread(temp_filename, 'tcon');
%                         all_m2_amp(1:size(temp_tcon,1), 1:size(temp_tcon,2), temp_testind)=temp_tcon(:,:,tide_info.index(1),1).*100;
                        temp_diffname=[temp_filedir, temp_testname, '_diff_4con_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
                        clear con4_amp_diff
                        load(temp_diffname)
                        all_con4_amp_diff(:,:,temp_testind)=con4_amp_diff;
                    end
                    mean_con4_amp_diff=mean(all_con4_amp_diff,3);
%                      all_con4_amp_diff(100,100,:)
%                     mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%                     mask_model(mask_model==0)=NaN;    
%                     con4_amp_diff=con4_amp_diff.*mask_model;
                    
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
%                     ax1=axes;
                    plot1=m_pcolor(double(lon_rho)',lat_rho',squeeze(mean_con4_amp_diff'));
%                     plot2=m_pcolor(double(lon_rho)',lat_rho',squeeze(con4_amp_diff'));

%                     m_surface(double(lon_rho)',lat_rho',squeeze(con4_amp_diff'),[0.8 0.8 0.8])

    % %                     [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(m2_amp_diff'), [0:20:240], 'color','k', ...
    % %                         'linewidth', 1.5, 'linestyle', '-');
    % %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    % %                         'Rotation', 0,'fontweight', 'bold');
                    shading(gca,m_pcolor_shading_method);
%                     set(plot1, 'facecolor', 'r');
%                     m_gshhs_i('color',m_gshhs_line_color);
%                     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%                     set(gca1, 'color', [0.8,0.8,0.8]);

%                     surface('Tag','m_grid_color','Parent',Parent1,'ZData',ZData1,'YData',YData1,...
%                         'XData',XData1,...
%                         'LineStyle','none',...
%                         'FaceColor',[0.8 0.8 0.8]);
%                     m_patch([130 140  140 130], [30 30 35 35], 'w', 'LineStyle', 'none')

                    titlename = strcat('4con amp diff, ',testname, ',(',num2str(inputyear(end),'%04i'),'-',num2str(inputyear(1),'%04i'), ')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([-15 15]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);
%                     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');
%                     save([filedir,testname,'_diff_4con_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'], 'con4_amp_diff');
                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end             
            
% start-------------------- k1 amp diff
            fig_flag=fig_flags{24,2};
            figname=fig_flags{24,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_k1_amp_diff_', ...
                    num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    tcon_end=ncread(endfilename, 'tcon');
                    startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    tcon_start=ncread(startfilename, 'tcon');

                    k1_amp_end=tcon_end(:,:,tide_info.index(3),1).*100;
                    k1_amp_start=tcon_start(:,:,tide_info.index(3),1).*100;
                    k1_amp_diff=k1_amp_end-k1_amp_start;

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(k1_amp_diff'));
    % %                     [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(k1_amp_diff'), [0:20:240], 'color','k', ...
    % %                         'linewidth', 1.5, 'linestyle', '-');
    % %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    % %                         'Rotation', 0,'fontweight', 'bold');

                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('k1 amp diff, ',testname, ',(',num2str(inputyear(end),'%04i'),'-',num2str(inputyear(1),'%04i'), ')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([-10 10]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');
                    save([filedir,testname,'_diff_k1_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'], 'k1_amp_diff');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 
            
% start-------------------- o1 amp diff
            fig_flag=fig_flags{25,2};
            figname=fig_flags{25,1};
            while (fig_flag)
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_', varname, '_o1_amp_diff_', ...
                    num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    tcon_end=ncread(endfilename, 'tcon');
                    startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    tcon_start=ncread(startfilename, 'tcon');

                    o1_amp_end=tcon_end(:,:,tide_info.index(4),1).*100;
                    o1_amp_start=tcon_start(:,:,tide_info.index(4),1).*100;
                    o1_amp_diff=o1_amp_end-o1_amp_start;

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(o1_amp_diff'));
    % %                     [C,h2]=m_contour(double(lon_rho)',lat_rho', squeeze(o1_amp_diff'), [0:20:240], 'color','k', ...
    % %                         'linewidth', 1.5, 'linestyle', '-');
    % %                     clabel(C,h2,'FontSize',13,'Color','k', ...
    % %                         'Rotation', 0,'fontweight', 'bold');

                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('o1 amp diff, ',testname, ',(',num2str(inputyear(end),'%04i'),'-',num2str(inputyear(1),'%04i'), ')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([-10 10]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'backcolor', [0.8 0.8 0.8]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');
                    save([filedir,testname,'_diff_o1_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'], 'o1_amp_diff');

                    disp(' ')
                    disp([figname, ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
                fig_flag=0;
            end 
            
            
            
        end
    end
end