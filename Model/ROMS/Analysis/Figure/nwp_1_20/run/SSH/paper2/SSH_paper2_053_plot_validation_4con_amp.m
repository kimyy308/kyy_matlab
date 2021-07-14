close all; clear all; clc;

dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\t_tide_v1.3beta']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path


all_region2 ={'AKP4'};

all_testname2 = {'test66'};
inputyear=[2006, 2100];
% inputyear=[2006];

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
                fig_flags{5,1}='amplitude difference(M2) (end-start)';
                fig_flags{6,1}='phase difference(M2) (end-start)';
            end
            for flagi=1:13
                fig_flags{flagi,2}=0;
            end

            fig_flags{1,2}=0;
            fig_flags{2,2}=1;
            fig_flags{3,2}=1;
            fig_flags{4,2}=2;
            fig_flags{5,2}=1;
            fig_flags{6,2}=1;
            fig_flags{7,2}=1;
            fig_flags{8,2}=0;
            fig_flags{9,2}=0;
            fig_flags{10,2}=0;
            fig_flags{11,2}=0;
            fig_flags{12,2}=0;
            fig_flags{13,2}=0;
            
            testname=all_testname2{testnameind2};
            tempyear=inputyear(yearind2);
            outputdir=['G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear)];
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

                lonfilename=['G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear),'\ocean_his_lon_rho_AKP4', '.nc'];
                latfilename=['G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear),'\ocean_his_lat_rho_AKP4', '.nc'];

                if (exist('lon_rho' , 'var') ~= 1)
                    lon_rho=ncread(lonfilename, 'lon_rho');
                    lat_rho=ncread(latfilename, 'lat_rho');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                    cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
            end

            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
            filedir = strcat('G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            
            run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
            wrmap = bwrmap(51:100,:);
            
            run(param_script);
        
% % % %             ncoutfilename=[outputdir,'\harmonic_analysis_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
% % % % 
% % % %             tidevarid=netcdf.defVar(ncid, 'tname', 'NC_CHAR', [tcon_dimid tide_dimid ]);
% % % %             netcdf.putAtt(ncid,tidevarid,'long_name','tide_name');
% % % % 
% % % %             lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
% % % %             netcdf.putAtt(ncid,lon_rhovarid,'long_name','longitude');
% % % %             netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');
% % % % 
% % % %             lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
% % % %             netcdf.putAtt(ncid,lat_rhovarid,'long_name','latitude');
% % % %             netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
% % % % 
% % % %             m2_ampvarid=netcdf.defVar(ncid, 'm2_amp', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
% % % %             netcdf.putAtt(ncid,m2_ampvarid,'long_name','m2_amplitude');
% % % %             netcdf.putAtt(ncid,m2_ampvarid,'units','m');
% % % % 
% % % %             tconvarid=netcdf.defVar(ncid, 'tcon', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid tide_dimid tcon_dimid]);
% % % %             netcdf.putAtt(ncid,tconvarid,'long_name','tcon(amp, amp_err, phase, phase_err)');
            
            filename=[filedir, num2str(tempyear,'%04i'), '\', 'harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
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
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

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
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

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
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

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
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

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
                    endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', 'harmonic_analysis_', varname, '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    tcon_end=ncread(endfilename, 'tcon');
                    startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', varname, '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
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
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([-10 10]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

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
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

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
                    endfilename=[filedir, num2str(max(inputyear),'%04i'), '\', 'harmonic_analysis_', varname, '_', regionname, '_', num2str(max(inputyear), '%04i'), '.nc'];
                    tcon_end=ncread(endfilename, 'tcon');
                    startfilename=[filedir, num2str(min(inputyear),'%04i'), '\', 'harmonic_analysis_', varname, '_', regionname, '_', num2str(min(inputyear), '%04i'), '.nc'];
                    tcon_start=ncread(startfilename, 'tcon');
                    
                    m2_amp_end=tcon_end(:,:,tide_info.index(2),1).*100;
                    m2_amp_start=tcon_start(:,:,tide_info.index(2),1).*100;
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
                    titlename = strcat('S2 amp diff, ',testname, ',(',num2str(inputyear(end),'%04i'),'-',num2str(inputyear(1),'%04i'), ')');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'cm','fontsize',colorbar_title_fontsize);
                    caxis([-10 10]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

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
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
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
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
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
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
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
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
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
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
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
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
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
        
        
        
        
        end  
    end
end