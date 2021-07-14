close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test53', 'test54','test55','ens06','test56'};
% all_testname2 = {'GODAS', 'MyOcean'};
all_testname2 = {'ORAS5'};

% all_region2 ={'NWP', 'YS', 'AKP2'}
% all_region2 ={'AKP4'}

all_region2 ={'AKP4'};

% all_var2 = {'SST', 'SSS', 'SSH', 'BT'};
% all_var2 = {'SST', 'SSH', 'SSS'};

all_var2 = {'SSH'};

% all_region2 ={'NWP'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2
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
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end

%         shadlev = [0 35];
%         rms_shadlev = [0 4];
%     %     trendlev = [-3 3];  %% trend lev
%         trendlev = [-10 10];  %% trend lev
%         abstrendlev =[4 7];
%         reltrendlev =[-5 5];
%         conlev  = 0:5:35;
%         meanplotlev =[-0.3 0.3];
%         trendplotlev = [3 7];
%         sshlev =[-0.7 1.3];
%         sshdifflev = [40 70];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear1 = [1994:2017]; % % should put year from 1994 to 2017
        inputyear2 = [1994:2017]; % % should put year from 1994 to 2017
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]

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
            case('NES') %% Northern East Sea
                refpolygon=nespolygon;
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

        % % % for EKB
        % regionname='EKB';
        % lonlat = [127, 129.5, 38, 40.5];

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
%         filename = ['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_soda_', regionname, '.m']
            filedir = strcat('E:\Data\Reanalysis\',testname,'\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
        variable='SSH';
        run(param_script);

        figdir=[figrawdir,'CLIM\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);

% % start-------------------- later decadal current plot
%         if (length(inputmonth)==1)
%             pngname=strcat(outfile, '_', testname,'_',regionname, '_',calendarname{inputmonth}(1:3),'_uv_',num2str(min(inputyear2),'%04i'), ...
%                 '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
%         else
%             pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_uv_',num2str(min(inputyear2),'%04i'), ...
%                 '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
%         end
%         if (exist(pngname , 'file') ~= 2)        
%             for yearij=1:length(inputyear2)
%                 tempyear=inputyear2(yearij);
%                 yearstr=num2str(tempyear, '%04i');
%                 for monthij=1:length(inputmonth)
%                     tempmonth=inputmonth(monthij);
%                     monthstr=num2str(tempmonth, '%02i');
%                     ufilename=[filedir, 'NWP3soda_u_trend_1994_2017.nc'];
%                     vfilename=[filedir, 'NWP3soda_v_trend_1994_2017.nc'];
%                     
%                     if (exist('lon_rho' , 'var') ~= 1)
%                         lon= ncread(ufilename, 'lon');
%                         lat= ncread(ufilename, 'lat');
%                         [lat_rho, lon_rho]=meshgrid(lat,lon);
%                         [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
%                         cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                         cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                     end
%                     data_info = ncinfo(ufilename, 'soda_u'); 
%                     indij=(tempyear-1994)*12+tempmonth;
%                     u = ncread(ufilename,'soda_u',[lon_min(1) lat_min(1) indij], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%                     v = ncread(vfilename,'soda_v',[lon_min(1) lat_min(1) indij], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
% 
%                     if (exist('mean_u' , 'var') ~= 1)
%                         mean_u=zeros(size(u));
%                         mean_v=zeros(size(v));
%                     end
%                     mean_u=mean_u + (u / (length(inputyear2) * length(inputmonth)));
%                     mean_v=mean_v + (v / (length(inputyear2) * length(inputmonth)));
%                 end
%             end
% %             u_rho = u2rho_2d(mean_u')';
% %             v_rho = v2rho_2d(mean_v')';
%               u_rho = mean_u;
%               v_rho = mean_v;
%               
%             mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%             mask_model(mask_model==0)=NaN;
%             u_rho= u_rho .* mask_model;
%             v_rho= v_rho .* mask_model;
%             
%             if (exist('ref_vec_x_range' , 'var') ~= 1)
%                 ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
%                 ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
% %                 ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
% %                 ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
%                 ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2));
%                 ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));
%             end
%             u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
%             v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;     
% 
%             m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             hold on;
%             m_gshhs_i('color',m_gshhs_line_color)  
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% 
%             uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
%                             cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
%                             u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
%                             v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
%                             'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
% 
%             m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
%             if (length(inputmonth)==1) 
%                 titlename = strcat('UV mean, ',calendarname{inputmonth}(1:3),', ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution                                
%             else
%                 titlename = strcat('UV mean, ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution                
%             end
%             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%             saveas(gcf,pngname,'tif');
%             close all;
%             clear lon_rho mean_u ref_vec_x_range
%         end
% % end-------------------- later decadal current plot

% start-------------------- earlier decadal SST, SSS plot
        for varind2=1:length(all_var2)
            variable=all_var2{varind2};
            if (length(inputmonth)==1)
                pngname=strcat(outfile, '_', testname,'_',regionname, '_',calendarname{inputmonth}(1:3),'_', variable,'_',num2str(min(inputyear1),'%04i'), ...                
                    '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
            else
                pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', variable,'_',num2str(min(inputyear1),'%04i'), ...                
                    '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
            end
            if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                for yearij=1:length(inputyear1)
                    tempyear=inputyear1(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        filename=[filedir, 'NWP',testname,'_',varname,'_trend_1994_2017.nc'];
                        soda_varname=[testname,'_',varname];
%                         if (exist('lon_rho' , 'var') ~= 1)
%                             lon_rho=ncread(filename, 'lon_rho');
%                             lat_rho=ncread(filename, 'lat_rho');
%                             [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
%                             cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                             cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                         end
                        if (exist('lon_rho' , 'var') ~= 1)
                            lon= ncread(filename, 'lon');
                            lat= ncread(filename, 'lat');
                            [lat_rho, lon_rho]=meshgrid(lat,lon);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        data_info = ncinfo(filename, soda_varname); 
                        
                        indij=(tempyear-1994)*12+tempmonth;
                        if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
                            data = ncread(filename,soda_varname,[lon_min(1) lat_min(1) indij], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        elseif (strcmp(variable,'SSH')==1)
                            data = ncread(filename,soda_varname,[lon_min(1) lat_min(1) indij], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                        elseif (strcmp(variable,'BT')==1)
                            data = ncread(filename,soda_varname,[lon_min(1) lat_min(1) indij], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                        end
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(inputyear1) * length(inputmonth)));
                    end
                end

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data= mean_data .* mask_model;

%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                if (strcmp(variable,'SSH')==1)
                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                else
                    m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                end
                shading(gca,m_pcolor_shading_method);   
%                 median(mean_data(:),'omitnan')
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                if (length(inputmonth)==1) 
                    titlename = strcat(variable, ' mean, ',calendarname{inputmonth}(1:3),', ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution                                        
                else
                    titlename = strcat(variable, ' mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution                    
                end
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
                caxis(colorbar_lev);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
        end
% end-------------------- earlier decadal SST SSS plot




    end
end