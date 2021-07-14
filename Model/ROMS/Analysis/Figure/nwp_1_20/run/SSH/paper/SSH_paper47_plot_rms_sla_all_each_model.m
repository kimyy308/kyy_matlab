close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
all_testname2 = {'test11', 'test12'};
testname = 'test11';
testname2 = 'test12';
sodatestname = 'SODA_3_4_2';
% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

all_region2 ={'NWP'};
% all_region2 ={'AKP4'};
% all_region2 ={'AKP4','YSECS', 'ECS2', 'ES', 'YS', 'NES', 'SES'}



% start-------------------- configuration
for regionind2=1:length(all_region2)
    for folding=1:1

    close all;
    clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 testname testname2 sodatestname
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

%     regionind2 = 2;

    shadlev = [0 35];
    rms_shadlev = [0 4];
    %     trendlev = [-3 3];  %% trend lev
    trendlev = [-10 10];  %% trend lev
    abstrendlev =[2 7];
    reltrendlev =[-5 5];
    conlev  = 0:5:35;
    meanplotlev =[-0.2 0.2];
    trendplotlev = [0 7];
    trenddifflev = [-10 10];
    sshlev =[-0.3 0.3];
    sshdifflev = [0 20];

    % for snu_desktopd
    % testname=all_testname2{testnameind2}    % % need to change
    inputyear = [1994:2014]; % % put year which you want to plot [year year ...]
    inputmonth = [1:12]; % % put month which you want to plot [month month ...]

    varname ='zeta';
    variable='SSH';
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
        case('SES') %% Southern East Sea
            refpolygon=sespolygon;
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        case('ECS2') %% East China Sea
            refpolygon=ecs2polygon;
        case('YSECS') %% East China Sea
            refpolygon=ysecspolygon;
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
    filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

    confilename = ['E:\Data\Model\ROMS\nwp_1_10\',testname2,'\run\',testname2,'_',regionname, ...
                '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

    sodafilename =  ['E:\Data\Reanalysis\SODA\', sodatestname, '\',sodatestname,'_',regionname, ...
                '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];



    valnum=0;
    run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
    wrmap = bwrmap(51:100,:);


    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\all\',testname,'\',regionname,'\'); % % where figure files will be saved
    %             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
        param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_', regionname, '.m']
        filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
    elseif (strcmp(system_name,'GLNXA64'))
    end
    %         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

    run(param_script);

    figdir=[figrawdir,'Trend\SSH\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);

    cmems_trend=ncread(filename, 'cmems_trend');
    cmems_mask=ones(size(cmems_trend));
    cmems_mask(isnan(cmems_trend))=NaN;

    lon_rho=ncread(filename,'lon_rho');
    lat_rho=ncread(filename,'lat_rho');

    end
    % end-------------------- configuration

    lon_rho=ncread(filename,'lon_rho');
    lat_rho=ncread(filename,'lat_rho');

    cmems_sla = ncread(filename, 'cmems_sla');
    cmems_2y_lowpass_sla = ncread(filename, 'cmems_2y_lowpass');
    xlen=size(cmems_sla,1); ylen = size(cmems_sla,2); tlen = size(cmems_sla,3);

    % start-------------------- rms_ lowpassed get
    for folding=1:1
            nyears=[1:5];
            nc_varname_prefixes={'interped_sla'};
            nc_titlename_prefixes={'lowpass'};

    %         for nc_varnameij=1:length(nc_varname_prefixes)
    %             nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
    %             nc_varname=[nc_varname_prefix];
    %             eval(['corr=ncread(filename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr(1,nc_varnameij)=mean(corr(:),'omitnan');
    %             eval(['corr=ncread(confilename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr_con(1,nc_varnameij)=mean(corr(:),'omitnan');
    %             eval(['corr=ncread(sodafilename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr_soda(1,nc_varnameij)=mean(corr(:),'omitnan');
    %         end
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, '_', num2str(nyear),'y_lowpass'];
    %                 jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
    %                     '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    %                 if (exist(jpgname , 'file') ~= 2)



                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname, '''', ',' '''', '\run\','''', ',', '''', testname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_sat(1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname2, '''', ',' '''', '\run\','''', ',', '''', testname2, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_con(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_glo(1,nc_varnameij,:,:,:)=sq_diff_all;

                        eval(['fname=[','''','E:\Data\Reanalysis\SODA\', '''', ',', '''', sodatestname, '''', ',' '''', '\','''', ',', '''', sodatestname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_soda(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_gcm_soda(1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname, '''', ',' '''', '\run\','''', ',', '''', testname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_sat(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname2, '''', ',' '''', '\run\','''', ',', '''', testname2, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_con(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_glo(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Reanalysis\SODA\', '''', ',', '''', sodatestname, '''', ',' '''', '\','''', ',', '''', sodatestname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_soda(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_gcm_soda(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
    %                     eval(['corr=ncread(confilename,','''','rms_', nc_varname,'''',');']);
    %                     corr=corr.*cmems_mask;
    %                     mcorr_con(nyeari+1,nc_varnameij)=mean(corr(:),'omitnan');
    %                     eval(['corr=ncread(sodafilename,','''','rms_', nc_varname,'''',');']);
    %                     corr=corr.*cmems_mask;
    %                     mcorr_soda(nyeari+1,nc_varnameij)=mean(corr(:),'omitnan');

                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

    %                     eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        tit_prefix=nc_titlename_prefixes{1};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mrms,2)));  

    %                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    %                     saveas(gcf,jpgname,'tif');

                        disp(' ')
%                         disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
    %                     disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
    %                 end
                end
            end
    end
% end-------------------- rms_ lowpassed get



% start-------------------- rms_ lowpassed corrected get
    for folding=1:1
            nyears=[1:5];
            nc_varname_prefixes={'corrected_interped_sla'};
            nc_titlename_prefixes={'lowpass'};

            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, '_', num2str(nyear),'y_lowpass'];


                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname, '''', ',' '''', '\run\','''', ',', '''', testname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_corrected(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_sat(1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname2, '''', ',' '''', '\run\','''', ',', '''', testname2, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_con_corrected(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_glo(1,nc_varnameij,:,:,:)=sq_diff_all;

                        eval(['fname=[','''','E:\Data\Reanalysis\SODA\', '''', ',', '''', sodatestname, '''', ',' '''', '\','''', ',', '''', sodatestname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_soda_corrected(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_gcm_soda(1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname, '''', ',' '''', '\run\','''', ',', '''', testname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_corrected(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_sat(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname2, '''', ',' '''', '\run\','''', ',', '''', testname2, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_con_corrected(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_glo(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Reanalysis\SODA\', '''', ',', '''', sodatestname, '''', ',' '''', '\','''', ',', '''', sodatestname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_soda_corrected(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_gcm_soda(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;

                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

                        tit_prefix=nc_titlename_prefixes{1};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mrms,2)));  

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 


                        disp(' ')
                        disp(' ')
                        disp(' ')

                        hold off
                        close all;
    %                 end
                end
            end
    end
% end-------------------- rms_ lowpassed corrected get


% % % start-------------------- RMS perfomance by lowpass filter window
 % (RCM-GLO)
    for folding=1:1
         jpgname=strcat(outfile, '_', testname, '_',regionname, '_rcm_glo_rms_lp_perfomance_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%             if (exist(jpgname , 'file') ~= 2)
            mslplot=plot(0:length(nyears),mrms_soda(:,1),'b')
            hold on
            msl2plot=plot(0:length(nyears),mrms_soda_corrected(:,1),'b')
            msl3plot=plot(0:length(nyears),mrms_con(:,1),'g')
            msl4plot=plot(0:length(nyears),mrms_con_corrected(:,1),'g')

            xlabel('Low pass filter window(year)')
            ylabel('RMS')
            title([regionname, ', RMS performance '])
    %         datetick('x','yyyy','keepticks')
            axis tight;
            ylim([1.5 6.5])
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(msl4plot,'LineWidth',2);
            set(msl2plot,'LineStyle','--');
            set(msl4plot,'LineStyle','--');

            set(gca,'FontSize',15);
            grid on
            lgd=legend('GCM-SODA','GCM-SODA-corrected','RCM-GLO','RCM-GLO-corrected');
            set(lgd,'FontSize',10);
            % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(gcf,'PaperPosition', [0 0 24 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            set(gca,'XTick', (0:1:length(nyears)))
            saveas(gcf,jpgname,'jpg');
            grid off
%             end
        close all;
    end
    % % % end-------------------- RMS perfomance by lowpass filter window
    % (RCM-GLO)
    
 
% % % start-------------------- RMS perfomance by lowpass filter window
% (RCM-SAT)
    for folding=1:1
         jpgname=strcat(outfile, '_', testname, '_',regionname, '_rcm_sat_rms_lp_perfomance_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%             if (exist(jpgname , 'file') ~= 2)
            mslplot=plot(0:length(nyears),mrms_con(:,1),'g')
            hold on
            msl2plot=plot(0:length(nyears),mrms_con_corrected(:,1),'g')
            msl3plot=plot(0:length(nyears),mrms(:,1),'r')
%             msl4plot=plot(0:length(nyears),mrms_corrected(:,1),'y')

            xlabel('Low pass filter window(year)')
            ylabel('RMS')
            title([regionname, ', RMS performance '])
    %         datetick('x','yyyy','keepticks')
            axis tight;
            ylim([1.5 6.5])
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(msl2plot,'LineStyle','--');

%             set(msl4plot,'LineWidth',2);

            set(gca,'FontSize',15);
            grid on
            lgd=legend('RCM-SAT','RCM-GLO','RCM-GLO-corrected');
%             lgd=legend('RCM-GLO','RCM-SAT','RCM-GLO-corrected', 'RCM-SAT-corrected');
            set(lgd,'FontSize',10);
            % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(gcf,'PaperPosition', [0 0 24 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            set(gca,'XTick', (0:1:length(nyears)))
            saveas(gcf,jpgname,'jpg');
            grid off
%             end
        close all;
    end
    % % % end-------------------- RMS perfomance by lowpass filter window
% (RCM-SAT)
end

