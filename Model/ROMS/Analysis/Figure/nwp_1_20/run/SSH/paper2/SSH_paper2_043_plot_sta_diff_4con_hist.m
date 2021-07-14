close all; clear all; clc;

dropboxpath='C:\Users\User\Dropbox';
% addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
% addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\t_tide_v1.3beta']));
% addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
% addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
% addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

pastyear=1989;
futureyear=2005;


% % set info of station
for folding=1:1
station.name{1}='Anheung';
station.name{2}='Gunsan';
station.name{3}='Mokpo';
station.name{4}='Heuksando';
station.name{5}='Chujado';
station.name{6}='Wando';
station.name{7}='Yeosu';
station.name{8}='Tongyeong';
station.name{9}='Gadeokdo';
station.name{10}='Busan';
station.name{11}='Ieodo';
station.name{12}='Jeju';
station.name{13}='Seogwipo';
station.name{14}='Geomundo';
station.name{15}='Ulsan';
station.name{16}='Pohang';
station.name{17}='Mukho';
station.name{18}='Sokcho';
station.name{19}='Ulleungdo';

num_sta=length(station.name);

station.name{1}='Anheung';
station.name{2}='Gunsan';
station.name{3}='Mokpo';
station.name{4}='Heuksando';
station.name{5}='Chujado';
station.name{6}='Wando';
station.name{7}='Yeosu';
station.name{8}='Tongyeong';
station.name{9}='Gadeokdo';
station.name{10}='Busan';
station.name{11}='Ieodo';
station.name{12}='Jeju';
station.name{13}='Seogwipo';
station.name{14}='Geomundo';
station.name{15}='Ulsan';
station.name{16}='Pohang';
station.name{17}='Mukho';
station.name{18}='Sokcho';
station.name{19}='Ulleungdo';

load('Z:\내 드라이브\MEPL\project\SSH\4th_year\09_backup\UST\99_재생산자료\01_매시별\tide_all.mat')

Obs.lat(1)=36.67;
Obs.lat(2)=35.97;
Obs.lat(3)=34.77;
Obs.lat(4)=34.68;
Obs.lat(5)=33.95;
Obs.lat(6)=34.30;
Obs.lat(7)=34.73;
Obs.lat(8)=34.82;
Obs.lat(9)=35.02;
Obs.lat(10)=35.08;
Obs.lat(11)=NaN;
Obs.lat(12)=33.52;
Obs.lat(13)=33.23;
Obs.lat(14)=34.02;
Obs.lat(15)=35.5;
Obs.lat(16)=36.03;
Obs.lat(17)=37.55;
Obs.lat(18)=38.2;
Obs.lat(19)=37.48;

end
% % get sum of 4 major tidal constitute (Observation)

% % get py data
for stai=1:num_sta
    if (stai~=11)
        data_1y_ind=find(Obs.data{stai}(:,1)==pastyear);
        Obs_py.data_1y{stai}=Obs.data{stai}(data_1y_ind,6);
        
        [tname,tfreq,tcon,tout]=t_tide(Obs_py.data_1y{stai},...
               'interval',1, ...                     % hourly data
               'start',datenum(pastyear,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
               'latitude',Obs.lat(stai),...               % Latitude of Obs_py
               'rayleigh',1, 'error','wboot');                       % Use SNR=1 for synthesis. 

        Obs_py.tcon(stai,:,:)=tcon;

        tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio
        
        tide_info.name{1}='M2  ';
        tide_info.name{2}='S2  ';
        tide_info.name{3}='K1  ';
        tide_info.name{4}='O1  ';
        
        num_tide_all=size(tname,1);
        num_tide_tgt=length(tide_info.name);
        for coni=1:num_tide_all
            for tide_namei=1:num_tide_tgt
                if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                    tide_info.index(tide_namei)=coni;
                end
            end
        end

        Obs_py.amp_4con(stai)= sum(tcon(tide_info.index(:),1));
        Obs_py.amp_M2(stai)=tcon(tide_info.index(1),1);
        Obs_py.amp_S2(stai)=tcon(tide_info.index(2),1);
        Obs_py.amp_K1(stai)=tcon(tide_info.index(3),1);
        Obs_py.amp_O1(stai)=tcon(tide_info.index(4),1);

    else
        Obs_py.data_1y{stai}=NaN;
        Obs_py.amp_4con(stai)= NaN;
        Obs_py.amp_M2(stai)=NaN;
        Obs_py.amp_S2(stai)=NaN;
        Obs_py.amp_K1(stai)=NaN;
        Obs_py.amp_O1(stai)=NaN;
    end  
end

% %  get future year data
for stai=1:num_sta
    if (stai~=11)
        data_1y_ind=find(Obs.data{stai}(:,1)==futureyear);
        Obs_fy.data_1y{stai}=Obs.data{stai}(data_1y_ind,6);
        
        [tname,tfreq,tcon,tout]=t_tide(Obs_fy.data_1y{stai},...
               'interval',1, ...                     % hourly data
               'start',datenum(futureyear,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
               'latitude',Obs.lat(stai),...               % Latitude of Obs_fy
               'rayleigh',1, 'error','wboot');                       % Use SNR=1 for synthesis. 

        Obs_fy.tcon(stai,:,:)=tcon;

        tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio
        
        tide_info.name{1}='M2  ';
        tide_info.name{2}='S2  ';
        tide_info.name{3}='K1  ';
        tide_info.name{4}='O1  ';
        
        num_tide_all=size(tname,1);
        num_tide_tgt=length(tide_info.name);
        for coni=1:num_tide_all
            for tide_namei=1:num_tide_tgt
                if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                    tide_info.index(tide_namei)=coni;
                end
            end
        end

        Obs_fy.amp_4con(stai)= sum(tcon(tide_info.index(:),1));
        Obs_fy.amp_M2(stai)=tcon(tide_info.index(1),1);
        Obs_fy.amp_S2(stai)=tcon(tide_info.index(2),1);
        Obs_fy.amp_K1(stai)=tcon(tide_info.index(3),1);
        Obs_fy.amp_O1(stai)=tcon(tide_info.index(4),1);

    else
        Obs_fy.data_1y{stai}=NaN;
        Obs_fy.amp_4con(stai)= NaN;
        Obs_fy.amp_M2(stai)=NaN;
        Obs_fy.amp_S2(stai)=NaN;
        Obs_fy.amp_K1(stai)=NaN;
        Obs_fy.amp_O1(stai)=NaN;
    end  
end

all_testname2 = {'test53', 'test54', 'test55', 'test56'};

for testnameind2=1:length(all_testname2)
    
    % % %         flag configuration
    for folding=1:1
        fig_flags{1,1}='M2 diff';
        fig_flags{1,1}='S2 diff';
        fig_flags{1,1}='K1 diff';
        fig_flags{1,1}='O1 diff';
    end
    for flagi=1:8
        fig_flags{flagi,2}=0;
    end

    fig_flags{1,2}=1;
    fig_flags{2,2}=1;
    fig_flags{3,2}=1;
    fig_flags{4,2}=1;

    testname=all_testname2{testnameind2};
% % get py 
    filename=['J:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(pastyear),'\sta.nc'];
    Model_py.zeta=ncread(filename, 'zeta')*100;
    Model_py.lon=ncread(filename, 'lon_rho');
    Model_py.lat=ncread(filename, 'lat_rho');
    Model_py.ocean_time=ncread(filename, 'ocean_time');
    
    for stai=1:num_sta    
        [tname,tfreq,tcon,tout]=t_tide(Model_py.zeta(stai,:),...
               'interval',1, ...                     % hourly data
               'start',datenum(pastyear,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
               'latitude',Model_py.lat(stai),...               % Latitude of Model_py
               'rayleigh',1, 'error','wboot');                       % Use SNR=1 for synthesis. 

           Model_py.tcon(stai,:,:)=tcon;

        tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio

        num_tide_all=size(tname,1);
        num_tide_tgt=length(tide_info.name);
        for coni=1:num_tide_all
            for tide_namei=1:num_tide_tgt
                if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                    tide_info.index(tide_namei)=coni;
                end
            end
        end

        Model_py.amp_4con(stai)= sum(tcon(tide_info.index(:),1));
        Model_py.amp_M2(stai)=tcon(tide_info.index(1),1);
        Model_py.amp_S2(stai)=tcon(tide_info.index(2),1);
        Model_py.amp_K1(stai)=tcon(tide_info.index(3),1);
        Model_py.amp_O1(stai)=tcon(tide_info.index(4),1);
    end
    
% % get fy 
    filename=['J:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(futureyear),'\sta.nc'];
    Model_fy.zeta=ncread(filename, 'zeta')*100;
    Model_fy.lon=ncread(filename, 'lon_rho');
    Model_fy.lat=ncread(filename, 'lat_rho');
    Model_fy.ocean_time=ncread(filename, 'ocean_time');
    
    for stai=1:num_sta    
        [tname,tfreq,tcon,tout]=t_tide(Model_fy.zeta(stai,:),...
               'interval',1, ...                     % hourly data
               'start',datenum(futureyear,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
               'latitude',Model_fy.lat(stai),...               % Latitude of Model_fy
               'rayleigh',1, 'error','wboot');                       % Use SNR=1 for synthesis. 

           Model_fy.tcon(stai,:,:)=tcon;

        tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio

        num_tide_all=size(tname,1);
        num_tide_tgt=length(tide_info.name);
        for coni=1:num_tide_all
            for tide_namei=1:num_tide_tgt
                if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                    tide_info.index(tide_namei)=coni;
                end
            end
        end

        Model_fy.amp_4con(stai)= sum(tcon(tide_info.index(:),1));
        Model_fy.amp_M2(stai)=tcon(tide_info.index(1),1);
        Model_fy.amp_S2(stai)=tcon(tide_info.index(2),1);
        Model_fy.amp_K1(stai)=tcon(tide_info.index(3),1);
        Model_fy.amp_O1(stai)=tcon(tide_info.index(4),1);
    end
    
    
    Obs.M2_diff=Obs_fy.amp_M2-Obs_py.amp_M2;
    Obs.S2_diff=Obs_fy.amp_S2-Obs_py.amp_S2;
    Obs.K1_diff=Obs_fy.amp_K1-Obs_py.amp_K1;
    Obs.O1_diff=Obs_fy.amp_O1-Obs_py.amp_O1;
    Model.M2_diff=Model_fy.amp_M2-Model_py.amp_M2;
    Model.S2_diff=Model_fy.amp_S2-Model_py.amp_S2;
    Model.K1_diff=Model_fy.amp_K1-Model_py.amp_K1;
    Model.O1_diff=Model_fy.amp_O1-Model_py.amp_O1;

    figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\station\'); % % where figure files will be saved

    figdir=[figrawdir,'CLIM\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir);
        
    fig_flag=fig_flags{1,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_M2_diff_',num2str(pastyear), '-', num2str(futureyear), '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)

            mslplot{1}=plot(Model.M2_diff, 'b');
            hold on
            mslplot{2}=plot(Obs.M2_diff, 'k');
            hold off
            xticks(1:num_sta)
            xticklabels(station.name)
            xtickangle(45)

            xlabel('Station')
            ylabel('Tidal amplitude difference (cm)')
            title([num2str(min(pastyear),'%04i'),'-',num2str(max(futureyear),'%04i')])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);


            set(gca,'FontSize',15);
            grid on
            lgd=legend('Model','TG-UST');
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
            saveas(gcf,pngname,'tif');
            grid off
            close all;
        end
        fig_flag=0;
    end
     
    fig_flag=fig_flags{2,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_S2_diff_',num2str(pastyear), '-', num2str(futureyear), '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)

            mslplot{1}=plot(Model.S2_diff, 'b');
            hold on
            mslplot{2}=plot(Obs.S2_diff, 'k');
            hold off
            xticks(1:num_sta)
            xticklabels(station.name)
            xtickangle(45)

            xlabel('Station')
            ylabel('Tidal amplitude difference (cm)')
            title([num2str(min(pastyear),'%04i'),'-',num2str(max(futureyear),'%04i')])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);


            set(gca,'FontSize',15);
            grid on
            lgd=legend('Model','TG-UST');
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
            saveas(gcf,pngname,'tif');
            grid off
            close all;
        end
        fig_flag=0;
    end
     
    fig_flag=fig_flags{3,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_K1_diff_',num2str(pastyear), '-', num2str(futureyear), '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)

            mslplot{1}=plot(Model.K1_diff, 'b');
            hold on
            mslplot{2}=plot(Obs.K1_diff, 'k');
            hold off
            xticks(1:num_sta)
            xticklabels(station.name)
            xtickangle(45)

            xlabel('Station')
            ylabel('Tidal amplitude difference (cm)')
            title([num2str(min(pastyear),'%04i'),'-',num2str(max(futureyear),'%04i')])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);


            set(gca,'FontSize',15);
            grid on
            lgd=legend('Model','TG-UST');
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
            saveas(gcf,pngname,'tif');
            grid off
            close all;
        end
        fig_flag=0;
    end
     
    fig_flag=fig_flags{4,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_O1_diff_',num2str(pastyear), '-', num2str(futureyear), '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)

            mslplot{1}=plot(Model.O1_diff, 'b');
            hold on
            mslplot{2}=plot(Obs.O1_diff, 'k');
            hold off
            xticks(1:num_sta)
            xticklabels(station.name)
            xtickangle(45)

            xlabel('Station')
            ylabel('Tidal amplitude difference (cm)')
            title([num2str(min(pastyear),'%04i'),'-',num2str(max(futureyear),'%04i')])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);


            set(gca,'FontSize',15);
            grid on
            lgd=legend('Model','TG-UST');
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
            saveas(gcf,pngname,'tif');
            grid off
            close all;
        end
        fig_flag=0;
     end

end