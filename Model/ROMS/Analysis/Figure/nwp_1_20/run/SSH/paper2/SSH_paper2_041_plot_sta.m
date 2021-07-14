close all; clear all; clc;

dropboxpath='C:\Users\User\Dropbox';
% addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
% addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\t_tide_v1.3beta']));
% addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
% addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
% addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path




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

Obs.amp_4con_choi(1)=NaN;
Obs.amp_4con_choi(2)=209.2+74.8+33.2+24.6;
Obs.amp_4con_choi(3)=122.6+36.7+31.6+23.9;
Obs.amp_4con_choi(4)=NaN;
Obs.amp_4con_choi(5)=NaN;
Obs.amp_4con_choi(6)=108+48+28+21;
Obs.amp_4con_choi(7)=102+47+21+13;
Obs.amp_4con_choi(8)=80+37+16+12;
Obs.amp_4con_choi(9)=56+28+8+4;
Obs.amp_4con_choi(10)=39.8+18.8+4.4+1.5;
Obs.amp_4con_choi(11)=NaN;
Obs.amp_4con_choi(12)=71.3+30.1+23.0+17.3;
Obs.amp_4con_choi(13)=77+34+24+18;
Obs.amp_4con_choi(14)=92+42+24+18;
Obs.amp_4con_choi(15)=17.8+9.4+4.3+3.2;
Obs.amp_4con_choi(16)=NaN;
Obs.amp_4con_choi(17)=6.9+2.3+5.6+4.7;
Obs.amp_4con_choi(18)=7.1+2.6+5.2+4.7;
Obs.amp_4con_choi(19)=4+2+4+4;

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
for stai=1:num_sta
    if (stai~=11)
        data_1y_ind=find(Obs.data{stai}(:,1)==2005);
        Obs.data_1y{stai}=Obs.data{stai}(data_1y_ind,6);
        
        [tname,tfreq,tcon,tout]=t_tide(Obs.data_1y{stai},...
               'interval',1, ...                     % hourly data
               'start',datenum(2005,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
               'latitude',Obs.lat(stai),...               % Latitude of obs
               'rayleigh',1, 'error','wboot');                       % Use SNR=1 for synthesis. 

        Obs.tcon(stai,:,:)=tcon;

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

        Obs.amp_4con(stai)= sum(tcon(tide_info.index(:),1));
    else
        Obs.data_1y{stai}=NaN;
        Obs.amp_4con(stai)= NaN;
    end  
end

% % Obs.tcon(18,tide_info.index(1:4),1) % Ulsan amp

all_testname2 = {'test53', 'test54', 'test55', 'test56'};

for testnameind2=1:length(all_testname2)
    
    % % %         flag configuration
    for folding=1:1
        fig_flags{1,1}='sum of 4 major con';
    end
    for flagi=1:8
        fig_flags{flagi,2}=0;
    end

    fig_flags{1,2}=2;

    testname=all_testname2{testnameind2};
    filename=['J:\Data\Model\ROMS\nwp_1_20\', testname, '\run\1989\sta.nc'];
    Model.zeta=ncread(filename, 'zeta')*100;
    Model.lon=ncread(filename, 'lon_rho');
    Model.lat=ncread(filename, 'lat_rho');
    Model.ocean_time=ncread(filename, 'ocean_time');
    
    for stai=1:num_sta    
        [tname,tfreq,tcon,tout]=t_tide(Model.zeta(stai,:),...
               'interval',1, ...                     % hourly data
               'start',datenum(2006,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
               'latitude',Model.lat(stai),...               % Latitude of obs
               'rayleigh',1, 'error','wboot');                       % Use SNR=1 for synthesis. 

           Model.tcon(stai,:,:)=tcon;
        %    
        % [tname,tfreq,tcon,tout] =t_tide(anassh(5289:14048), ...
        %                             'interval', 1, ...
        %                             'start time', [2003,1,1,0,0,0], ...
        %                             'latitude', lat, ...
        %                             'rayleigh', 1);


        tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio
        % as long as the SNR > 10; and is
        % probably not bad for SNR as low as 2 or 3. The nonlinear
        % procedure gives similar results to the linearized
        % procedure at high SNR, and is more accurate at low
        % SNR.        

        num_tide_all=size(tname,1);
        num_tide_tgt=length(tide_info.name);
        for coni=1:num_tide_all
            for tide_namei=1:num_tide_tgt
                if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                    tide_info.index(tide_namei)=coni;
                end
            end
        end

        Model.amp_4con(stai)= sum(tcon(tide_info.index(:),1));
    end

%     if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\station\'); % % where figure files will be saved
%         param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', 'AKP4', '.m']
%         if (strcmp(scenname,'rcp26')==1)
%             drivename='H';
%         elseif (strcmp(scenname,'rcp85')==1)
%             drivename='G';
%         end
%         filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
%     elseif (strcmp(system_name,'GLNXA64'))
%     end

%     run(param_script);

    figdir=[figrawdir,'CLIM\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir);
        
    fig_flag=fig_flags{1,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_sum_4con_1989', '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)

            mslplot{1}=plot(Model.amp_4con, 'b')
            hold on
            mslplot{2}=plot(Obs.amp_4con_choi, 'r')
            mslplot{3}=plot(Obs.amp_4con, 'k')
            hold off
            xticks(1:num_sta)
            xticklabels(station.name)
            xtickangle(45)

            xlabel('Station')
            ylabel('Tidal amplitude[4 major sum] (cm)')
            % title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);
            set(mslplot{3},'LineWidth',2);


            set(gca,'FontSize',15);
            grid on
            lgd=legend('Model','Report-Choi','TG-UST');
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