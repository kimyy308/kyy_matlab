clc; close all; clear all;

RCM_info.years = [2015];
% RCM_info.years = [2015, 2020, 2025, 2030];

% tmp.variable = 'wstrcurl';
% % RCM_info.name = {'test2127', 'test2201'};
RCM_info.name = {'test2201'};
% RCM_info.name = {'test2127'};

dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '/source/matlab/Common/t_tide_v1.3beta']));

tide_info.name{1}='M2  ';
tide_info.name{2}='S2  ';
tide_info.name{3}='K1  ';
tide_info.name{4}='O1  ';
% tide_info.name{5}='MU2 ';
  
for testnameind2=1:length(RCM_info.name)
    for yeari=1:length(RCM_info.years)
        tmp.tempyear=RCM_info.years(yeari);
        tmp.yearstr=num2str(tmp.tempyear);
        tmp.testname=RCM_info.name{testnameind2};   % % need to change
%         dirs.stadir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\',tmp.yearstr, filesep);
%         dirs.stadir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\',tmp.yearstr, '_test_05', filesep);
        dirs.stadir = strcat('D:\MEPL\project\SSH\7th_year(2022)\tide_byun\test_05', filesep);
        tmp.staname = [dirs.stadir, 'sta.nc'];
        
%         ncinfo(tmp.staname)
%         mean_data=ncread(tmp.staname, 'zeta');
        mean_data=ncread(tmp.staname, 'zeta', [1, 1], [inf, 745]);

        lon=ncread(tmp.staname,'lon_rho');
        lat=ncread(tmp.staname,'lat_rho');
        datas.([tmp.testname, '_', tmp.yearstr])=mean_data(1,:); % Incheon -> 1st station
        
        [Model_harmonic.([tmp.testname, '_tname']), ...
            Model_harmonic.([tmp.testname, '_tfreq']), ...
            Model_harmonic.([tmp.testname, '_tcon']), ...
            Model_harmonic.([tmp.testname, '_tout'])]=t_tide_noprint(datas.([tmp.testname, '_', tmp.yearstr]),...
                               'interval',1, ...                     % hourly data
                               'start',datenum(tmp.tempyear,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
                               'latitude',lat(1),...               % Latitude of Model
                               'rayleigh',1, 'error','wboot', 'output', 'none');

        Model_harmonic.([tmp.testname, '_tsnr'])=(Model_harmonic.([tmp.testname, '_tcon'])(:,1)./Model_harmonic.([tmp.testname, '_tcon'])(:,2)).^2;  % signal to noise ratio

        num_tide_all=size(Model_harmonic.([tmp.testname, '_tname']),1);
        num_tide_tgt=length(tide_info.name);
        for coni=1:num_tide_all
            for tide_namei=1:num_tide_tgt
                if (strcmp(tide_info.name{tide_namei}, Model_harmonic.([tmp.testname, '_tname'])(coni,:))==1)
                    tide_info.index(tide_namei)=coni;
                end
            end
        end
        Model_amp.([tmp.testname, '_', tmp.yearstr]).con4= sum(Model_harmonic.([tmp.testname, '_tcon'])(tide_info.index(:),1));
        Model_amp.([tmp.testname, '_', tmp.yearstr]).M2=Model_harmonic.([tmp.testname, '_tcon'])(tide_info.index(1),1);
        Model_amp.([tmp.testname, '_', tmp.yearstr]).S2=Model_harmonic.([tmp.testname, '_tcon'])(tide_info.index(2),1);
        Model_amp.([tmp.testname, '_', tmp.yearstr]).K1=Model_harmonic.([tmp.testname, '_tcon'])(tide_info.index(3),1);
        Model_amp.([tmp.testname, '_', tmp.yearstr]).O1=Model_harmonic.([tmp.testname, '_tcon'])(tide_info.index(4),1);
%         Model_amp.([tmp.testname, '_', tmp.yearstr]).MU2=Model_harmonic.([tmp.testname, '_tcon'])(tide_info.index(5),1);
        a=1;
    end
end

save('D:\MEPL\project\SSH\7th_year(2022)\tide_byun\test_05\output.mat');

% tmp.ampm2(1)=Model_amp.test2127_2015.M2;
% tmp.ampm2(2)=Model_amp.test2127_2020.M2;
% tmp.ampm2(3)=Model_amp.test2127_2025.M2;
% tmp.ampm2(4)=Model_amp.test2127_2030.M2;
% 
% tmp.ampm22(1)=Model_amp.test2201_2015.M2;
% tmp.ampm22(2)=Model_amp.test2201_2020.M2;
% tmp.ampm22(3)=Model_amp.test2201_2025.M2;
% tmp.ampm22(4)=Model_amp.test2201_2030.M2;
% 
% tmp.ampm23(1)=Model_amp.test2201_2015.MU2;
% tmp.ampm23(2)=Model_amp.test2201_2020.MU2;
% tmp.ampm23(3)=Model_amp.test2201_2025.MU2;
% tmp.ampm23(4)=Model_amp.test2201_2030.MU2;