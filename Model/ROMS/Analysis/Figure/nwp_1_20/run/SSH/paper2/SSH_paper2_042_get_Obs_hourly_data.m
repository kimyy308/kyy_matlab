close all; clear all; clc;

dropboxpath='C:\Users\User\Dropbox';
% addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
% addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\t_tide_v1.3beta']));
% addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
% addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
% addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

% Define inference parameters.
%        infername=['P1';'K2'];
%        inferfrom=['K1';'S2'];
%        infamp=[.33093;.27215];
%        infphase=[-7.07;-22.40];

% 
% hold on
% for i=1:19
%     plot(zeta(i,:))
% end
% hold off

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

Obs.filename{1}='KHOA_OBS-L2_Anheung_SSH_Hourly_19890101_20181231.txt';
Obs.filename{2}='KHOA_OBS-L2_Gunsan_SSH_Hourly_19830101_20181231.txt';
Obs.filename{3}='KHOA_OBS-L2_Mokpo_SSH_Hourly_19600101_20181231.txt';
Obs.filename{4}='KHOA_OBS-L2_Heuksando_SSH_Hourly_19820101_20181231.txt';
Obs.filename{5}='KHOA_OBS-L2_Chujado_SSH_Hourly_19860101_20181231.txt';
Obs.filename{6}='KHOA_OBS-L2_Wando_SSH_Hourly_19850101_20181231.txt';
Obs.filename{7}='KHOA_OBS-L2_Yeosu_SSH_Hourly_19690101_20181231.txt';
Obs.filename{8}='KHOA_OBS-L2_Tongyeong_SSH_Hourly_19790101_20181231.txt';
Obs.filename{9}='KHOA_OBS-L2_Gadeokdo_SSH_Hourly_19830101_20181231.txt';
Obs.filename{10}='KHOA_OBS-L2_Busan_SSH_Hourly_19750101_20181231.txt';
Obs.filename{11}='Ieodo';
Obs.filename{12}='KHOA_OBS-L2_Busan_SSH_Hourly_19750101_20181231.txt';
Obs.filename{13}='KHOA_OBS-L2_Seogwipo_SSH_Hourly_19850101_20181231.txt';
Obs.filename{14}='KHOA_OBS-L2_Geomundo_SSH_Hourly_19850101_20181231.txt';
Obs.filename{15}='KHOA_OBS-L2_Ulsan_SSH_Hourly_19770101_20181231.txt';
Obs.filename{16}='KHOA_OBS-L2_Pohang_SSH_Hourly_19770101_20181231.txt';
Obs.filename{17}='KHOA_OBS-L2_Mukho_SSH_Hourly_19770101_20181231.txt';
Obs.filename{18}='KHOA_OBS-L2_Sokcho_SSH_Hourly_19770101_20181231.txt';
Obs.filename{19}='KHOA_OBS-L2_Ulleungdo_SSH_Hourly_19830101_20181231.txt';


for stai=1:num_sta
    tic
    if (stai~=11)
%         filename = ['Z:\내 드라이브\MEPL\project\SSH\4th_year\09_backup\UST\99_재생산자료\01_매시별\',Obs.filename{stai}];
        filename = ['F:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\09_backup\UST\99_재생산자료\01_매시별\',Obs.filename{stai}];

        startRow = 10;

        formatSpec = '%4s%3s%3s%3s%3s%9s%s%[^\n\r]';

        fileID = fopen(filename,'r');

        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

        fclose(fileID);

        raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
        for col=1:length(dataArray)-1
            raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
        end
        numericData = NaN(size(dataArray{1},1),size(dataArray,2));

        for col=[1,2,3,4,5,6,7]
            % 입력 셀형 배열의 텍스트를 숫자로 변환합니다. 숫자형이 아닌 텍스트를 NaN으로 바꿨습니다.
            rawData = dataArray{col};
            for row=1:size(rawData, 1)
                % 숫자형이 아닌 접두사 및 접미사를 검색하고 제거하는 정규 표현식을 만듭니다.
                regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                try
                    result = regexp(rawData(row), regexstr, 'names');
                    numbers = result.numbers;

                    % 천 단위가 아닌 위치에서 쉼표를 검색했습니다.
                    invalidThousandsSeparator = false;
                    if numbers.contains(',')
                        thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                        if isempty(regexp(numbers, thousandsRegExp, 'once'))
                            numbers = NaN;
                            invalidThousandsSeparator = true;
                        end
                    end
                    % 숫자형 텍스트를 숫자로 변환합니다.
                    if ~invalidThousandsSeparator
                        numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                        numericData(row, col) = numbers{1};
                        raw{row, col} = numbers{1};
                    end
                catch
                    raw{row, col} = rawData{row};
                end
            end
        end
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
        raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기
        Obs.data{stai} = cell2mat(raw);
        clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;
    end
    disp(num2str(stai)); toc
end

save('Z:\내 드라이브\MEPL\project\SSH\4th_year\09_backup\UST\99_재생산자료\01_매시별\tide_all.mat', 'station', 'Obs');
