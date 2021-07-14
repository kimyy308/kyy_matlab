clear all; clc; close all;
% etopo_lon=ncread('C:\Users\kyy\Desktop\nwp\Roms_tools\Topo\etopo1.nc','lon');
% etopo_lat=double(ncread('C:\Users\kyy\Desktop\nwp\Roms_tools\Topo\etopo1.nc','lat'));
% etopo_depth=ncread('C:\Users\kyy\Desktop\nwp\Roms_tools\Topo\etopo1.nc','z');
% [Lat Lon] =meshgrid(etopo_lat(5950:9094)',etopo_lon(17317:21057)');
% xx=(110+1/240 : 1/120 : 170-1/240);
% yy=(10+1/240 : 1/120 : 60-1/240);
% depth(:,:)=griddata(Lon,Lat,etopo_depth(17317:21057,5950:9094),xx,yy');
load('etopo_interp.mat')

cd C:\Users\kyy\Desktop\KorBathy30s\KorBathy30s
list=ls;
cd C:\Users\kyy\Desktop\KorBathy30s
i=3;
for i=3:288
    templat=list(i,2:3)
    templon=list(i,5:7)
    filename=strcat('C:\Users\kyy\Desktop\KorBathy30s\KorBathy30s\',list(i,:))
%     filename = 'C:\Users\kyy\Desktop\KorBathy30s\KorBathy30s\N30E120.dpt';
    startRow = 2;
    endRow = 2;
    formatSpec = '%7s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
    fclose(fileID);

    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121]
        % 입력 셀형 배열의 문자열을 숫자로 변환합니다. 숫자형이 아닌 문자열을 NaN으로 바꿨습니다.
        rawData = dataArray{col};
        for row=1:size(rawData, 1);
            % 숫자형이 아닌 접두사 및 접미사를 검색하고 제거하는 정규 표현식을 만듭니다.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;

                % 천 단위가 아닌 위치에서 쉼표를 검색했습니다.
                invalidThousandsSeparator = false;
                if any(numbers==',');
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'));
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % 숫자 문자열을 숫자로 변환합니다.
                if ~invalidThousandsSeparator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch me
            end
        end
    end
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
    raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기
    var=cell2mat(raw(:, 1:120));
    clearvars filename formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
%     for iii=1:120
%         for jjj=1:120
% % % %        i -> lat , j -> lon
% templat='30'; templon='130';
            istart=(str2num(templat)-10)*120+1;
            iend=(str2num(templat)-9)*120;
            jstart=(str2num(templon)-110)*120+1;
            jend=(str2num(templon)-109)*120;
            templat
            templon
            xx(jstart)
            xx(jend)
            yy(istart)
            yy(iend)
            depth(istart:iend,jstart:jend) = 0-var(120:-1:1,1:120);
%         end
%     end
end

depth(find(depth(:,:)>10000))=5000;
depth(find(depth(:,:)<-10000))=depth(find(depth(:,:)<-10000))/1000.;