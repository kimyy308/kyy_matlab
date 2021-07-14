clear all;

%% 텍스트 파일에서 데이터를 가져옵니다.
% 다음 텍스트 파일에서 데이터를 가져오기 위한 스크립트:
%
%    D:\MEPL\project\SSH\중간보고\smooth13_vtvs\color_test.txt
%
% 선택한 다른 데이터나 텍스트 파일로 코드를 확장하려면 스크립트 대신 함수를 생성하십시오.

% MATLAB에서 다음 날짜에 자동 생성됨: 2017/08/31 15:15:17

%% 변수를 초기화합니다.
filename = 'D:\MEPL\project\SSH\중간보고\smooth13_vtvs\color_test.txt';
delimiter = ' ';
startRow = 3;

%% 데이터 열을 텍스트로 읽음:
% 자세한 내용은 도움말 문서에서 TEXTSCAN을 참조하십시오.
formatSpec = '%*q%*q%*q%*q%*q%*q%q%[^\n\r]';

%% 텍스트 파일을 엽니다.
fileID = fopen(filename,'r');

%% 형식에 따라 데이터 열을 읽습니다.
% 이 호출은 이 코드를 생성하는 데 사용되는 파일의 구조체를 기반으로 합니다. 다른 파일에 대한 오류가 발생하는 경우 가져오기 툴에서
% 코드를 다시 생성하십시오.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% 텍스트 파일을 닫습니다.
fclose(fileID);

%% 숫자형 텍스트가 있는 열의 내용을 숫자로 변환합니다.
% 숫자형이 아닌 텍스트를 NaN으로 바꿉니다.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));


%% 데이터를 숫자 및 셀 열로 분할합니다.
rawNumericColumns = {};
rawCellColumns = raw(:, 1);


%% 가져온 배열을 열 변수 이름으로 할당
VarName1 = rawCellColumns(:, 1);


%% 임시 변수 지우기
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawNumericColumns rawCellColumns;

fid=fopen('test.txt','w')
for i=1:length(VarName1)
fprintf(fid,'%s\n',char(VarName1(i)));
end
fclose(fid)

%% 텍스트 파일에서 데이터를 가져옵니다.
% 다음 텍스트 파일에서 데이터를 가져오기 위한 스크립트:
%
%    D:\MEPL\project\SSH\중간보고\smooth13_vtvs\test.txt
%
% 선택한 다른 데이터나 텍스트 파일로 코드를 확장하려면 스크립트 대신 함수를 생성하십시오.

% MATLAB에서 다음 날짜에 자동 생성됨: 2017/08/31 15:26:37

%% 변수를 초기화합니다.
filename = 'D:\MEPL\project\SSH\중간보고\smooth13_vtvs\test.txt';
delimiter = ' ';

%% 데이터 열을 텍스트로 읽음:
% 자세한 내용은 도움말 문서에서 TEXTSCAN을 참조하십시오.
formatSpec = '%s%s%s%s%[^\n\r]';

%% 텍스트 파일을 엽니다.
fileID = fopen(filename,'r');

%% 형식에 따라 데이터 열을 읽습니다.
% 이 호출은 이 코드를 생성하는 데 사용되는 파일의 구조체를 기반으로 합니다. 다른 파일에 대한 오류가 발생하는 경우 가져오기 툴에서
% 코드를 다시 생성하십시오.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

%% 텍스트 파일을 닫습니다.
fclose(fileID);

%% 숫자형 텍스트가 있는 열의 내용을 숫자로 변환합니다.
% 숫자형이 아닌 텍스트를 NaN으로 바꿉니다.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4]
    % 입력 셀형 배열의 텍스트를 숫자로 변환합니다. 숫자형이 아닌 텍스트를 NaN으로 바꿨습니다.
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
            % 숫자형 텍스트를 숫자로 변환합니다.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% 가져온 배열을 열 변수 이름으로 할당
R = cell2mat(raw(:, 1));
G = cell2mat(raw(:, 2));
B = cell2mat(raw(:, 3));
A = cell2mat(raw(:, 4));


%% 임시 변수 지우기
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me;

jet_mod(:,1)=R(31:-1:1)./255;
jet_mod(:,2)=G(31:-1:1)./255;
jet_mod(:,3)=B(31:-1:1)./255;
figure
colorbar
colormap(jet_mod)