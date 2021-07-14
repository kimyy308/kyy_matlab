clear all;

%% �ؽ�Ʈ ���Ͽ��� �����͸� �����ɴϴ�.
% ���� �ؽ�Ʈ ���Ͽ��� �����͸� �������� ���� ��ũ��Ʈ:
%
%    D:\MEPL\project\SSH\�߰�����\smooth13_vtvs\color_test.txt
%
% ������ �ٸ� �����ͳ� �ؽ�Ʈ ���Ϸ� �ڵ带 Ȯ���Ϸ��� ��ũ��Ʈ ��� �Լ��� �����Ͻʽÿ�.

% MATLAB���� ���� ��¥�� �ڵ� ������: 2017/08/31 15:15:17

%% ������ �ʱ�ȭ�մϴ�.
filename = 'D:\MEPL\project\SSH\�߰�����\smooth13_vtvs\color_test.txt';
delimiter = ' ';
startRow = 3;

%% ������ ���� �ؽ�Ʈ�� ����:
% �ڼ��� ������ ���� �������� TEXTSCAN�� �����Ͻʽÿ�.
formatSpec = '%*q%*q%*q%*q%*q%*q%q%[^\n\r]';

%% �ؽ�Ʈ ������ ���ϴ�.
fileID = fopen(filename,'r');

%% ���Ŀ� ���� ������ ���� �н��ϴ�.
% �� ȣ���� �� �ڵ带 �����ϴ� �� ���Ǵ� ������ ����ü�� ������� �մϴ�. �ٸ� ���Ͽ� ���� ������ �߻��ϴ� ��� �������� ������
% �ڵ带 �ٽ� �����Ͻʽÿ�.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% �ؽ�Ʈ ������ �ݽ��ϴ�.
fclose(fileID);

%% ������ �ؽ�Ʈ�� �ִ� ���� ������ ���ڷ� ��ȯ�մϴ�.
% �������� �ƴ� �ؽ�Ʈ�� NaN���� �ٲߴϴ�.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));


%% �����͸� ���� �� �� ���� �����մϴ�.
rawNumericColumns = {};
rawCellColumns = raw(:, 1);


%% ������ �迭�� �� ���� �̸����� �Ҵ�
VarName1 = rawCellColumns(:, 1);


%% �ӽ� ���� �����
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawNumericColumns rawCellColumns;

fid=fopen('test.txt','w')
for i=1:length(VarName1)
fprintf(fid,'%s\n',char(VarName1(i)));
end
fclose(fid)

%% �ؽ�Ʈ ���Ͽ��� �����͸� �����ɴϴ�.
% ���� �ؽ�Ʈ ���Ͽ��� �����͸� �������� ���� ��ũ��Ʈ:
%
%    D:\MEPL\project\SSH\�߰�����\smooth13_vtvs\test.txt
%
% ������ �ٸ� �����ͳ� �ؽ�Ʈ ���Ϸ� �ڵ带 Ȯ���Ϸ��� ��ũ��Ʈ ��� �Լ��� �����Ͻʽÿ�.

% MATLAB���� ���� ��¥�� �ڵ� ������: 2017/08/31 15:26:37

%% ������ �ʱ�ȭ�մϴ�.
filename = 'D:\MEPL\project\SSH\�߰�����\smooth13_vtvs\test.txt';
delimiter = ' ';

%% ������ ���� �ؽ�Ʈ�� ����:
% �ڼ��� ������ ���� �������� TEXTSCAN�� �����Ͻʽÿ�.
formatSpec = '%s%s%s%s%[^\n\r]';

%% �ؽ�Ʈ ������ ���ϴ�.
fileID = fopen(filename,'r');

%% ���Ŀ� ���� ������ ���� �н��ϴ�.
% �� ȣ���� �� �ڵ带 �����ϴ� �� ���Ǵ� ������ ����ü�� ������� �մϴ�. �ٸ� ���Ͽ� ���� ������ �߻��ϴ� ��� �������� ������
% �ڵ带 �ٽ� �����Ͻʽÿ�.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

%% �ؽ�Ʈ ������ �ݽ��ϴ�.
fclose(fileID);

%% ������ �ؽ�Ʈ�� �ִ� ���� ������ ���ڷ� ��ȯ�մϴ�.
% �������� �ƴ� �ؽ�Ʈ�� NaN���� �ٲߴϴ�.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4]
    % �Է� ���� �迭�� �ؽ�Ʈ�� ���ڷ� ��ȯ�մϴ�. �������� �ƴ� �ؽ�Ʈ�� NaN���� �ٲ���ϴ�.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % �������� �ƴ� ���λ� �� ���̻縦 �˻��ϰ� �����ϴ� ���� ǥ������ ����ϴ�.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % õ ������ �ƴ� ��ġ���� ��ǥ�� �˻��߽��ϴ�.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % ������ �ؽ�Ʈ�� ���ڷ� ��ȯ�մϴ�.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% ������ �迭�� �� ���� �̸����� �Ҵ�
R = cell2mat(raw(:, 1));
G = cell2mat(raw(:, 2));
B = cell2mat(raw(:, 3));
A = cell2mat(raw(:, 4));


%% �ӽ� ���� �����
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me;

jet_mod(:,1)=R(31:-1:1)./255;
jet_mod(:,2)=G(31:-1:1)./255;
jet_mod(:,3)=B(31:-1:1)./255;
figure
colorbar
colormap(jet_mod)