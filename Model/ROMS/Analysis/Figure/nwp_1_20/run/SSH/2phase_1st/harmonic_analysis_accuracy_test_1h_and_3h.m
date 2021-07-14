close all; clc; clear all;

dropboxpath='C:\Users\User\Dropbox';
% addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
% addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\t_tide_v1.3beta']));

%% 텍스트 파일에서 데이터를 가져옴
% 다음 텍스트 파일에서 데이터를 가져오기 위한 스크립트:
%
%    파일 이름: Z:\내 드라이브\Data\Observation\Tidal gauge\인천_DT_78_2020_KR.txt
%
% MATLAB에서 2021-03-23 10:09:33에 자동 생성됨

%% 가져오기 옵션을 설정하고 데이터 가져오기
opts = delimitedTextImportOptions("NumVariables", 14, "Encoding", "UTF-8");

% 범위 및 구분 기호 지정
opts.DataLines = [5, Inf];
opts.Delimiter = "\t";

% 열 이름과 유형 지정
opts.VariableNames = ["VarName1", "cm", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14"];
opts.SelectedVariableNames = ["VarName1", "cm"];
opts.VariableTypes = ["char", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% 파일 수준 속성 지정
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 변수 속성 지정
opts = setvaropts(opts, ["VarName1", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14"], "EmptyFieldRule", "auto");

% 데이터 가져오기
DT782020KR = readtable("Z:\내 드라이브\Data\Observation\Tidal gauge\인천_DT_78_2020_KR.txt", opts);

%% 출력 유형으로 변환
DT782020KR = table2cell(DT782020KR);
numIdx = cellfun(@(x) ~isnan(str2double(x)), DT782020KR);
DT782020KR(numIdx) = cellfun(@(x) {str2double(x)}, DT782020KR(numIdx));

%% 임시 변수 지우기
clear opts

tide_rawdata=cell2mat(DT782020KR(:,2));
tide_1h_data_temp=reshape(tide_rawdata, [60, length(tide_rawdata)/60]);
tide_1h_data=tide_1h_data_temp(1,:);
tide_3h_data_temp=reshape(tide_rawdata, [480, length(tide_rawdata)/480]);
tide_3h_data=tide_3h_data_temp(3,:);

[tname_1h,tfreq_1h,tcon_1h,tout_1h]=t_tide(tide_1h_data,...
       'interval',1, ...                     % hourly data
       'start',datenum(2020,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
       'latitude',37.451944,...               % Latitude of Model_fy
       'rayleigh',1, 'error','wboot');                       % Use SNR=1 for synthesis. 

[tname_3h,tfreq_3h,tcon_3h,tout_3h]=t_tide(tide_3h_data,...
       'interval',3, ...                     % hourly data
       'start',datenum(2020,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
       'latitude',37.451944,...               % Latitude of Model_fy
       'rayleigh',1, 'error','wboot');                       % Use SNR=1 for synthesis. 


tide_info.name{1}='M2  ';
tide_info.name{2}='S2  ';
tide_info.name{3}='K1  ';
tide_info.name{4}='O1  ';

num_tide_all=size(tname_1h,1);
num_tide_tgt=length(tide_info.name);
for coni=1:num_tide_all
    for tide_namei=1:num_tide_tgt
        if (strcmp(tide_info.name{tide_namei}, tname_1h(coni,:))==1)
            tide_info.index(tide_namei)=coni;
        end
    end
end

tide_1h_amp_M2=tcon_1h(tide_info.index(1),1);
tide_1h_amp_S2=tcon_1h(tide_info.index(2),1);
tide_1h_amp_K1=tcon_1h(tide_info.index(3),1);
tide_1h_amp_O1=tcon_1h(tide_info.index(4),1);


num_tide_all=size(tname_3h,1);
num_tide_tgt=length(tide_info.name);
for coni=1:num_tide_all
    for tide_namei=1:num_tide_tgt
        if (strcmp(tide_info.name{tide_namei}, tname_3h(coni,:))==1)
            tide_info.index(tide_namei)=coni;
        end
    end
end

tide_3h_amp_M2=tcon_3h(tide_info.index(1),1);
tide_3h_amp_S2=tcon_3h(tide_info.index(2),1);
tide_3h_amp_K1=tcon_3h(tide_info.index(3),1);
tide_3h_amp_O1=tcon_3h(tide_info.index(4),1);