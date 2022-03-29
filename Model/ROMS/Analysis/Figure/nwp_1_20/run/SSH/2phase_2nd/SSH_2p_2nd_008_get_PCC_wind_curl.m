clc; close all; clear all;


RCM_info.years = 1985:2014;
% RCM_info.season='February';
RCM_info.season='all';
[RCM_info.months,error_status] = Func_0019_get_month_from_season(RCM_info.season);
tmp.regionname = 'ES_KHOA';
tmp.variable = 'wstrcurl';
RCM_info.name = {'test2117', 'test2118', 'test2119', 'test2120', 'test2121', 'v04'};

for testnameind2=1:length(RCM_info.name)
    tmp.testname=RCM_info.name{testnameind2};   % % need to change
    dirs.matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\mean\');

    tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
                '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
            '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
    load(tmp.matname)
    datas.(tmp.testname)=mean_data;
end

aaa= datas.v04;
% bbb= datas.test2117;
% bbb= datas.test2118;
% bbb= datas.test2119;
% bbb= datas.test2120;
bbb= datas.test2121;

aaa2=aaa(logical(~isnan(aaa).*~isnan(bbb)));
bbb2=bbb(logical(~isnan(aaa).*~isnan(bbb)));

corrcoef(aaa2, bbb2)


         