function [abbname, error_status] = Func_0018_abbname(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [abbname, error_status] = Func_0018_abbname(testname)
%
% get the abbreviation name corresponding to the testname
%
%  input:
%  testname             RCM or GCM testname (string)
%
%  output:
%  abbname             abbreviation name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    17-Jan-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(testname)
    case{'test06'}
        abbname = 'Reana';
    case {'test2117', 'test2127'}
        abbname = 'RCM-CNE';
    case {'test2118', 'test2128'}
        abbname = 'RCM-ECV';
    case {'test2119', 'test2129'}
        abbname = 'RCM-ACC';
    case {'test2120', 'test2130'}
        abbname = 'RCM-CNH';
    case {'test2121', 'test2131'}
        abbname = 'RCM-CMC';
    case {'ens2201', 'ens2202', 'ENS4_hist', 'ENS4_fut'}
        abbname = 'RCM-ENS';
    case {'CNRM-ESM2-1'}
        abbname = 'GCM-CNE';
    case {'EC-Earth3-Veg'}
        abbname = 'GCM-ECV';
    case {'ACCESS-CM2'}
        abbname = 'GCM-ACC';
    case {'CNRM-CM6-1-HR'}
        abbname = 'GCM-CNH';
    case {'CMCC-ESM2'}
        abbname = 'GCM-CMC';
    case {'v04'}
        abbname = 'RCM-SOD';
    case {'v05'}
        abbname = 'RCM-GRS';
    case {'SODA'}
        abbname = 'GCM-SOD';
    case {'GLORYS'}
        abbname = 'GCM-GRS';
end

error_status=1;
end

