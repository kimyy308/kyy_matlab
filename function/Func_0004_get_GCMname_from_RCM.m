function [gcmtestname] = Func_0004_get_GCMname_from_RCM(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=get_GCMname_from_RCM(testname);
%
% get the GCM testname corresponding to RCM test name
%
%  input:
%  testname             ROMS RCM test name (string)
%
%  output:
%  gcmtestname          CMIP5 GCM test name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    14-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(testname)
    case {'test53', 'test57', 'test61', 'test65'}
        gcmtestname='IPSL-CM5A-LR';
    case {'test54', 'test58', 'test62', 'test66'}
        gcmtestname='IPSL-CM5A-MR';
    case {'test55', 'test59', 'test63', 'test67'}
        gcmtestname='NorESM1-M';
    case {'test56', 'test60', 'test64', 'test68'}
        gcmtestname='MPI-ESM-LR';
    case {'test2102', 'test2107', 'test2117', 'test2127', 'test2210'}
        gcmtestname='CNRM-ESM2-1';
    case {'test2103', 'test2108', 'test2118', 'test2128', 'test2211', 'test2213'}
        gcmtestname='EC-Earth3-Veg';
    case {'test2104', 'test2109', 'test2119', 'test2129', 'test2208'}
        gcmtestname='ACCESS-CM2';
    case {'test2105', 'test2110', 'test2120', 'test2130', 'test2212'}
        gcmtestname='CNRM-CM6-1-HR';
    case {'test2106', 'test2111', 'test2121', 'test2131', 'test2209'}
        gcmtestname='CMCC-ESM2';
end

end

