function [gcmtestname] = get_GCMname_from_RCM(testname)
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
end

end

