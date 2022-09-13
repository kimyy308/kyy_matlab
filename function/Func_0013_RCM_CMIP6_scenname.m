function [scenname, error_status] = Func_0013_RCM_CMIP6_scenname(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=RCM_CMIP6_testname_his(testname);
%
% get the scenario name corresponding to SSP test name
%
%  input:
%  testname             ROMS SSP test name (string)
%
%  output:
%  scenname             CMIP6 scenario name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    6-Jul-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(testname)
    case {'test2102', 'test2103', 'test2104', 'test2105', 'test2106', 'v04', 'v05'}
        scenname='historical';
    case {'test2107', 'test2108', 'test2109', 'test2110', 'test2111'}
        scenname='ssp585';
    case {'test2117', 'test2118', 'test2119', 'test2120', 'test2121'}
        scenname='historical';
    case {'test2127', 'test2128', 'test2129', 'test2130', 'test2131'}
        scenname='ssp585';
end

error_status = 1;
end

