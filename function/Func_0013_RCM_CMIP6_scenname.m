function [scenname, error_status] = RCM_CMIP6_scenname(testname)
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
    case {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
        scenname='historical';
end

error_status = 1;
end

