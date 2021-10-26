function [scenname, error_status] = RCM_CMIP5_scenname(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=RCM_CMIP5_testname_his(testname);
%
% get the scenario name corresponding to RCP test name
%
%  input:
%  testname             ROMS RCP test name (string)
%
%  output:
%  scenname             CMIP5 scenario name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    3-May-2021 by Yong-Yub Kim
%  Updated    17-Oct-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(testname)
    case {'test53', 'test54', 'test55', 'test56', 'ens03'}
        scenname='historical';
    case {'test57', 'test58', 'test59', 'test60', 'ens08'}
        scenname='rcp45';
    case {'test61', 'test62', 'test63', 'test64', 'ens09'}
        scenname='rcp26';
    case {'test65', 'test66', 'test67', 'test68', 'ens10'}
        scenname='rcp85';
end

error_status=1;
end

