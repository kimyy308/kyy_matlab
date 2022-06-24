function [testname_his, error_status] = Func_0023_RCM_CMIP6_testname_his(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=RCM_CMIP6_testname_his(testname);
%
% get the historical testname corresponding to SSP test name
%
%  input:
%  testname             ROMS SSP test name (string)
%
%  output:
%  testname_his         ROMS historical test name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    18-Apr-2022 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch(testname)
        case {'test2117', 'test2127', 't17s_t27c'}
            testname_his='test2117';
        case {'test2118', 'test2128', 't18s_t28c'}
            testname_his='test2118';
        case {'test2119', 'test2129', 't19s_t29c'}
            testname_his='test2119';
        case {'test2120', 'test2130', 't20s_t30c'}
            testname_his='test2120';
        case {'test2121', 'test2131', 't21s_t31c'}
            testname_his='test2121';
        case {'t27s_t17c'}
            testname_his='test2127';
        case {'t28s_t18c'}
            testname_his='test2128';
        case {'t29s_t19c'}
            testname_his='test2129';
        case {'t30s_t20c'}
            testname_his='test2130';
        case {'t31s_t21c'}
            testname_his='test2131';
        case {'ens2201'}
            testname_his='ens2202';
        case {'ens2203'}
            testname_his='ens2204';
        case {'prob_ens2203'}
            testname_his='prob_ens2204';
           
    end
    error_status=1;
end

