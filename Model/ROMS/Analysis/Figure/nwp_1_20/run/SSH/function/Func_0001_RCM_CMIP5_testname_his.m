function [testname_his] = RCM_CMIP5_testname_his(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=RCM_CMIP5_testname_his(testname);
%
% get the historical testname corresponding to RCP test name
%
%  input:
%  testname             ROMS RCP test name (string)
%
%  output:
%  testname_his         ROMS historical test name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    1-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(testname)
    case {'test57', 'test61', 'test65'}
        testname_his='test53';
    case {'test58', 'test62', 'test66'}
        testname_his='test54';
    case {'test59', 'test63', 'test67'}
        testname_his='test55';
    case {'test60', 'test64', 'test68'}
        testname_his='test56';
end

end

