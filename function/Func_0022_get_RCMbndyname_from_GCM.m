function [rcmtestname, error_status] = Func_0022_get_RCMbndyname_from_GCM(testname, scenname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=get_RCMname_from_GCM(testname, scenname);
%
% get the GCM testname corresponding to RCM test name
%
%  input:
%  testname             ROMS RCM test name (string)
%  scenname             scenario name (string)
%
%  output:
%  rcmtestname          ROMS RCM test name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    15-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(testname)
    case {'CNRM-ESM2-1'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test2102';
            case {'ssp585'}
                rcmtestname='test2107';
        end
    case {'EC-Earth3-Veg'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test2103';
            case {'ssp585'}
                rcmtestname='test2108';
        end
    case {'ACCESS-CM2'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test2104';
            case {'ssp585'}
                rcmtestname='test2109';
        end
    case {'CNRM-CM6-1-HR'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test2105';
            case {'ssp585'}
                rcmtestname='test2110';
        end
    case {'CMCC-ESM2'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test2106';
            case {'ssp585'}
                rcmtestname='test2111';
        end
end

error_status = 1;
end

