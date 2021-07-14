function [rcmtestname, error_status] = Func_0006_get_RCMname_from_GCM(testname, scenname)
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
    case {'IPSL-CM5A-LR'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test53';
            case {'rcp26'}
                rcmtestname='test61';
            case {'rcp45'}
                rcmtestname='test57';
            case {'rcp85'}
                rcmtestname='test65';
        end
    case {'IPSL-CM5A-MR'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test54';
            case {'rcp26'}
                rcmtestname='test62';
            case {'rcp45'}
                rcmtestname='test58';
            case {'rcp85'}
                rcmtestname='test66';
        end
    case {'NorESM1-M'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test55';
            case {'rcp26'}
                rcmtestname='test63';
            case {'rcp45'}
                rcmtestname='test59';
            case {'rcp85'}
                rcmtestname='test67';
        end
    case {'MPI-ESM-LR'}
        switch(scenname)
            case {'historical'}
                rcmtestname='test56';
            case {'rcp26'}
                rcmtestname='test64';
            case {'rcp45'}
                rcmtestname='test60';
            case {'rcp85'}
                rcmtestname='test68';
        end
end

error_status = 1;
end

