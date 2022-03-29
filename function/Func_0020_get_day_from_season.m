function [days,error_status] = Func_0020_get_day_from_season(months, flag_leap_year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [abbname, error_status] = Func_0019_get_month_from_season(testname)
%
% get the abbreviation name corresponding to the testname
%
%  input:
%  season             season or month (string)
%
%  output:
%  months             month (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    04-Feb-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%flag_leap_year from years
% flag_leap_year=year2day(year)-365;

%% calculation
for mind=1:length(months)
    tempmonth=months(mind);
    if flag_leap_year==1
        mdays_s=datenum(0,tempmonth,1);
        mdays_e=datenum(0,tempmonth+1,1)-1;
    else
        mdays_s=datenum(1,tempmonth,1)-366;
        mdays_e=datenum(1,tempmonth+1,1)-1-366;
    end
    if mind==1
        days = mdays_s : mdays_e;
    else
        days = [days, mdays_s:mdays_e];
    end
end

error_status=1;
end

