function [months,error_status] = Func_0019_get_month_from_season(season)
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
%  Updated    20-Jan-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(season)
    case 'all'
        months =1:12;  
    case 'spring'
        months =[3,4,5];  
    case 'summer'
        months =[6,7,8];  
    case 'fall'
        months =[9,10,11];  
    case 'winter'
        months =[1,2,12];  
    case 'January'
        months =1;
    case 'February'
        months =2;
    case 'March'
        months =3;
    case 'April'
        months =4;
    case 'May'
        months =5;
    case 'June'
        months =6;
    case 'July'
        months =7;
    case 'August'
        months =8;
    case 'September'
        months =9;
    case 'October'
        months =10;
    case 'November'
        months =11;
    case 'December'
        months =12;
end
error_status=1;
end

