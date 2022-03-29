disp('subroutine SSH_2p_2nd_000_sub_000_time_set') 

tmp.tind=1;
for yearij = 1:length(RCM_info.years)
    for month=1:length(RCM_info.months) 
        tmp.year = RCM_info.years(yearij);
%                     tmp.month = RCM_info.months(monthij);
        RCM_time.ftime(tmp.tind) = datenum(tmp.year,month,15) - datenum(1900,12,31);
        tmp.tind=tmp.tind+1;
    end
end
for month=1:length(RCM_info.months)  
    tmp.year = RCM_info.years(yearij);
    RCM_time.climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
end

for i =1:length(RCM_info.years) 
    tmp.year = RCM_info.years(yearij);
    for month=1:length(RCM_info.months)  
        RCM_time.xData((12*(i-1))+month) = datenum([num2str(tmp.year),'-',num2str(month,'%02i'),'-01',]); 
    end
end

RCM_time.trendtime=RCM_info.years(1):1/length(RCM_info.months) : RCM_info.years(end)+1-1/length(RCM_info.months) ;
RCM_time.trendtime_yearly=RCM_info.years(1) : RCM_info.years(end);