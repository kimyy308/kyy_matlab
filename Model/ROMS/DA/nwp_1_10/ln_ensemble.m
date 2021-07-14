close all;clear all;clc
% system('ln -s /dev/ttyACM0 /dev/ttyS101')
filepath = '/HDD1/kimyy/reanalysis_data/ensemble/ens_v01/';
ens = 1;
j = 1;
% len_day = 12412;
len_day = 13448;
% len_day = 365;
for i = 12412:7:len_day
    initdaynum=datenum(1981,12,31);
    nowdate=datestr(initdaynum+i,'yyyymmdd');
    year=nowdate(1:4);
    numyear=str2num(year);
    month=nowdate(5:6);
    day2=nowdate(7:8);
    d=datetime([numyear str2num(month) str2num(day2)]);
    yearday=day(d,'dayofyear');
    if (yearday==366)
        yearday=365;
    end
    %     if (mod(i-365*(ens-1),3)==1)
    ens =1; 
    for j=1982:1982+33
        file1 = ['/HDD1/kimyy/Model/ROMS/roms_nwp/nwp_1_10/output/test06/run/', ...
            num2str(j),'/ocean_avg_',num2str(yearday,'%4.4d'),'.nc'];
        file2 = [filepath,'ocean_ens',num2str(ens,'%2.2d'),'_',num2str(i,'%5.5d'),'.nc'];
        eval(['system(''ln -sf ',file1,' ',file2,''')'])
%         j = j+1;
        ens=ens+1;
    end
    
%     end
%     if (mod(i,365) ==0)
%         ens = ens+1;
%         j=1;
%     end
end

% % j=1;
% % ens = 11;
% for i = 1:3650
%     if (i<21)
%         ii=i+3650;
%     else
%         ii=i;
%     end
% %     if (mod(i-365*(ens-11),3)==1)
%         file1 = ['ocean_avg_',num2str(ii-20,'%4.4d'),'.nc'];
%         file2 = [filepath,'ocean_ens',num2str(ens,'%2.2d'),'_',num2str(j,'%5.5d'),'.nc'];
%         eval(['system(''ln ',file1,' ',file2,''')']);
%         j = j+1;
% %     end
%     if (mod(i,365) ==0)
%         ens = ens+1;
%         j=1;
%     end
% end
% 
% for i = 1:3650
%     if (i>3630)
%         ii=i-3650;
%     else
%         ii=i;
%     end
% %     if (mod(i-365*(ens-21),3)==1)
%         file1 = ['ocean_avg_',num2str(ii+20,'%4.4d'),'.nc'];
%         file2 = [filepath,'ocean_ens',num2str(ens,'%2.2d'),'_',num2str(j,'%3.3d'),'.nc'];
%         eval(['system(''ln ',file1,' ',file2,''')']);
%         j = j+1;
% %     end
%     if (mod(i,365) ==0)
%         ens = ens+1;
%         j=1;
%     end
% end
% 
% %% static ensemble
% % close all;clear all;clc
% % % system('ln -s /dev/ttyACM0 /dev/ttyS101')
% % filepath = '/data1/kimyy/Data_Assimilation/EnOI/ensemble/ens05/';
% % ens = 1;
% % for i = 1825:60:3650
% %     file1 = ['ocean_avg_',num2str(i,'%4.4d'),'.nc'];
% %     file2 = [filepath,'ocean_ens',num2str(ens,'%2.2d'),'.nc'];
% %     eval(['system(''ln ',file1,' ',file2,''')']);
% %     ens = ens+1;
% % end
