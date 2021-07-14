close all; clear all; clc;

% 
% uniqline = unique(line_st,'stable');
% for i=1:length(uniqline)
%     uniqind=find(line_st==uniqline(i),1);
%     uniqlon(i)=lon(uniqind);
%     uniqlat(i)=lat(uniqind);
% end
% 
% 
% YS -> Line 311  --> find(line_st>=31100 & line_st<31200)
% SS -> Line 314
% KS -> Line 206
% ES -> Line 105
% 
% (line,station,depth,year,month)

year=1980:1985;
month=2:2:12;
line=[107 106 105 104 103 102 209 208 ...
    207 206 400 205 204 203 317 316 315 314 313 ...
    312 311 310 309 308 307]; 
% line=[ ...
%     207 206 400 205 204 203 314 313]; 
% maxdepth=75;
stddepth=[0 10 20 30 50 75 100 125 150 200 250 300 400 500];

for numyear=1:length(year)
    tempyear=year(numyear)
    [~, ~, NIFS_raw] = xlsread(['E:\Data\Observation\NIFS\xls\NIFS_',num2str(tempyear,'%04i'),'.xls']);

    NIFS_raw = string(NIFS_raw);
    NIFS_raw(ismissing(NIFS_raw)) = '';

    for i=3:length(NIFS_raw)
        NIFS.line(i-2)=str2double(NIFS_raw{i,2});
        NIFS.station(i-2)=str2double(NIFS_raw{i,3});
        NIFS.lat(i-2)=str2double(NIFS_raw{i,5});
        NIFS.lon(i-2)=str2double(NIFS_raw{i,6});
        NIFS.year(i-2)=str2double(NIFS_raw{i,7}(1:4));
        NIFS.month(i-2)=str2double(NIFS_raw{i,7}(6:7));
        NIFS.day(i-2)=str2double(NIFS_raw{i,7}(9:10));
        if (mod(NIFS.month(i-2),2) == 1 && NIFS.day(i-2) <= 15)
            NIFS.month2(i-2) = NIFS.month(i-2) - 1;
            if (NIFS.month2(i-2) == 0)
                NIFS.month2(i-2) = 2;
            end
        else
            NIFS.month2(i-2) = NIFS.month(i-2);
        end
        NIFS.depth(i-2)=str2double(NIFS_raw{i,8});

        NIFS.temp_flag(i-2)=str2double(NIFS_raw{i,10});
        if (NIFS.temp_flag(i-2)==1)
            NIFS.temp(i-2)=str2double(NIFS_raw{i,9});
        else
            NIFS.temp(i-2)=NaN;
        end
        NIFS.salt_flag(i-2)=str2double(NIFS_raw{i,12});
        if (NIFS.salt_flag(i-2)==1)
            NIFS.salt(i-2)=str2double(NIFS_raw{i,11});
        else
            NIFS.salt(i-2)=NaN;
        end
    end
    for numline=1:length(line)
        templine=line(numline);
        
        line_ind=find(NIFS.line==templine);
        NIFS2.line=NIFS.line(line_ind);
        NIFS2.station=NIFS.station(line_ind);
        NIFS2.lat=NIFS.lat(line_ind);
        NIFS2.lon=NIFS.lon(line_ind);
        NIFS2.year=NIFS.year(line_ind);
        NIFS2.month2=NIFS.month2(line_ind);
        NIFS2.day=NIFS.day(line_ind);
        NIFS2.depth=NIFS.depth(line_ind);
        NIFS2.temp=NIFS.temp(line_ind);
        NIFS2.salt=NIFS.salt(line_ind);
        
        
%         line=[107 106 105 104 103 102 209 208 ...
%     207 206 400 205 204 203 314 313 ...
%     312 311 310 309 308 307]; 
%         switch(templine)
%             case(107)
%                 station=1:7;
%             case(106)
%                 station=2:11;
%             case(105)
%                 station=3:15;
%             case(311)
%                 station=5:10;
%             case(314)
%                 station=1:10;
%             case(206)
%                 station=3:5;
%             case(105)
%                 station=4:15;
%         end
        station=1:30;
        
        for numsta=1:length(station)
            tempstation=station(numsta);
            sta_ind=find(NIFS2.station==tempstation);
            NIFS3.line=NIFS2.line(sta_ind);
            NIFS3.station=NIFS2.station(sta_ind);
            NIFS3.lat=NIFS2.lat(sta_ind);
            NIFS3.lon=NIFS2.lon(sta_ind);
            NIFS3.year=NIFS2.year(sta_ind);
            NIFS3.month2=NIFS2.month2(sta_ind);
            NIFS3.day=NIFS2.day(sta_ind);
            NIFS3.depth=NIFS2.depth(sta_ind);
            NIFS3.temp=NIFS2.temp(sta_ind);
            NIFS3.salt=NIFS2.salt(sta_ind);

%             if (max(NIFS3.depth)<75)
        % %         skip station
%             else
                for nummon=1:length(month)
                    tempmonth=month(nummon);
                    mon_ind=find(NIFS3.month2==tempmonth);
                    NIFS4.line=NIFS3.line(mon_ind);
                    NIFS4.station=NIFS3.station(mon_ind);
                    NIFS4.lat=NIFS3.lat(mon_ind);
                    NIFS4.lon=NIFS3.lon(mon_ind);
                    NIFS4.year=NIFS3.year(mon_ind);
                    NIFS4.month2=NIFS3.month2(mon_ind);
                    NIFS4.day=NIFS3.day(mon_ind);
                    NIFS4.depth=NIFS3.depth(mon_ind);
                    NIFS4.temp=NIFS3.temp(mon_ind);
                    NIFS4.salt=NIFS3.salt(mon_ind);
                    
                    uniqdepth = unique(NIFS4.depth,'stable');
                    for i=1:length(uniqdepth)
                        depth_ind=find(NIFS4.depth==uniqdepth(i),1);
                        NIFS5.line(i)=NIFS4.line(depth_ind);
                        NIFS5.station(i)=NIFS4.station(depth_ind);
                        NIFS5.lat(i)=NIFS4.lat(depth_ind);
                        NIFS5.lon(i)=NIFS4.lon(depth_ind);
                        NIFS5.year(i)=NIFS4.year(depth_ind);
                        NIFS5.month2(i)=NIFS4.month2(depth_ind);
                        NIFS5.day(i)=NIFS4.day(depth_ind);
                        NIFS5.depth(i)=NIFS4.depth(depth_ind);
                        NIFS5.temp(i)=NIFS4.temp(depth_ind);
                        NIFS5.salt(i)=NIFS4.salt(depth_ind);
                    end
                    
                    if (exist('NIFS5')==1)
                        depthlen=length(stddepth(stddepth<=max(NIFS5.depth)));
                        NIFS6.depth(numyear,numline,numsta,nummon,1:depthlen)=stddepth(stddepth<=max(NIFS5.depth));
                        NIFS6.temp(numyear,numline,numsta,nummon,1:depthlen)=interp1(NIFS5.depth,NIFS5.temp,NIFS6.depth(numyear,numline,numsta,nummon,1:depthlen));
                        NIFS6.temp(NIFS6.temp==0)=NaN;
                        NIFS6.salt(numyear,numline,numsta,nummon,1:depthlen)=interp1(NIFS5.depth,NIFS5.salt,NIFS6.depth(numyear,numline,numsta,nummon,1:depthlen));
                        NIFS6.salt(NIFS6.salt==0)=NaN;
                        NIFS6.salt75_sta(numyear,numline,numsta,nummon)=mean(NIFS6.salt(numyear,numline,numsta,nummon,1:depthlen),5,'omitnan');
                        NIFS6.salt_surf_sta(numyear,numline,numsta,nummon)=NIFS6.salt(numyear,numline,numsta,nummon,1);
                        NIFS6.year(numyear,numline,numsta,nummon,1:depthlen)=NIFS5.year(1);
                        NIFS6.month2(numyear,numline,numsta,nummon,1:depthlen)=NIFS5.month2(1);
                        NIFS6.line(numyear,numline,numsta,nummon,1:depthlen)=NIFS5.line(1);
                        NIFS6.station(numyear,numline,numsta,nummon,1:depthlen)=NIFS5.station(1);
                        NIFS6.lon(numyear,numline,numsta,nummon,1:depthlen)=NIFS5.lon(1);
                        NIFS6.lat(numyear,numline,numsta,nummon,1:depthlen)=NIFS5.lat(1);

                        clear NIFS5;
                    else
                        NIFS6.depth(numyear,numline,numsta,nummon,:)=NaN;
                        NIFS6.temp(numyear,numline,numsta,nummon,:)=NaN;
                        NIFS6.salt(numyear,numline,numsta,nummon,:)=NaN;
                        NIFS6.salt75_sta(numyear,numline,numsta,nummon)=NaN;
                        NIFS6.salt_surf_sta(numyear,numline,numsta,nummon)=NaN;
                        NIFS6.year(numyear,numline,numsta,nummon,:)=NaN;
                        NIFS6.month2(numyear,numline,numsta,nummon,:)=NaN;
                        NIFS6.line(numyear,numline,numsta,nummon,:)=NaN;
                        NIFS6.station(numyear,numline,numsta,nummon,:)=NaN;
                        if (isnan(sum(NIFS3.lon,'omitnan'))==1 || sum(NIFS3.lon,'omitnan')==0)
                            NIFS6.lon(numyear,numline,numsta,nummon,1:length(stddepth))=NaN;
                            NIFS6.lat(numyear,numline,numsta,nummon,1:length(stddepth))=NaN;
                        else
                            NIFS6.lon(numyear,numline,numsta,nummon,1:length(stddepth))=NIFS3.lon(1);
                            NIFS6.lat(numyear,numline,numsta,nummon,1:length(stddepth))=NIFS3.lat(1);
                        end
                    end
                end
%             end
        end
    end
    % % year, line, station, month, depth
end

NIFS6.salt75_sta(NIFS6.salt75_sta==0)=NaN;
NIFS6.salt_surf_sta(NIFS6.salt_surf_sta==0)=NaN;
NIFS6.salt75_line=squeeze(mean(NIFS6.salt75_sta,3,'omitnan'));
NIFS6.salt_surf_line=squeeze(mean(NIFS6.salt_surf_sta,3,'omitnan'));

numi=1;
for i=1:size(NIFS6.salt75_line,1)  %% year
    for j=1:size(NIFS6.salt75_line,3)  %% month
%         NIFS6.salt75_line_comb=reshape(NIFS6.salt75_line,size(NIFS6.salt75_line,1)*size(NIFS6.salt75_line,3),size(NIFS6.salt75_line,2));
        NIFS6.salt75_line_comb(numi,:)=NIFS6.salt75_line(i,:,j);
        NIFS6.salt_surf_line_comb(numi,:)=NIFS6.salt_surf_line(i,:,j);
        numi=numi+1;
    end
end
NIFS6.salt75_line_comb(NIFS6.salt75_line_comb==0)=NaN;

% plot(squeeze(NIFS6.salt75_line_comb(:,3)))  %%206 line

for numline=1:size(NIFS6.salt75_line,2)
    ttt=1/6:1/6:length(NIFS6.salt75_line_comb(:,numline))/6;
    idx=isnan(NIFS6.salt75_line_comb(:,numline)');
    p=polyfit(ttt(~idx),squeeze(NIFS6.salt75_line_comb(~idx,numline))',1);
    NIFS6.salt75_line_comb_trend(numline)=p(1);
    idx=isnan(NIFS6.salt_surf_line_comb(:,numline)');
    p=polyfit(ttt(~idx),squeeze(NIFS6.salt_surf_line_comb(~idx,numline))',1);
    NIFS6.salt_surf_line_comb_trend(numline)=p(1);
end
NIFS6.salt75_line_comb_trend(NIFS6.salt75_line_comb_trend==0)=NaN;
NIFS6.salt_surf_line_comb_trend(NIFS6.salt_surf_line_comb_trend==0)=NaN;

% % % [year line sta mon dep]
for i=1:size(NIFS6.salt,2)  %% line
    for j=1:size(NIFS6.salt,3)  %% sta
        numi=1;
        for k=1:size(NIFS6.salt,1) %% year
            for l=1:size(NIFS6.salt,4) %% month
        %         NIFS6.salt75_line_comb=reshape(NIFS6.salt75_line,size(NIFS6.salt75_line,1)*size(NIFS6.salt75_line,3),size(NIFS6.salt75_line,2));
                NIFS6.salt_surf_comb(numi,i,j)=NIFS6.salt(k,i,j,l,1);
                numi=numi+1;
            end
        end
    end
end
NIFS6.salt_surf_comb(NIFS6.salt_surf_comb==0)=NaN;

% % % [year line sta mon dep]
for i=1:size(NIFS6.salt75_sta,2)  %% line
    for j=1:size(NIFS6.salt75_sta,3)  %% sta
        numi=1;
        for k=1:size(NIFS6.salt75_sta,1) %% year
            for l=1:size(NIFS6.salt75_sta,4) %% month
                NIFS6.salt75_comb(numi,i,j)=NIFS6.salt75_sta(k,i,j,l);
                numi=numi+1;
            end
        end
    end
end
NIFS6.salt75_comb(NIFS6.salt75_comb==0)=NaN;


for i=1:size(NIFS6.salt,2)  %% line
    for j=1:size(NIFS6.salt,3)  %% sta
        idx=isnan(NIFS6.salt_surf_comb(:,i,j)');
        p=polyfit(ttt(~idx),squeeze(NIFS6.salt_surf_comb(~idx,i,j))',1);
        NIFS6.salt_surf_comb_trend(i,j)=p(1);
    end
end
NIFS6.salt_surf_comb_trend(NIFS6.salt_surf_comb_trend==0)=NaN;

for i=1:size(NIFS6.salt75_comb,2)  %% line
    for j=1:size(NIFS6.salt75_comb,3)  %% sta
        idx=isnan(NIFS6.salt75_comb(:,i,j)');
        p=polyfit(ttt(~idx),squeeze(NIFS6.salt75_comb(~idx,i,j))',1);
        NIFS6.salt75_comb_trend(i,j)=p(1);
    end
end
NIFS6.salt75_comb_trend(NIFS6.salt75_comb_trend==0)=NaN;


save(['E:\Data\Observation\NIFS\xls\','eco_NIFS_salt_trend_',num2str(min(year),'%04i'),'_',num2str(max(year),'%04i'),'.mat'], 'NIFS6');

