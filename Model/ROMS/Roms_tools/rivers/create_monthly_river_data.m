close all; clc; clear all;


startyear=1961;
endyear=2019;
filedir=['Z:\내 드라이브\Data\Observation\CJD'];
fid=fopen([filedir, '\datong_1961_2018.txt']);
data=textscan(fid, '%f %f %f', 'headerlines', 1);
fclose(fid);
oldCJD=data{2};  % 1961 ~ 2018 CJD

fout= [filedir, '\', 'datong_1961_', num2str(endyear), '.dat' ];
fid=fopen(fout);
for yeari=startyear:2018
    for monthi=1:12
        time1((yeari-1961)*12+monthi)=yeari+monthi/100;
        time2((yeari-1961)*12+monthi)=yeari+(monthi-1)/12;
    end
end

time1(end+1)=2019.01;
time1(end+1)=2019.02;
time1(end+1)=2019.03;
time1(end+1)=2019.04;
time1(end+1)=2019.05;

time2(end+1)=2019.00;
time2(end+1)=2019.08;
time2(end+1)=2019.17;
time2(end+1)=2019.25;
time2(end+1)=2019.33;

oldCJD(end+1)= 20101; %Jan
oldCJD(end+1)= 19200; %Feb
oldCJD(end+1)= 30100; %Mar
oldCJD(end+1)= 26400; %Apr
oldCJD(end+1)= 37600; %May

yeari=2019;
for monthi=6:12
    monthi
    mfile =[filedir, '\', num2str(endyear), '\', 'datong_archive_', num2str(yeari), num2str(monthi, '%02i'), '.dat'];
    fid=fopen(mfile);
    data=textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f');
    if monthi==6
        year = data{1}';
        month = data{2}';
        day = data{3}';
        cjd = data{7}';
    else
        year = [year data{1}'];
        month = [month data{2}'];
        day = [day data{3}'];
        cjd = [cjd data{7}'];
    end
    fclose(fid);
end
for monthi=6:12
    mind=find(month==monthi);
    tempcjd=cjd(mind);
    tempcjd(tempcjd==0)=NaN;
    oldCJD(end+1)=mean(tempcjd, 'omitnan');
    time1(end+1)=time1(end)+0.01;
    time2(end+1)=2019+(monthi-1)/12;
end

fid=fopen([filedir, '\datong_1961_', num2str(endyear), '.txt'], 'w+');
for yeari=startyear:endyear
    for monthi=1:12
        fprintf(fid, '%6.2f %11.2f %11.2f\n',[time1((yeari-1961)*12+monthi), oldCJD((yeari-1961)*12+monthi), time2((yeari-1961)*12+monthi)]);
    end
end
fclose(fid);

