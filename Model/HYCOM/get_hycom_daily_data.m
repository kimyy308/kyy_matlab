% % 2021-01-05 yykim
clc;clear all;close all
year = 2002;
filepath = ['F:\hycom\',num2str(year),'\'];
if exist(filepath)==0
    mkdir(filepath)
end
cd(filepath);
list=datenum(year,01,02):datenum(year,12,31);
op=weboptions('timeout',14400);
for i=1:length(list)
    version='expt_53.X';
    tim = [00 03 06 09 12 15 18 21];
%     tim = [00];
    tic;
    for k = 1 : length(tim)
        j = tim(k);
        [yy,mm,dd]=datevec(list(i));
        disp([sprintf('%04d',yy),'-',sprintf('%02d',mm),'-',sprintf('%02d',dd),'-',sprintf('%02d',j)]);
        websave([sprintf('NEAMS_%04d',yy),sprintf('%02d',mm),sprintf('%02d',dd),sprintf('%02d',j),'.nc4'],...
            ['http://ncss.hycom.org/thredds/ncss/GLBv0.08/',sprintf('%s',version),'/data/',num2str(yy),'?var=salinity&var=water_temp&var=water_u&var=water_v&north=55&west=125&east=150&south=30&disableProjSubset=on&horizStride=1&time=',num2str(yy),'-',sprintf('%02d',mm),'-',sprintf('%02d',dd),'T',sprintf('%02d',j),'%3A00%3A00Z&vertStride=1&addLatLon=true&accept=netcdf4'],op);
%              http://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_53.X/data/2015?var=salinity_bottom&var=surf_el&var=water_temp_bottom&var=water_u_bottom&var=water_v_bottom&var=salinity&var=water_temp&var=water_u&var=water_v&disableLLSubset=on&disableProjSubset=on&horizStride=1&time_start=2015-01-01T12%3A00%3A00Z&time_end=2015-12-31T09%3A00%3A00Z&timeStride=1&vertCoord=&accept=netcdf
    end
    toc;
end