expname='oman_restore_04_18';

readname='D:\need_to_presentation\';
readname2='\ke.nc';
rname=strcat(readname,expname,readname2);
kine=ncread(rname,'KINE');

readname='D:\need_to_presentation\';
readname2='\te.nc';
rname=strcat(readname,expname,readname2);
tempe=ncread(rname,'TEMPE');
'reading data for figure 4 is completed'