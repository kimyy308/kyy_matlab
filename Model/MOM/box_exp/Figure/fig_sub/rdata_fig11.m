bwr_map;


readname='D:\need_to_presentation\';
readname2='\isodepth.nc';
rname=strcat(readname,expname,readname2);

isod5=ncread(rname,'isodepth');
vort5=ncread(rname,'sfc_rel_vort');
f5=ncread(rname,'f');
pot_vort5=ncread(rname,'pot_vort');
lon5=ncread(rname,'XT_OCEAN');
lat5=ncread(rname,'YT_OCEAN');



vort5mean=mean(vort5(:,:,109:120),3);
isod5mean=mean(isod5(:,:,109:120),3);
pot_vort5mean=mean(pot_vort5(:,:,109:120),3);


'reading data for figure 11 is completed'