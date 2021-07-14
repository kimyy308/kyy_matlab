% % read upperlayer thickness, vorticity, coriolis frequency

bwr_map;


readname='D:\need_to_presentation\';
readname2='\isodepth.nc';
rname=strcat(readname,expname1,readname2);

isod5_1=ncread(rname,'isodepth');
vort5_1=ncread(rname,'sfc_rel_vort');
f5_1=ncread(rname,'f');
pot_vort5_1=ncread(rname,'pot_vort');
lon5_1=ncread(rname,'XT_OCEAN');
lat5_1=ncread(rname,'YT_OCEAN');

vort5mean_1=mean(vort5_1(:,:,109:120),3);
isod5mean_1=mean(isod5_1(:,:,109:120),3);
pot_vort5mean_1=mean(pot_vort5_1(:,:,109:120),3);


rname=strcat(readname,expname2,readname2);

isod5_2=ncread(rname,'isodepth');
vort5_2=ncread(rname,'sfc_rel_vort');
f5_2=ncread(rname,'f');
pot_vort5_2=ncread(rname,'pot_vort');
lon5_2=ncread(rname,'XT_OCEAN');
lat5_2=ncread(rname,'YT_OCEAN');
vort5mean_2=mean(vort5_2(:,:,109:120),3);
isod5mean_2=mean(isod5_2(:,:,109:120),3);
pot_vort5mean_2=mean(pot_vort5_2(:,:,109:120),3);


rname=strcat(readname,expname3,readname2);

isod5_3=ncread(rname,'isodepth');
vort5_3=ncread(rname,'sfc_rel_vort');
f5_3=ncread(rname,'f');
pot_vort5_3=ncread(rname,'pot_vort');
lon5_3=ncread(rname,'XT_OCEAN');
lat5_3=ncread(rname,'YT_OCEAN');
vort5mean_3=mean(vort5_3(:,:,109:120),3);
isod5mean_3=mean(isod5_3(:,:,109:120),3);
pot_vort5mean_3=mean(pot_vort5_3(:,:,109:120),3);


rname=strcat(readname,expname4,readname2);

isod5_4=ncread(rname,'isodepth');
vort5_4=ncread(rname,'sfc_rel_vort');
f5_4=ncread(rname,'f');
pot_vort5_4=ncread(rname,'pot_vort');
lon5=ncread(rname,'XT_OCEAN');
lat5=ncread(rname,'YT_OCEAN');
vort5mean_4=mean(vort5_4(:,:,109:120),3);
isod5mean_4=mean(isod5_4(:,:,109:120),3);
pot_vort5mean_4=mean(pot_vort5_4(:,:,109:120),3);


'reading data for figure 19_9 is completed'