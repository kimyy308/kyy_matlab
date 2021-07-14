bwr_map;
% % % 
% % % exp 1.
% % % 


readname='D:\need_to_presentation\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname1,readname2);
u2_1=ncread(rname,'u');
v2_1=ncread(rname,'v');

temp2_1=ncread(rname,'temp');
lon2=ncread(rname,'xt_ocean');
lat2=ncread(rname,'yt_ocean');
level2=ncread(rname,'st_ocean');
ilevel2 = -level2;

readname2='\197912\ocean_diag_197912.nc';
rname=strcat(readname,expname1,readname2);
w2_1=ncread(rname,'wt');

'reading data for figure 17_exp1 is completed'

readname='D:\need_to_presentation\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname2,readname2);
u2_2=ncread(rname,'u');
v2_2=ncread(rname,'v');
temp2_2=ncread(rname,'temp');

readname2='\197912\ocean_diag_197912.nc';
rname=strcat(readname,expname1,readname2);
w2_2=ncread(rname,'wt');

'reading data for figure 17_exp2 is completed'

readname='D:\need_to_presentation\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname3,readname2);
u2_3=ncread(rname,'u');
v2_3=ncread(rname,'v');
temp2_3=ncread(rname,'temp');

readname2='\197912\ocean_diag_197912.nc';
rname=strcat(readname,expname1,readname2);
w2_3=ncread(rname,'wt');

'reading data for figure 17_exp3 is completed'

readname='D:\need_to_presentation\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname4,readname2);
u2_4=ncread(rname,'u');
v2_4=ncread(rname,'v');
temp2_4=ncread(rname,'temp');

readname2='\197912\ocean_diag_197912.nc';
rname=strcat(readname,expname1,readname2);
w2_4=ncread(rname,'wt');

'reading data for figure 17_exp4 is completed'

readname='D:\need_to_presentation\';
readname2='\197912\ocean_gridinfo_197912.nc';
rname=strcat(readname,expname1,readname2);
ht2=ncread(rname,'ht');
ht3=ht2;

[lon, lat] = meshgrid(lon2',lat2);

for i=1:120
    for j=1:110
        for k=1:30
            for l=1:12
                if (u2_1(i,j,k,l)< -1.0e+5)
                    u2_1(i,j,k,l)=NaN;
                    u2_2(i,j,k,l)=NaN;
                    u2_3(i,j,k,l)=NaN;
                    u2_4(i,j,k,l)=NaN;
                end
                if (v2_1(i,j,k,l)< -1.0e+5)
                    v2_1(i,j,k,l)=NaN;
                    v2_2(i,j,k,l)=NaN;
                    v2_3(i,j,k,l)=NaN;
                    v2_4(i,j,k,l)=NaN;
                end
                if (w2_1(i,j,k,l)< -1.0e+5)
                    w2_1(i,j,k,l)=NaN;
                    w2_2(i,j,k,l)=NaN;
                    w2_3(i,j,k,l)=NaN;
                    w2_4(i,j,k,l)=NaN;
                end
                if (temp2_1(i,j,k,l)< -1.0e+5)
                    temp2_1(i,j,k,l)=NaN;
                    temp2_2(i,j,k,l)=NaN;
                    temp2_3(i,j,k,l)=NaN;
                    temp2_4(i,j,k,l)=NaN;
                end
            end
        end
        if (ht3(i,j)< -1.0e+5)
          ht3(i,j)=NaN;
        else
          ht3(i,j)=0;
        end
    end
end


for i=1:12
    u2_1(44:46,20:22,1:30,i)=0.2;
    v2_1(44:46,20:22,1:30,i)=0.00001;
    u2_2(44:46,20:22,1:30,i)=0.2;
    v2_2(44:46,20:22,1:30,i)=0.00001;
    u2_3(44:46,20:22,1:30,i)=0.2;
    v2_3(44:46,20:22,1:30,i)=0.00001;
    u2_4(44:46,20:22,1:30,i)=0.2;
    v2_4(44:46,20:22,1:30,i)=0.00001;
    v2_1(1:20,5:9,16,i)=0.2;
    w2_1(1:20,5:9,16,i)=0.00001;
    v2_2(1:20,5:9,16,i)=0.2;
    w2_2(1:20,5:9,16,i)=0.00001;
    v2_3(1:20,5:9,16,i)=0.2;
    w2_3(1:20,5:9,16,i)=0.00001;
    v2_4(1:20,5:9,16,i)=0.2;
    w2_5(1:20,5:9,16,i)=0.00001;
end


for i=1:120
  for j=1:110
    for k=1:30
      u2mean_1(i,j,k,1)=mean(u2_1(i,j,k,1:12));
      v2mean_1(i,j,k,1)=mean(v2_1(i,j,k,1:12));
      w2mean_1(i,j,k,1)=mean(w2_1(i,j,k,1:12));
      u2mean_2(i,j,k,1)=mean(u2_2(i,j,k,1:12));
      v2mean_2(i,j,k,1)=mean(v2_2(i,j,k,1:12));
      w2mean_2(i,j,k,1)=mean(w2_1(i,j,k,1:12));
      u2mean_3(i,j,k,1)=mean(u2_3(i,j,k,1:12));
      v2mean_3(i,j,k,1)=mean(v2_3(i,j,k,1:12));
      w2mean_3(i,j,k,1)=mean(w2_1(i,j,k,1:12));
      u2mean_4(i,j,k,1)=mean(u2_4(i,j,k,1:12));
      v2mean_4(i,j,k,1)=mean(v2_4(i,j,k,1:12));
      w2mean_4(i,j,k,1)=mean(w2_1(i,j,k,1:12));
      temp2mean_1(i,j,k,1)=mean(temp2_1(i,j,k,1:12));
      temp2mean_2(i,j,k,1)=mean(temp2_2(i,j,k,1:12));
      temp2mean_3(i,j,k,1)=mean(temp2_3(i,j,k,1:12));
      temp2mean_4(i,j,k,1)=mean(temp2_4(i,j,k,1:12));
    end
  end
end

% maxekwc_1=max(max(v2mean_1(1:15,21:80,1,1)))
% maxekwc_2=max(max(v2mean_2(1:15,21:80,1,1)))
% maxekwc_3=max(max(v2mean_3(1:15,21:80,1,1)))
% maxekwc_4=max(max(v2mean_4(1:15,21:80,1,1)))
% [i1, i2, i3, i4] = ind2sub(size(v2mean_1), find(v2mean_1==maxekwc_1));
% [i1, i2, i3, i4] = ind2sub(size(v2mean_2), find(v2mean_2==maxekwc_2));
% [i1, i2, i3, i4] = ind2sub(size(v2mean_3), find(v2mean_3==maxekwc_3));
% [i1, i2, i3, i4] = ind2sub(size(v2mean_4), find(v2mean_4==maxekwc_4));
'reading data for figure 21 is completed'


