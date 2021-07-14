bwr_map;
expname='restore_f_plane_07_01';

readname='D:\need_to_presentation\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname,readname2);
u2=ncread(rname,'u');
v2=ncread(rname,'v');
temp2=ncread(rname,'temp');
lon2=ncread(rname,'xt_ocean');
lat2=ncread(rname,'yt_ocean');
level2=ncread(rname,'st_ocean');
ilevel2 = -level2;

readname='D:\need_to_presentation\';
readname2='\197912\ocean_gridinfo_197912.nc';
rname=strcat(readname,expname,readname2);
ht2=ncread(rname,'ht');
ht3=ht2;

[lon, lat] = meshgrid(lon2',lat2);

for i=1:120
    for j=1:110
        for k=1:30
            for l=1:12
                if (u2(i,j,k,l)< -1.0e+5)
                    u2(i,j,k,l)=NaN;
                end
                if (v2(i,j,k,l)< -1.0e+5)
                    v2(i,j,k,l)=NaN;
                end
                if (temp2(i,j,k,l)< -1.0e+5)
                    temp2(i,j,k,l)=NaN;
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
    u2(82:84,14:17,1:30,i)=0.2;
    v2(82:84,14:17,1:30,i)=0.00001;
end

for i=1:120
  for j=1:110
    for k=1:30
      u2mean(i,j,k,1)=mean(u2(i,j,k,1:12));
      v2mean(i,j,k,1)=mean(v2(i,j,k,1:12));
      temp2mean(i,j,k,1)=mean(temp2(i,j,k,1:12));
    end
  end
end
'reading data for figure 7 is completed'