% % % making isodepth_3dd.nc

% % % expname='restore_f_plane_07_01';
% % % 
% % % readname='D:\need_to_presentation\';
% % % readname2='\197912\ocean_snap_197912.nc';
% % % rname=strcat(readname,expname,readname2);
% % % lon2=ncread(rname,'xt_ocean');
% % % lat2=ncread(rname,'yt_ocean');
% % % level2=ncread(rname,'st_ocean');
% % % ilevel2 = -level2;
% % % temp=ncread(rname,'temp');
% % % temp=squeeze(mean(temp,4));
% % % [lon, lat] = meshgrid(lon2',lat2);
% % % 
% % % % for 19.5 C isosurface. put NaN value in the maximum depth
% % % for i=1:120
% % %     for j=1:110
% % %         k=1;
% % %         while(k>0)
% % %             if (temp(i,j,k)>19.5);
% % %                 k=k+1;
% % %             elseif (temp(i,j,1)<0)
% % %                 isodepth(i,j)=NaN;
% % %                 k=0;
% % %             else
% % %                 isodepth(i,j)=(k-1)*10;
% % %                 if (temp(i,j,k)<0)
% % %                     isodepth(i,j)=NaN;
% % %                 end
% % %                 k=0;
% % %             end
% % %         end
% % %     end
% % %     i
% % % end
% % % 
% % % % for 19.5 C isosurface, + put maximum depth value
% % % for i=1:120
% % %     for j=1:110
% % %         k=1;
% % %         while(k>0)
% % %             if (temp(i,j,k)>19.5);
% % %                 k=k+1;
% % %             elseif (temp(i,j,1)<0)
% % %                 isodepth2(i,j)=NaN;
% % %                 k=0;
% % %             else
% % %                 isodepth2(i,j)=(k-1)*10;
% % %                 k=0;
% % %             end
% % %         end
% % %     end
% % %     i
% % % end
% % % 
% % % % for isodepth, 10.0 C isosurface.
% % % for i=1:120
% % %     for j=1:110
% % %         k=1;
% % %         while(k>0)
% % %             if (temp(i,j,k)>10.0);
% % %                 k=k+1;
% % %             elseif (temp(i,j,1)<0)
% % %                 isodepth3(i,j)=NaN;
% % %                 k=0;
% % %             else
% % %                 isodepth3(i,j)=(k-1)*10;
% % %                 if (temp(i,j,k)<0)
% % %                     isodepth3(i,j)=NaN;
% % %                 end
% % %                 k=0;
% % %             end
% % %         end
% % %     end
% % %     i
% % % end
% % % 
% % % % for 10.0 C isosurface, + put maximum depth value
% % % for i=1:120
% % %     for j=1:110
% % %         k=1;
% % %         while(k>0)
% % %             if (temp(i,j,k)>10.0);
% % %                 k=k+1;
% % %             elseif (temp(i,j,1)<0)
% % %                 isodepth4(i,j)=NaN;
% % %                 k=0;
% % %             else
% % %                 isodepth4(i,j)=(k-1)*10;
% % %                 k=0;
% % %             end
% % %         end
% % %     end
% % %     i
% % % end
% % % 
% % % ncid= netcdf.create('D:\need_to_presentation\restore_f_plane_07_01\isodepth_3dd.nc','NOCLOBBER');
% % % xdimid = netcdf.defDim(ncid,'XT_OCEAN',120);
% % % ydimid = netcdf.defDim(ncid,'YT_OCEAN',110);
% % % xvarid = netcdf.defVar(ncid,'XT_OCEAN','double',xdimid);
% % % yvarid = netcdf.defVar(ncid,'YT_OCEAN','double',ydimid);
% % % isovarid = netcdf.defVar(ncid,'isodepth','double',[xdimid,ydimid]);
% % % iso2varid = netcdf.defVar(ncid,'isodepth2','double',[xdimid,ydimid]);
% % % iso3varid = netcdf.defVar(ncid,'isodepth3','double',[xdimid,ydimid]);
% % % iso4varid = netcdf.defVar(ncid,'isodepth4','double',[xdimid,ydimid]);
% % % netcdf.putAtt(ncid,isovarid,'info','19.5^o isosurface, + NaN values are used in the warm area');
% % % netcdf.putAtt(ncid,iso2varid,'info','19.5^o isosurface, + Maximum depth values are used in the warm area');
% % % netcdf.putAtt(ncid,iso3varid,'info','10.0^o isosurface, + NaN values are used in the warm are');
% % % netcdf.putAtt(ncid,iso4varid,'info','10.0^o isosurface, + Maximum depth values are used in the warm area');
% % % netcdf.endDef(ncid);
% % % netcdf.putVar(ncid,xvarid,lon2);
% % % netcdf.putVar(ncid,yvarid,lat2);
% % % netcdf.putVar(ncid,isovarid,isodepth);
% % % netcdf.putVar(ncid,iso2varid,isodepth2);
% % % netcdf.putVar(ncid,iso3varid,isodepth3);
% % % netcdf.putVar(ncid,iso4varid,isodepth4);
% % % netcdf.close(ncid);


bwr_map;
expname='restore_f_plane_07_01';
readname='D:\need_to_presentation\';
readname2='\isodepth_3dd.nc';
rname=strcat(readname,expname,readname2);
lon2=ncread(rname,'XT_OCEAN');
lat2=ncread(rname,'YT_OCEAN');
isod1=ncread(rname,'isodepth');
isod2=ncread(rname,'isodepth2');
isod3=ncread(rname,'isodepth3');
isod4=ncread(rname,'isodepth4');
isod1=(0-isod1);
isod2=(0-isod2);
isod3=(0-isod3);
isod4=(0-isod4);

for i=1:120
    for j=1:110
        if(isnan(isod1(i,j))==0)
            isod5(i,j)=-200.0;
        else
            isod5(i,j)=NaN;
        end
        if(isnan(isod3(i,j))==0)
            isod6(i,j)=-200.0;
        else
            isod6(i,j)=NaN;
        end
    end
end

for i=4:23
    for j=1:10
        if(isnan(isod5(i,j))==0)
            isod5(i,j)=-100;
            isod6(i,j)=-100;
        end
    end
end
for i=4:23
    for j=11:35
        if(isnan(isod5(i,j))==0)
            isod5(i,j)= -( 200 - 100*((35-j)/24.0) );
        end
        if(isnan(isod6(i,j))==0)
            isod6(i,j)= -( 200 - 100*((35-j)/24.0) );
        end
    end
end



'reading data for figure 9 is completed'