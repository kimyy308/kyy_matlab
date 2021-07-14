

rname ='D:\need_to_presentation\WOA2013/T_1_7.nc';
temp2_woa=ncread(rname,'T_MOD');
lon2_woa=ncread(rname,'LON1240_1240');
lat2_woa=ncread(rname,'LAT500_516');
level2_woa=ncread(rname,'DEPTH1_25');
ilevel2_woa = -level2_woa/20.0;

rname ='D:\need_to_presentation\WOA2013/T_1_7_FNR.nc';
temp2_woa_fnr=ncread(rname,'T_MOD_FNR');
% lon2_woa=ncread(rname,'LON1240_1240');
% lat2_woa=ncread(rname,'LAT500_516');
% level2_woa=ncread(rname,'DEPTH1_25');
% ilevel2_woa = -level2_woa/20.0;

for i=1:12
    temp2_woa_fnr2(:,:,i)=griddata(lat2_woa,ilevel2_woa,squeeze(temp2_woa_fnr(1,:,:,i))',(34.875:0.025:38.875)',(0:-0.125:-10));  %% [81 161 12]
%     temp2_woa_fnr2(:,:,i)=griddata(ilevel2_woa,lat2_woa,squeeze(temp2_woa_fnr(1,:,:,i)),(0:-0.125:-10),(34.875:0.025:38.875)');
end
lat2_woa_fnr=(34.875:0.025:38.875)';
ilevel2_woa_fnr=(0:-0.125:-10);

% 38.125 ~ 38.875 = NaN (4 point)
% -10~-8.75 =NaN (38.125)  -> 131~~, 81 ~71 -8.75
% -10~-7.5 = NaN (38.375)  -> 141~~, 81 ~61 -7.5
% -10~-6.25 = NaN (38.625) -> 151~~, 81 ~51 -6.25
% -10~-6.25 = NaN (38.875) -> 161~~, 81 ~51 -6.25

% % % make interpolated land
% % for i=31:-1:11
% %     temp2_woa_fnr2(81:-1:71-(-i+31),i,:) = NaN;
% % end
% % for i=10:-1:1
% %     temp2_woa_fnr2(81:-1:51,i,:) = NaN;
% % end

yfilled(1:10)=-6.25;
for i=11:31
    yfilled(i)=-6.25-(i-10)*0.125;
end
yfilled(32)=-10;
yfilled(33)=-10;
yfilled(34)=-6.25;
% for i=1:31
%     yfilled2(i)=-6.25-(i-1)*(2.5/30.0)
% end
xfilled(1:31)=lat2_woa_fnr(1:31);
xfilled(32)=lat2_woa_fnr(31);
xfilled(33)=lat2_woa_fnr(1);
xfilled(34)=lat2_woa_fnr(1);
'reading data for figure_woa is completed'


