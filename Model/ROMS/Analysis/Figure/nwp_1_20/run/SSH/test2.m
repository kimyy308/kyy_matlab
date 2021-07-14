

%comb_temp_raw_data; % 306 x 240 x 50 x 300

clear all; close all;

regionname = 'AKP2'
filedir = '/data1/stlee/ext_hdd/MyOcean/';
load([filedir,regionname,'_New_Myocean_comb_raw.mat']);
%% interpolation ÀÛ¾÷
% ¼ö½É ÀÚ·á¿¡ Ç¥Ãþ(0m)Ãß°¡
MyOcean_st3=[0;MyOcean_st];
v5=zeros(length(MyOcean_lon),length(MyOcean_lat),300);% vertical mean of temp
v6=zeros(300,1); % total spatial mean
end_depth=zeros(length(MyOcean_lon),length(MyOcean_lat),300);
for cc=1:300
     disp([num2str(cc)])
    num_v=zeros(length(MyOcean_lon),length(MyOcean_lat));
    v4=zeros(length(MyOcean_lon),length(MyOcean_lat));
    for aa=1:length(MyOcean_lon)
        for bb=1:length(MyOcean_lat)
            %Æ¯Á¤Áö¿ª ¼ö¿Â»Ì¾Æ¼­ Ç¥Ãþ(0m)¼ö¿Â Ãß°¡
            point_temp = comb_temp_raw_data(aa,bb,:,cc);
            point_temp2 = squeeze(point_temp);%1(Ç¥Ãþ)~50(½ÉÇØ)
            if (~isnan(point_temp2(1))) % NaNÀÌ ¾Æ´Ï¶ó¸é
                point_temp3 = zeros(51,1);
                point_temp3(1) = point_temp2(1);
                point_temp3(2:length(point_temp2)+1)=point_temp2(1:end);
                %
                a = find(isnan(point_temp3));
                b = find(~isnan(point_temp3)); % NaNÀÌ ¾Æ´Ñ°÷ Ã£±â
                std_dep = 0:1:fix(MyOcean_st3(b(end)));
                std_dep = std_dep';
                point_temp4 = interp1(MyOcean_st3(1:b(end)),point_temp3(1:b(end)),std_dep,'spline');
                v4(aa,bb) = sum(point_temp4); % sum of vertical temp
                num_v(aa,bb) = length(point_temp4); % ÇÑ pointÀÇ ºÐÇÒ°³¼ö
                end_depth(aa,bb,cc)=MyOcean_st3(b(end));
                clear a b point_temp point_temp2 std_dep point_temp3 point_temp4
            else % NaNÀÌ¶ó¸é
                v4(aa,bb)=NaN; % sum of vertical temp
                num_v(aa,bb)=NaN; % ÇÑ pointÀÇ ºÐÇÒ°³¼ö
                end_depth(aa,bb,cc) = 0;
            end
        end
    end
    v5(:,:,cc) = v4./num_v; % vertical mean of temp
    v6(cc) = sum(nansum(v4))/sum(nansum(num_v)); % total spatial mean
    clear v4 num_v
end

  save([filedir,regionname,'akp2_new_interpol_data','.mat'],'v5','v6','varname','regionname',...
      'MyOcean_lon','MyOcean_lat','MyOcean_st3','end_depth');


