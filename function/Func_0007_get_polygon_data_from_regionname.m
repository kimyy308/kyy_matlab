function [refpolygon, lonlat, error_status] = Func_0007_get_polygon_data_from_regionname(regionname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [refpolygon, lonlat]=get_polygon_data_from_regionname(regionname);
%
% get the polygon data from nwp_polygon_point.m corresponding to regionname
%
%  input:
%  regionname             ROMS RCM region name (string)
%
%  output:
%  refpolygon             Reference polygon point (2-D array, [lon; lat])
%  lonlat                 edge of the domain (lon_min, lon_max, lat_min, lat_max)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    17-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('nwp_polygon_point.m');

switch(regionname)
    case('NWP') %% North western Pacific
        lonlat = [115, 164, 15, 52];  %% whole data area
        refpolygon(1,1)=lonlat(1);
        refpolygon(2,1)=lonlat(2);
        refpolygon(1,2)=lonlat(3);
        refpolygon(2,2)=lonlat(4);
    case('NWP2') %% North western Pacific
        lonlat = [115, 145, 25, 52];  %% whole data area
        refpolygon(1,1)=lonlat(1);
        refpolygon(2,1)=lonlat(2);
        refpolygon(1,2)=lonlat(3);
        refpolygon(2,2)=lonlat(4);
    case('SNWP') %% southern northwest Pacific
        refpolygon=snwppolygon;
    case('ES') %% East Sea
        refpolygon=espolygon;
    case('EKWC2') %% East Korean Warm Current
        refpolygon=ekwc2polygon;
    case('NES') %% Northern East Sea
        refpolygon=nespolygon;
    case('NES2') %% Northern East Sea
        refpolygon=nes2polygon;
    case('SES') %% Southern East Sea
        refpolygon=sespolygon;
    case('SS') %% South Sea
        refpolygon=sspolygon;
    case('YS') %% Yellow Sea
        refpolygon=yspolygon;
    case('ECS') %% East China Sea
        refpolygon=ecspolygon;
     case('ECS2') %% East China Sea
        refpolygon=ecs2polygon;
    case('YSECS') %% East China Sea
        refpolygon=ysecspolygon;
    case('AKP') %% Around Korea Peninsula
        refpolygon=akppolygon;
    case('AKP2') %% Around Korea Peninsula
        refpolygon=akp2polygon;
    case('AKP3') %% Around Korea Peninsula
        refpolygon=akp3polygon;
    case('AKP4') %% Around Korea Peninsula
        refpolygon=akp4polygon;
    case('CA') %% Around Korea Peninsula
        refpolygon=capolygon;
    case('EKB') %% Around Korea Peninsula
        refpolygon=akp2polygon;
    case('BOH') %% Around Korea Peninsula
        refpolygon=bohpolygon;
    case('pollock_egg')
        refpolygon=pollock_eggpolygon;
    case('pollock_egg2')
        refpolygon=pollock_egg2polygon;
    case('pollock_egg3')
        refpolygon=pollock_egg3polygon;
    case('pollock_egg4')
        refpolygon=pollock_egg4polygon;
    case('CA') %% Coastal Area around korea peninsula
        refpolygon=capolygon;
    case('EKB') %% Coastal Area around korea peninsula
        refpolygon=ekbpolygon;
    case('EKB2') %% Coastal Area around korea peninsula
        refpolygon=ekb2polygon;
    case('ES_KHOA') %% East Sea
        refpolygon=es_khoapolygon;
    case('YS_KHOA') %% East Sea
        refpolygon=ys_khoapolygon;
    case('SS_KHOA') %% East Sea
        refpolygon=ss_khoapolygon;
    case('TEST') %% for debugging
        refpolygon=testpolygon;
    case('SK_EEZ') %%south korea EEZ (real fishing area)
        refpolygon=SK_EEZ_polygon;
    case('SK_coastal') %%south korea coastal area
        refpolygon=SK_coastal_polygon;
    case('adm_div_all') %% south korea administrative ocean division
        refpolygon=adm_div_all_polygon;
    case('adm_div_YS') %% south korea administrative ocean division (Yellow Sea)
        refpolygon=adm_div_YS_polygon;
    case('adm_div_SS') %% south korea administrative ocean division (South Sea)
        refpolygon=adm_div_SS_polygon;
    case('adm_div_ES') %% south korea administrative ocean division (East Sea)
        refpolygon=adm_div_ES_polygon;
    case('adm_div_GGD') %% south korea administrative ocean division (Gyeonggi-do)
        refpolygon=adm_div_GGD_polygon;
    case('adm_div_CCND') %% south korea administrative ocean division (Chungcheongnam-do)
        refpolygon=adm_div_CCND_polygon;
    case('adm_div_JBD') %% south korea administrative ocean division (Jeollabuk-do)
        refpolygon=adm_div_JBD_polygon;
    case('adm_div_JND') %% south korea administrative ocean division (Jeollabuk-do)
        refpolygon=adm_div_JND_polygon;
    case('adm_div_JJD') %% south korea administrative ocean division (Jeju-do)
        refpolygon=adm_div_JJD_polygon;
    case('adm_div_GSND') %% south korea administrative ocean division (Gyeongsangnam-do)
        refpolygon=adm_div_GSND_polygon;
    case('adm_div_GSBD') %% south korea administrative ocean division (Gyeongsangbuk-do)
        refpolygon=adm_div_GSBD_polygon;
    case('adm_div_GSBD_coastal') %% south korea administrative ocean division (Gyeongsangbuk-do, ulleung-do and dok-do are excluded)
        refpolygon=adm_div_GSBD_coastal_polygon;
    case('adm_div_ULD') %% south korea administrative ocean division (Ulleung-do and dok-do)
        refpolygon=adm_div_ULD_polygon;
    case('adm_div_ULD_only') %% south korea administrative ocean division (Ulleung-do only)
        refpolygon=adm_div_ULD_only_polygon;
    case('adm_div_DD_only') %% south korea administrative ocean division (Dok-do only)
        refpolygon=adm_div_DD_only_polygon;
    case('adm_div_GWD') %% south korea administrative ocean division (Ulleung-do and dok-do)
        refpolygon=adm_div_GWD_polygon;
    case('ref_sal') %% south korea administrative ocean division (Ulleung-do and dok-do)
        refpolygon=ref_sal_polygon;
    otherwise
        ('?')
end
lonlat(1)=min(refpolygon(:,1));
lonlat(2)=max(refpolygon(:,1));
lonlat(3)=min(refpolygon(:,2));
lonlat(4)=max(refpolygon(:,2));

error_status=1;
end

