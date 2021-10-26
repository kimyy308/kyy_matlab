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
    case('ES') %% East Sea
        refpolygon=espolygon;
    case('NES') %% Northern East Sea
        refpolygon=nespolygon;
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
    case('CA') %% Coastal Area around korea peninsula
        refpolygon=capolygon;
    case('EKB') %% Coastal Area around korea peninsula
        refpolygon=ekbpolygon;
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
    otherwise
        ('?')
end
lonlat(1)=min(refpolygon(:,1));
lonlat(2)=max(refpolygon(:,1));
lonlat(3)=min(refpolygon(:,2));
lonlat(4)=max(refpolygon(:,2));

error_status=1;
end

