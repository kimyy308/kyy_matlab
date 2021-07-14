clear all;close all;clc;

brytype='SODA';
for yy=2007
    ydays=year2day(yy);

    if ydays == 366
        time_month=[31 29 31 30 31 30 31 31 30 31 30 31];
    else
        time_month=[31 28 31 30 31 30 31 31 30 31 30 31];
    end


    for i_ot=1:1:12
        if i_ot ==1
            time(i_ot,1)=(time_month(i_ot)/2);
        else
            time(i_ot,1)=(sum(time_month(1:i_ot-1))+time_month(i_ot)/2);
        end
    end
    time=[15:30:345];
    bryname=['roms_NP_bry2_',brytype,'-Y',num2str(yy),'.nc'];
%     grdname='D:\add2_ini_bry_grd\grid\roms_grid2_ADD_08_2_ep.nc';
    grdname=['D:\data\roms_input\np\grid\NP_grid_30layer.nc'];
create_bryfile_t(bryname,brytype,grdname,yy,ydays,time,1);
end


