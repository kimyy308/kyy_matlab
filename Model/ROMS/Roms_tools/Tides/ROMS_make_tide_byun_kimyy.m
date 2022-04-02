clear all; close all;
%% updated 31st March, 2022, Y.-Y. Kim

dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
ele=1;uv=1;
gridname='D:\Data\Model\ROMS\nwp_1_20\input\test2117\roms_grid_nwp_1_20_test2117.nc';
%source='NaoJ';
source='TPXO7';
TPXO7name=['D:\Data\Model\ROMS\nwp_1_20\etc\Roms_tools\TPXO7\TPXO7.nc'];
% tidenum=10;
tidenum=2; %M2, S2
outputname=['D:\Data\Model\ROMS\nwp_1_20\input\tides\nwp_1_20_tides_byun_m2s2.nc'];

%  source='TPXO7.2';
% tidenum=16;

% source='large domain';
adj_depth='n';

%%

% % % % % % % [infile,inpath]=uigetfile('*.nc','Grid File');infile=char(infile);
% % nc=netcdf([inpath,infile],'write');

nc=netcdf([gridname],'write');
lon_rho=nc{'lon_rho'}(:);
lat_rho=nc{'lat_rho'}(:);
lon_u=nc{'lon_u'}(:);
lat_u=nc{'lat_u'}(:);
lon_v=nc{'lon_v'}(:);
lat_v=nc{'lat_v'}(:);
% angle=nc{'angle'}(:);
h=nc{'h'}(:);
mask_rho=nc{'mask_rho'}(:);
close(nc)
[N M]=size(h);
disp(['Source = ',source])
switch source
    case 'NaoJ'
        if tidenum>16;tidenum=16;end
        comp = {'M2';'S2';'K1';'O1';'N2';'K2';'P1';'Q1';'J1';'OO1';'M1';'2N2';'MU2';'T2';'NU2';'L2'};
        components = 'M2  S2  K1  O1  N2  K2  P1  Q1  J1  OO1  M1  2N2  MU2  T2  NU2  T2';
        period=[12.4206012; 12.0000000; 23.9344697; 25.8193417; 12.6583482; 11.9672348; 24.0658902; 26.8683567;...
        23.0984768; 22.3060742; 24.8412024; 12.9053745; 12.8717576; 12.0164492; 12.6260044; 12.1916202];
        for ij=1:tidenum
            if ele==1
                tide_data = load(['./Nao99Jb/tide/',char(comp(ij,:)),'_j_out.dat']);
                lon = tide_data(:,1);
                lat = tide_data(:,2);
                amp = tide_data(:,3)/100;
                pha = tide_data(:,4);
                index = find(lon >= min(min(lon_rho))-1 & lon <= max(max(lon_rho)+1) & lat >= min(min(lat_rho))-1 & lat <= max(max(lat_rho))+1);
                lon_1 = lon(index); clear lon
                lat_1 = lat(index); clear lat
                amp_1 = amp(index); clear amp
                pha_1 = pha(index); clear pha
                pha_2 = pha_1;
                chan=find(pha_2>180 & pha_2<=360);
                pha_2(chan)=pha_2(chan)-360;
                zz_amp = griddata(lat_1,lon_1,amp_1,lat_rho,lon_rho);
                zz_pha = griddata(lat_1,lon_1,pha_1,lat_rho,lon_rho);
                zz_pha2 = griddata(lat_1,lon_1,pha_2,lat_rho,lon_rho);
                chan2=find(zz_pha2<0);
                zz_pha2(chan2) = zz_pha2(chan2)+360;
                convert = find((zz_pha>10&zz_pha<350)&((zz_pha2>0&zz_pha2<10)|(zz_pha2>350&zz_pha2<360)));
                zz_pha(convert) = zz_pha2(convert);
                tide_amp(ij,:,:) = inpaint_nans(zz_amp(:,:),1);
                tide_pha(ij,:,:) = inpaint_nans(mod(zz_pha(:,:),360),1);
                clear zz_amp zz_pha zz_pha2 chan chan2
            end
            if uv==1
                curr_data = load(['./Nao99Jb/curr/',char(comp(ij,:)),'_ellipse.dat']);
                lon = curr_data(:,1);
                lat = curr_data(:,2);
                Cmax = curr_data(:,3)/100;
                Cmin = curr_data(:,4)/100;
                Cangle = curr_data(:,5);
                Cphase = curr_data(:,6);
                lon_1 = lon(index); clear lon
                lat_1 = lat(index); clear lat
                Cmax_1 = Cmax(index); clear Cmax
                Cmin_1 = Cmin(index); clear Cmin
                Cangle_1 = Cangle(index); clear Cangle 
                Cphase_1 = Cphase(index); clear Cphase
                zz_Cmax = griddata(lat_1,lon_1,Cmax_1,lat_rho,lon_rho);
                zz_Cmin = griddata(lat_1,lon_1,Cmin_1,lat_rho,lon_rho);
                zz_Cangle = griddata(lat_1,lon_1,Cangle_1,lat_rho,lon_rho);
                zz_Cphase = griddata(lat_1,lon_1,Cphase_1,lat_rho,lon_rho);
                tide_Cmax(ij,:,:) = inpaint_nans(zz_Cmax(:,:),1);
                tide_Cmin(ij,:,:) = inpaint_nans(zz_Cmin(:,:),1);
                tide_Cangle(ij,:,:) = inpaint_nans(zz_Cangle(:,:),1);
                tide_Cphase(ij,:,:) = inpaint_nans(zz_Cphase(:,:),1);
                clear zz_Cmax zz_Cmin zz_Cangle zz_Cphase
            end
            disp([char(comp(ij,:)),' making complete !'])
        end
    case 'TPXO7'
        if tidenum>10;tidenum=10;end
%         comp=['m2';'s2';'n2';'k2';'k1';'o1';'p1';'q1';'mf';'mm'];
%         components = 'm2  s2  n2  k2  k1  o1  p1  q1  mf  mm';
%         period=[12.4206012; 12.0000000; 12.6583482; 11.9672348; 23.9344697; 25.8193417; 24.0658902; 26.8683567;...
%             327.8589689; 661.3092049];
        comp=['m2';'s2'];
        components = 'm2  s2';
        period=[12.4206012; 12.0000000];
        for ij=1:tidenum
            if ele==1
                ur=ext_data_tpxo(TPXO7name,'ssh_r',ij,lon_rho,lat_rho,'r',0);
                ui=ext_data_tpxo(TPXO7name,'ssh_i',ij,lon_rho,lat_rho,'r',0);
                ei=complex(ur,ui);
                tide_pha(ij,:,:)=mod(-180./pi.*angle(ei),360.0);
                tide_amp(ij,:,:)=abs(ei);
            end
            if uv==1
                ur=ext_data_tpxo(TPXO7name,'u_r',ij,lon_rho,lat_rho,'u',0);
                ui=ext_data_tpxo(TPXO7name,'u_i',ij,lon_rho,lat_rho,'u',0);
                ei=complex(ur,ui);
                upha=mod(-180./pi.*angle(ei),360.0); 
                uamp=abs(ei);
                ur=ext_data_tpxo(TPXO7name,'v_r',ij,lon_rho,lat_rho,'v',0);
                ui=ext_data_tpxo(TPXO7name,'v_i',ij,lon_rho,lat_rho,'v',0);
                ei=complex(ur,ui);
                vpha=mod(-180./pi.*angle(ei),360.0); 
                vamp=abs(ei);
                [major,eccentricity,inclination,phase]=ap2ep(uamp,upha,vamp,vpha);
                tide_Cmin(ij,:,:)=major.*eccentricity;
                tide_Cmax(ij,:,:)=major;
                tide_Cangle(ij,:,:)=inclination;
                tide_Cphase(ij,:,:)=phase;
            end
        end
    case 'TPXO7.2'
        if tidenum>13;tidenum=13;end
        comp={'M2 ';'S2 ';'N2 ';'K2 ';'K1 ';'O1 ';'P1 ';'Q1 ';'mf ';'mm ';'m4 ';'ms4';'mn4'};
        components = 'M2  S2  N2  K2  K1  O1  P1  Q1  mf  mm  m4  ms4 mn4';
        period=[12.4206012; 12.0000000; 12.6583482; 11.9672348; 23.9344697; 25.8193417; 24.0658902; 26.8683567;...
            327.8589689; 661.3092049; 6.2103006; 6.1033393; 6.2691739];
        ncload ./TPXO7.2/h_tpxo7.2.nc;
        ncload ./TPXO7.2/u_tpxo7.2.nc;
        for ij=1:tidenum
            if ele==1
                lon = reshape(lon_z,[],1);
                lat = reshape(lat_z,[],1);
                amp = reshape(squeeze(ha(ij,:,:)),[],1);%/100;
                pha = reshape(squeeze(hp(ij,:,:)),[],1);
                index = find(lon >= min(min(lon_rho))-1 & lon <= max(max(lon_rho)+1) & lat >= min(min(lat_rho))-1 & lat <= max(max(lat_rho))+1);
                lon_1 = lon(index); clear lon
                lat_1 = lat(index); clear lat
                amp_1 = amp(index); clear amp
                pha_1 = pha(index); clear pha
                pha_2 = pha_1;
                chan=find(pha_2>180 & pha_2<=360);
                pha_2(chan)=pha_2(chan)-360;
                zz_amp = griddata(lat_1,lon_1,amp_1,lat_rho,lon_rho);
                zz_pha = griddata(lat_1,lon_1,pha_1,lat_rho,lon_rho);
                zz_pha2 = griddata(lat_1,lon_1,pha_2,lat_rho,lon_rho);
                chan2=find(zz_pha2<0);
                zz_pha2(chan2) = zz_pha2(chan2)+360;
                convert = find((zz_pha>10&zz_pha<350)&((zz_pha2>0&zz_pha2<10)|(zz_pha2>350&zz_pha2<360)));
                zz_pha(convert) = zz_pha2(convert);
                tide_amp(ij,:,:) = inpaint_nans(zz_amp(:,:),1);
                tide_pha(ij,:,:) = inpaint_nans(zz_pha(:,:),1);
                clear zz_amp zz_pha zz_pha2 chan chan2 convert
            end
            if uv==1
                lon = reshape(lon_u,[],1);
                lat = reshape(lat_u,[],1);
                Cua = reshape(squeeze(ua(ij,:,:)),[],1);%/100;
                Cup = reshape(squeeze(up(ij,:,:)),[],1);
                lon_1 = lon(index); clear lon
                lat_1 = lat(index); clear lat
                Cua_1 = Cua(index); clear Cua
                Cup_1 = Cup(index); clear Cup
                Cup_2 = Cup_1;
                chan=find(Cup_2>180 & Cup_2<=360);
                Cup_2(chan)=Cup_2(chan)-360;
                zz_Cua = griddata(lat_1,lon_1,Cua_1,lat_rho,lon_rho);
                zz_Cup = griddata(lat_1,lon_1,Cup_1,lat_rho,lon_rho);
                zz_Cup2 = griddata(lat_1,lon_1,Cup_2,lat_rho,lon_rho);
                chan2=find(zz_Cup2<0);
                zz_Cup2(chan2) = zz_Cup2(chan2)+360;
                convert = find((zz_Cup>10&zz_Cup<350)&((zz_Cup2>0&zz_Cup2<10)|(zz_Cup2>350&zz_Cup2<360)));
                zz_Cup(convert) = zz_Cup2(convert);
                tide_Cua(ij,:,:) = zz_Cua(:,:);
                tide_Cup(ij,:,:) = zz_Cup(:,:);
                clear zz_Cua zz_Cup zz_Cup2 chan chan2 convert
            
                lon = reshape(lon_v,[],1);
                lat = reshape(lat_v,[],1);
                Cva = reshape(squeeze(va(ij,:,:)),[],1);%/100;
                Cvp = reshape(squeeze(vp(ij,:,:)),[],1);
                lon_1 = lon(index); clear lon
                lat_1 = lat(index); clear lat
                Cva_1 = Cva(index); clear Cua
                Cvp_1 = Cvp(index); clear Cup
                Cvp_2 = Cvp_1;
                chan=find(Cvp_2>180 & Cvp_2<=360);
                Cvp_2(chan)=Cvp_2(chan)-360;
                zz_Cva = griddata(lat_1,lon_1,Cva_1,lat_rho,lon_rho);
                zz_Cvp = griddata(lat_1,lon_1,Cvp_1,lat_rho,lon_rho);
                zz_Cvp2 = griddata(lat_1,lon_1,Cvp_2,lat_rho,lon_rho);
                chan2=find(zz_Cvp2<0);
                zz_Cvp2(chan2) = zz_Cvp2(chan2)+360;
                convert = find((zz_Cvp>10&zz_Cvp<350)&((zz_Cvp2>0&zz_Cvp2<10)|(zz_Cvp2>350&zz_Cvp2<360)));
                zz_Cvp(convert) = zz_Cvp2(convert);
                tide_Cva(ij,:,:) = zz_Cva(:,:);
                tide_Cvp(ij,:,:) = zz_Cvp(:,:);
                clear zz_Cva zz_Cvp zz_Cvp2 chan chan2 convert
            end
            disp([char(comp(ij,:)),' making complete !'])
        end
        if uv==1
            [major,eccentricity,inclination,phase]=ap2ep(tide_Cua,tide_Cup,tide_Cva,tide_Cvp);
            tide_Cmin=inpaint_nans3(major.*eccentricity,1);
            tide_Cmax=inpaint_nans3(major,1);
            tide_Cangle=inpaint_nans3(inclination,1);
            tide_Cphase=inpaint_nans3(phase,1);
        end
    case 'large domain'
        if tidenum>16;tidenum=16;end
        comp = {'M2';'S2';'K1';'O1';'N2';'K2';'P1';'Q1';'J1';'OO1';'M1';'2N2';'MU2';'T2';'NU2';'L2'};
        components = 'M2  S2  K1  O1  N2  K2  P1  Q1  J1  OO1  M1  2N2  MU2  T2  NU2  L2';
        period=[12.4206012; 12.0000000; 23.9344697; 25.8193417; 12.6583482; 11.9672348; 24.0658902; 26.8683567;...
        23.0984768; 22.3060742; 24.8412024; 12.9053745; 12.8717576; 12.0164492; 12.6260044; 12.1916202];
        tide_dir = uigetdir('./', 'Tide Directory');
        for ij=1:tidenum
            if ele==1
                tide_data = load([tide_dir,'/ele/',char(comp(ij,:)),'.mat']);
                lon_1 = inpaint_nans(tide_data.lon,1);
                lat_1 = inpaint_nans(tide_data.lat,1);
                eval(['amp_1 = inpaint_nans(tide_data.amp',char(comp(ij,:)),',1);'])
                eval(['pha_1 = inpaint_nans(tide_data.pha',char(comp(ij,:)),',1);'])
                pha_2 = pha_1;
                chan=find(pha_2>180 & pha_2<=360);
                pha_2(chan)=pha_2(chan)-360;
                zz_amp = griddata(lat_1,lon_1,amp_1,lat_rho,lon_rho);
                zz_pha = griddata(lat_1,lon_1,pha_1,lat_rho,lon_rho);
                zz_pha2 = griddata(lat_1,lon_1,pha_2,lat_rho,lon_rho);
                chan2=find(zz_pha2<0);
                zz_pha2(chan2) = zz_pha2(chan2)+360;
                convert = find((zz_pha>10&zz_pha<350)&((zz_pha2>0&zz_pha2<10)|(zz_pha2>350&zz_pha2<360)));
                zz_pha(convert) = zz_pha2(convert);
                tide_amp(ij,:,:) = inpaint_nans(zz_amp(:,:),1);
                tide_pha(ij,:,:) = inpaint_nans(zz_pha(:,:),1);
                clear zz_amp zz_pha zz_pha2 chan chan2
            end
            if uv==1
                curr_data = load([tide_dir,'/uv/',char(comp(ij,:)),'.mat']);
                lon_1 = curr_data.lon;
                lat_1 = curr_data.lat;
                eval(['Cmax_1 = inpaint_nans(curr_data.maj_',char(comp(ij,:)),',1);'])
                eval(['Cmin_1 = inpaint_nans(curr_data.min_',char(comp(ij,:)),',1);'])
                eval(['Cangle_1 = inpaint_nans(curr_data.inc_',char(comp(ij,:)),',1);'])
                eval(['Cphase_1 = inpaint_nans(curr_data.pha_',char(comp(ij,:)),',1);'])
                zz_Cmax = griddata(lat_1,lon_1,Cmax_1,lat_rho,lon_rho);
                zz_Cmin = griddata(lat_1,lon_1,Cmin_1,lat_rho,lon_rho);
                zz_Cangle = griddata(lat_1,lon_1,Cangle_1,lat_rho,lon_rho);
                zz_Cphase = griddata(lat_1,lon_1,Cphase_1,lat_rho,lon_rho);
                tide_Cmax(ij,:,:) = inpaint_nans(zz_Cmax(:,:),1);
                tide_Cmin(ij,:,:) = inpaint_nans(zz_Cmin(:,:),1);
                tide_Cangle(ij,:,:) = inpaint_nans(zz_Cangle(:,:),1);
                tide_Cphase(ij,:,:) = inpaint_nans(zz_Cphase(:,:),1);
                clear zz_Cmax zz_Cmin zz_Cangle zz_Cphase
            end
            disp([char(comp(ij,:)),' making complete !'])
        end
    otherwise
%         cc;disp('Unknown method.'); break
end

% system(['copy ',gridname,' ',outputname])  %% make a tide file from copied grd file

% outfile=netcdf(outputname, 'write');
% outfile.source = source;
% outfile.date = date;
% outfile.date = ncchar(date);
% outfile.start_tide_mjd=0;
% outfile.components = ncchar(components);
% outfile.components = components;

nccreate(outputname, 'write')
ncwriteatt(outputname,'/','source', source);
ncwriteatt(outputname,'/','date',char(date));
ncwriteatt(outputname,'/','start_tide_mjd',0);
ncwriteatt(outputname,'/','components',components);
ncwriteatt(outputname,'/','author', 'Created by Y.Y. Kim');



% outfile('xi_u') = M-1;
% outfile('eta_u') = N;
% outfile('xi_v') = M;
% outfile('eta_v') = N-1;
% outfile('xi_rho') = M;
% outfile('eta_rho') = N;
% outfile('xi_psi') = M-1;
% outfile('eta_psi') = N-1;
% outfile('tide_period')=tidenum;
    
nccreate(outputname, 'tide_period', 'Dimensions', {'tide_period' tidenum}, 'Datatype', 'double');
ncwriteatt(outputname, 'tide_period', 'long_name', 'Tide angular period');
ncwriteatt(outputname, 'tide_period', 'units', 'Hours');
ncwrite(outputname, 'tide_period', period(1:tidenum));

if ele==1
    nccreate(outputname, 'tide_Ephase', 'Dimensions', {'xi_rho', M, 'eta_rho', N, 'tide_period', tidenum}, 'Datatype', 'double');
    ncwriteatt(outputname, 'tide_Ephase', 'long_name', 'Tidal elevation phase angle');
    ncwriteatt(outputname, 'tide_Ephase', 'units','Degrees');
    ncwrite(outputname, 'tide_Ephase', permute(tide_pha, [3 2 1]));

    nccreate(outputname, 'tide_Eamp', 'Dimensions', {'xi_rho', M, 'eta_rho', N, 'tide_period', tidenum}, 'Datatype', 'double');
    ncwriteatt(outputname, 'tide_Eamp', 'long_name', 'Tidal elevation amplitude');
    ncwriteatt(outputname, 'tide_Eamp', 'units', 'Meter');
    ncwrite(outputname, 'tide_Eamp', permute(tide_amp, [3 2 1]));
end
if uv==1
    nccreate(outputname, 'tide_Cmin', 'Dimensions', {'xi_rho', M, 'eta_rho', N, 'tide_period', tidenum}, 'Datatype', 'double');
    ncwriteatt(outputname, 'tide_Cmin', 'long_name', 'Tidal current ellipse semi-minor axis');
    ncwriteatt(outputname, 'tide_Cmin', 'units', 'Meter second-1');
    ncwrite(outputname, 'tide_Cmin', permute(tide_Cmin, [3 2 1]));

    nccreate(outputname, 'tide_Cmax', 'Dimensions', {'xi_rho', M, 'eta_rho', N, 'tide_period', tidenum}, 'Datatype', 'double');
    ncwriteatt(outputname, 'tide_Cmax', 'long_name', 'Tidal current ellipse semi-minor axis');
    ncwriteatt(outputname, 'tide_Cmax', 'units', 'Meter second-1');
    ncwrite(outputname, 'tide_Cmax', permute(tide_Cmax, [3 2 1]));

    nccreate(outputname, 'tide_Cangle', 'Dimensions', {'xi_rho', M, 'eta_rho', N, 'tide_period', tidenum}, 'Datatype', 'double');
    ncwriteatt(outputname, 'tide_Cangle', 'long_name', 'Tidal current inclination angle');
    ncwriteatt(outputname, 'tide_Cangle', 'units', 'Degrees between semi-major axis and East');
    ncwrite(outputname, 'tide_Cangle', permute(tide_Cangle, [3 2 1]));

    nccreate(outputname, 'tide_Cphase', 'Dimensions', {'xi_rho', M, 'eta_rho', N, 'tide_period', tidenum}, 'Datatype', 'double');
    ncwriteatt(outputname, 'tide_Cphase', 'long_name', 'Tidal current phase angle');
    ncwriteatt(outputname, 'tide_Cphase', 'units', 'Degrees');
    ncwrite(outputname, 'tide_Cphase', permute(tide_Cphase, [3 2 1]));
end


if adj_depth=='y'
    copyfile([inpath,infile],[inpath,infile(1:end-3),'_MSL.nc']);
    nc=netcdf([inpath,infile(1:end-3),'_MSL.nc'],'write');
    h=h+squeeze(nansum(tide_amp(1:4,:,:)));
    nc{'h'}(:)=h;
    disp('tide adjusting complete !')
    close(nc)
end