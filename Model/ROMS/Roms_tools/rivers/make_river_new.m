close all
clear all
clc

file=1992:2002;

grd_file='D:\Roms\06-12-25(1_4)\roms_grid_4degree.nc';
Fname='roms_river_4degree.nc';

grd=roms_get_grid(grd_file);
cycle=3960;
% cycle=3585;

create_river_mks(grd,Fname,cycle,grd_file);

ncin=netcdf('roms_1993_river.nc');
ncout=netcdf(Fname,'write');

theVarname = 'river';
ncout{theVarname}(:) = ncin{theVarname}(:);

theVarname = 'river_Xposition';
ncout{theVarname}(:) = ncin{theVarname}(:);

theVarname = 'river_Eposition';
ncout{theVarname}(:) = ncin{theVarname}(:);

theVarname = 'river_direction';
ncout{theVarname}(:) = ncin{theVarname}(:);

theVarname = 'river_flag';
ncout{theVarname}(:) = ncin{theVarname}(:);

theVarname = 'river_Vshape';
ncout{theVarname}(:,:) = ncin{theVarname}(:);

theVarname = 'river_time';
ncout{theVarname}(1:12) = ncin{theVarname}(:);

theVarname = 'river_transport';
ncout{theVarname}(1:12,:) = ncin{theVarname}(:,:);
    
theVarname = 'river_temp';
ncout{theVarname}(1:12,:,:) = ncin{theVarname}(:,:,:);
    
theVarname = 'river_salt';
ncout{theVarname}(1:12,:,:) = ncin{theVarname}(:,:,:);

close(ncin)
close(ncout)

f=1;
for i= 13:12:132
    Iname=['roms_',num2str(file(f+1)),'_river.nc'];
    disp(['in file is ',Iname]);
    ncin=netcdf(Iname);
    ncout=netcdf(Fname,'write');
    
    theVarname = 'river_time';
    ncout{theVarname}(i:i+11) = ncin{theVarname}(:);
        
    theVarname = 'river_transport';
    ncout{theVarname}(i:i+11,:) = ncin{theVarname}(:,:);
    
    theVarname = 'river_temp';
    ncout{theVarname}(i:i+11,:,:) = ncin{theVarname}(:,:,:);
    
    theVarname = 'river_salt';
    ncout{theVarname}(i:i+11,:,:) = ncin{theVarname}(:,:,:);
    close(ncin)
    close(ncout)
    f=f+1;
end


