function data=ext_data_OGCM_CMIP6(nc,X,Y,vname,tndx,lon,lat,k,missvalue,Roa,interp_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extrapole one horizontal ECCO (or Data) slice on a ROMS grid
%
%
%  Further Information:
%  http://www.brest.ird.fr/Roms_tools/
%
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven
%  e-mail:Pierrick.Penven@ird.fr
%
%  Contributions of P. Marchesiello (IRD) and J. Lefevre (IRD)
%
%  Updated    6-Sep-2006 by Pierrick Penven
%  Updated    18-Mar-2021 by Yong-Yub Kim (scale factor, add_offset)
%  Updated    16-Jun-2021 by Yong-Yub Kim (write NaN, if all data were default)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% extrapolation parameters
%
default = 0;
if strcmp(vname,'SAVE') | strcmp(vname,'salt') | strcmp(vname,'SALT')
    default = 34.6;
end
%
% Get the ROMS grid extension + a little margin (~ 2 data grid points)
%
dx = max(abs(gradient(X)));
dy = max(abs(gradient(Y)));
dl = 2*max([dx dy]);
%
lonmin = min(min(lon)) - dl;
lonmax = max(max(lon)) + dl;
latmin = min(min(lat)) - dl;
latmax = max(max(lat)) + dl;
%
% Extract a data subgrid
%
j = find(Y >= latmin & Y <= latmax);
i1 = find(X-360 >= lonmin & X-360 <= lonmax);
i2 = find(X >= lonmin & X <= lonmax);
i3 = find(X+360 >= lonmin & X+360 <= lonmax);
if ~isempty(i2)
    x = X(i2);
else
    x = [];
end
if ~isempty(i1)
    x = cat(2,X(i1)-360,x);
end
if ~isempty(i3)
    x = cat(2,x,X(i3)+360);
end
y = Y(j);
%
%  Get dimensions
%
ndims = length(dim(nc{vname}));
%
% Get data (Horizontal 2D matrix)
%

val_add_offset=nc{vname}.add_offset(:);
if isempty(val_add_offset)
    val_add_offset=0;
end
val_scale_factor=nc{vname}.scale_factor(:);
if isempty(val_scale_factor)
    val_scale_factor=1;
end

if ~isempty(i2)
    if ndims == 2
        data = squeeze(nc{vname}(j,i2));
    elseif ndims == 3
        data=squeeze(nc{vname}(tndx,j,i2));
    elseif ndims==4
        data = squeeze(nc{vname}(tndx,k,j,i2));
    else
        error(['Bad dimension number ',num2str(ndims)])
    end
else
    data = [];
end
if ~isempty(i1)
    if ndims == 2
        data = cat(2,squeeze(nc{vname}(j,i1)),data);
    elseif ndims == 3
        data=cat(2,squeeze(nc{vname}(tndx,j,i1)),data);
    elseif ndims == 4
        data = cat(2,squeeze(nc{vname}(tndx,k,j,i1)),data);
    else
        error(['Bad dimension number ',num2str(ndims)])
    end
end
if ~isempty(i3)
    if ndims == 2
        data = cat(2,data,squeeze(nc{vname}(j,i3)));
    elseif ndims == 3
        data=cat(2,data,squeeze(nc{vname}(tndx,j,i3)));
    elseif ndims == 4
        data = cat(2,data,squeeze(nc{vname}(tndx,k,j,i3)));
    else
        error(['Bad dimension number ',num2str(ndims)])
    end
end



%
% Perform the extrapolation
%
[data,interp_flag] = get_missing_val(x,y,data,missvalue,Roa,default);
if(sum(sum(data==0))==length(data(:)))
    data=NaN(size(data));
    disp(['missing data layer get to be NaN now. afterwards ',...
         ' these values will be filled by value of valid deepest layer']);
end

data = data .* val_scale_factor + val_add_offset;


%
% Interpolation on the ROMS grid
%
if(prod(prod(prod(isnan(data))))~=1)
    if interp_flag==0
        data = interp2(x,y,data,lon,lat,'nearest');
    else
        data = interp2(x,y,data,lon,lat,interp_method);
        %   data=griddata(x,y,data,lon,lat);
    end
end
%
return