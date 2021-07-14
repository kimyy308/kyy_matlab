function data=ext_data_OGCM_MyOcean(nc,X,Y,vname,tndx,lon,lat,NZ,missvalue,Roa,interp_method, direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extrapole sliced MyOcean data on a ROMS grid
%
%  Updated   17-Jul-2018 by Yong-Yub Kim
%
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
% dx = max(abs(gradient(X)));
% dy = max(abs(gradient(Y)));
% dl = 2*max([dx dy]);
%
% lonmin = min(min(lon)) - dl;
% lonmax = max(max(lon)) + dl;
% latmin = min(min(lat)) - dl;
% latmax = max(max(lat)) + dl;

% %
% % Extract a data subgrid
% %
% j = find(Y >= latmin & Y <= latmax);
% i1 = find(X-360 >= lonmin & X-360 <= lonmax);
% i2 = find(X >= lonmin & X <= lonmax);
% i3 = find(X+360 >= lonmin & X+360 <= lonmax);
% if ~isempty(i2)
%     x = X(i2);
% else
%     x = [];
% end
% if ~isempty(i1)
%     x = cat(2,X(i1)-360,x);
% end
% if ~isempty(i3)
%     x = cat(2,x,X(i3)+360);
% end
% y = Y(j);

%
%  Get dimensions
%
ndims = length(dim(nc{vname}));


%
% Get data (Horizontal 2D matrix)
%

i2=length(X);
j=length(Y);

if ~isempty(i2)
    add_offset=nc{vname}.add_offset(:);
    scale_factor=nc{vname}.scale_factor(:);
    if ndims == 2
        data = squeeze(nc{vname}(j,i2)) .* scale_factor + add_offset;
    elseif ndims == 3
        data=squeeze(nc{vname}(tndx,j,i2)) .* scale_factor + add_offset;
    elseif ndims==4
        data = squeeze(nc{vname}(tndx,NZ,j,i2)) .* scale_factor + add_offset;
    else
        error(['Bad dimension number ',num2str(ndims)])
    end
else
    data = [];
end

%
% Perform the extrapolation
%
x = X(i2);
y = Y(j);
missvalue = missvalue .*scale_factor + add_offset;
[data,interp_flag] = get_missing_val(x,y,data,missvalue,Roa,default);



%
% Interpolation on the ROMS grid
%

switch(direction)
    case(1)  %% south
        for k=1:NZ
            if interp_flag==0
                data(k,1,1:length(x)) = griddata(x,y,data,lon,lat,'nearest');
            else
                data(k,1,1:length(x)) = griddata(x,y,data,lon,lat,interp_method);
            end
        end
    case(2)  %% east
        for k=1:NZ
            if interp_flag==0
                data(k,1:length(y),1) = griddata(x,y,data,lon,lat,'nearest');
            else
                data(k,1:length(y),1) = griddata(x,y,data,lon,lat,interp_method);
            end
        end
    case(3)  %% north
       for k=1:NZ
            if interp_flag==0
                data(k,1,1:length(x)) = griddata(x,y,data,lon,lat,'nearest');
            else
                data(k,1,1:length(x)) = griddata(x,y,data,lon,lat,interp_method);
            end
        end
    case(4)  %% west
       for k=1:NZ
            if interp_flag==0
                data(k,1:length(y),1) = griddata(x,y,data,lon,lat,'nearest');
            else
                data(k,1:length(y),1) = griddata(x,y,data,lon,lat,interp_method);
            end
        end
    otherwise
            disp('??????????????')
end



return