function [figname] = Func_0002_get_jpgname(figrawdir, varname, typename, regionname, scenname, testname, timename, extensionname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function figname=Func_0002_get_jpgname(testname);
%
% get the historical testname corresponding to RCP test name
%
%  input:
%  figrawdir            main figure directory (string)
%  varname              variable name to get figure (string)
%  typename             type of the figure (string)
%  regionname           region of the figure (string)
%  scenname             scenario name (string)
%  testname             name of the test (string)
%  timename             period of the figure (string)
%  extensionname        extension of the file (string, ex) .tif, .jpg)
%
%  output:
%  figname              figure file name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    1-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirname = [figrawdir, '\', varname, '\', typename, '\' regionname, '\'];
figname = [dirname, varname, '_', typename, '_', regionname, '_', scenname, '_', testname, '_', timename, extensionname]; 
if (exist(figname , 'file') ~= 2)
    mkdir(strcat(dirname));
end
end

