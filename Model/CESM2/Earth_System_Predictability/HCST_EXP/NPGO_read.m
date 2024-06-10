%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Volumes/kyy_raid/kimyy/Observation/NPGO/NPGO_index_231019.txt
%
% Auto-generated by MATLAB on 19-Oct-2023 13:53:46

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [36, 920];
opts.Delimiter = ["  ", "\\"];

% Specify column names and types
opts.VariableNames = ["Var1", "deftab720", "VarName3", "VarName4"];
opts.SelectedVariableNames = ["deftab720", "VarName3", "VarName4"];
opts.VariableTypes = ["string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["deftab720", "VarName3", "VarName4"], "ThousandsSeparator", ",");

% Import the data
NPGOindex231019 = readtable("/Volumes/kyy_raid/kimyy/Observation/NPGO/NPGO_index_231019.txt", opts);

%% Convert to output type
NPGOindex231019 = table2array(NPGOindex231019);

%% Clear temporary variables
clear opts

NPGO_yearly=reshape(NPGOindex231019(1:876,3), [876/12, 12]);
NPGO_yearly=mean(NPGO_yearly,2);
plot(NPGO_yearly)
NPGO_3ym=movmean(NPGO_yearly,3);
plot(NPGO_3ym)
%1950:2022
plot(1960:2020,NPGO_3ym(11:71))