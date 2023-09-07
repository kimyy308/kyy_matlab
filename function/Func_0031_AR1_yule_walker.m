function [prog, cc] = Func_0031_AR1_yule_walker(init, tsteps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [psd, psd_sig, freq] = Func_0027_get_PSD_siglev(data, sig_level, DOF)
%
%  Estimation of AR1 parameters
%
%  input:
%  init         initial time-series
%  tstep        prognostic steps (lag)
%
%  output:
%  prog         time-series calculated by AR1 model estimated by yule-walker
%  cc           ACC skill of prognostic value 
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      22-Aug-2023 by Yong-Yub Kim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize for estimation (remove NaN)
init_finite = init(isfinite(init));

%% model parameter estimation
order=1;
[coef, noise]=aryule(init_finite, order);

%% integration
prog=init(1);
for i=1:tsteps
    prog = -prog*(coef(2)) + (noise);
end
% prog = -init*(coef(2)^(tsteps)) + (noise*tsteps);

% %% get ACC skill of prognostic value
% i1=init(tsteps:end);
% i2=prog(1:tsteps+1);
% 
% cc=corrcoef(i1(isfinite(i1)), i2(isfinite(i1)));
% cc=cc(1,2);
end


