function [mcamaps_l3d, pcs_l, mcamaps_r3d, pcs_r, lambdas, scf] = mca(data_l, data_r, varargin)
% [mcamaps_l3d, pcs_l, mcamaps_r3d, pcs_r, lambdas] = mca(data_l, data_r, varargin)
% Version 2.0
% MCA or SVD Toolbox for marine and atmospheric science.
% Function names MCA but not SVD here because the MATLAB build-in "svd" function
%%   Syntax
%       [mcamaps_l3d, pcs_l, mcamaps_r3d, pcs_r, lambdas, scf] = mca(data_l, data_r)
% 
%       [mcamaps_l3d, pcs_l, mcamaps_r3d, pcs_r, lambdas, scf] = mca(data_l, data_r, n_mca) 
%           Only calculate the first n_mca mode. 
%
%       This function will auto-omit nan values of data_l and data_r
%% Input:
%       data_l(lon1,lat1,time): 3d left data field
%       date_r(lon2,lat2,time): 3d right data field
%       The timesteps of two datasets should be the same.
%% Output
%       mcamaps_l3d: spatial modes of data_l
%       pcs_k: Priciple components time series of data_l
%       lambdas: eigenvalues of Covariance_Lr
%       scf: squared covariance fraction of every modes. SCF = lambdas ./
%       sum(lambdas);
%       RemarKs: To save the memory, this function estimate the Lambdas and
%           SCF only basing on the first n_mca modes.
% See also SVD
%% Author:
%	Zelun Wu,
%   Ph.D. student of Physical Oceanography,
%	Xiamen University & University of Delaware
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	6th, 09, 2020
%   2023 Sep 14th, minor revisions by Y.Y. Kim (covariance against time axis, getting 'squared' fraction)

  
%% Error checks 
narginchk(1,inf) 
%% Input parsing
[N_lon1, N_lat1, N_time1] = size(data_l);
[N_lon2, N_lat2, N_time2] = size(data_r);
assert(N_time1 == N_time2, 'Timesteps of two dataset should be the same');
if nargin > 2
    n_mca = varargin{1};
else
    error('Needs more than two dataset');
end
% core

[data_left_2d, in_notnan_l] = reshape3dto2d(data_l);
[data_right_2d, in_notnan_r] = reshape3dto2d(data_r);

[mcamaps_l, pcs_l, mcamaps_r, pcs_r, lambdas] = svdcore(data_left_2d, data_right_2d, n_mca);
% scf = lambdas./sum(lambdas);
scf = lambdas.^2./sum(lambdas.^2);

mcamaps_l3d = reshape2dto3d(mcamaps_l, [N_lon1, N_lat1, n_mca], in_notnan_l);
mcamaps_r3d = reshape2dto3d(mcamaps_r, [N_lon2, N_lat2, n_mca], in_notnan_r);

end
