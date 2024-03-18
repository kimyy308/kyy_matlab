function [U, pcs_X, V, pcs_Y, lambdas] = svdcore(X, Y, varargin)
% [U, pcs_X, V, pcs_Y, lambdas] = svdcore(X, Y, varargin)
% Version 1.0
% MCA or SVD analysis of two dataset X(lon1*lat1,time) and
% Y(lon2*lat2,time)
% The timesteps of two datasets should be the same.
%%   Syntax
%      U, pcs_X, V, pcs_Y, lambdas] = svdcore(X, Y, varargin)
% 
%      U, pcs_X, V, pcs_Y, lambdas] = svdcore(X, Y, varargin, n_mca) 
%           Only calculate the first n_mca mode. 
%% Author:
%	Zelun Wu,
%   Ph.D. student of Physical Oceanography,
%	Xiamen University & University of Delaware
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	15th May, 2020

%% Input parsing
[N_locX, N_time1] = size(X);
[N_locY, N_time2] = size(X);
assert(N_time1 == N_time2, 'Time dimension of the left field and right field should be the same');
N_time = N_time1;
% standardize
X = X - mean(X,2);
Y = Y - mean(Y,2);
%% core
Cov_XY = X*Y'/N_time;

% [U,S,V] = svd(Cov_XY, 'econ');
[U,S,V] = svds(Cov_XY, varargin{1});

n_svd = size(U,2);
if nargin>1
    n_svd = varargin{1};
end

U = U(:,1:n_svd);
V = V(:,1:n_svd);

pcs_X = U'*X;
pcs_Y = V'*Y;

lambdas = diag(S);
end