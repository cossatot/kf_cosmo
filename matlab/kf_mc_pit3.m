% Age profile calculator script for PV06 profile

clear

rng('default');
rng(69);

addpath('C:\research\cosmo\CRONUS-calc-June-23-2015\common\production')
addpath('C:\research\cosmo\CRONUS-calc-June-23-2015\sa\production\')
addpath('C:\research\cosmo\CRONUS-calc-June-23-2015\common\profilecalc\')

addpath('~/research/cosmo/CRONUS-calc-June-23-2015/common/production/')
addpath('~/research/cosmo/CRONUS-calc-June-23-2015/sa/production/')
addpath('~/research/cosmo/CRONUS-calc-June-23-2015/common/profilecalc/')


%addpath('~\research\cosmo\DepthCodeSAMarrero\profilecalc\')
%addpath('~\research\cosmo\DepthCodeSAMarrero\production\')

scaling_model = 'sa';

% Load data
load kf_pit3_chem.mat

% Calculate depth-to-middle of the samples
depths = nominal36(:,37) + 0.5*nominal36(:,6).*nominal36(:,5);
% ndepths=length(depths);

max_erosion_cm = 200;
max_erosion_g_cm2 = mean(nominal36(:,5)) * max_erosion_cm;

% Define the paramter space
er_min = 0;
er_max = 1;
er_minmax = [er_min er_max];

age_min = 0;
age_max = 150;
age_minmax = [age_min age_max];

inher_min = 0;
inher_max = 1e2;
inher_minmax = [inher_min inher_max];

n_iters = 1000;

[posterior_er, posterior_age, posterior_inher, MAP, rel_likes] ...
    = depth_profile_mc_36(nominal36, uncerts36(:,1), depths, ...
    er_minmax, age_minmax, inher_minmax, n_iters, scaling_model, ...
    max_erosion_g_cm2);

fprintf('MAP_erate = %.2f g/cm^2/kyr \n',MAP(1))
fprintf('MAP_age   = %.2f k years \n',MAP(2))
fprintf('MAP_inher = %.2f years of exposure \n \n',MAP(3))

% Create plots
%makeplots36(erates,ages,inhers,posterior_er,posterior_age,...
%    posterior_inher,nominal36,uncerts36,depths,MAP,jposterior);
save kf_p2_mc_results

figure;
scatter(posterior_age, posterior_er, [], rel_likes, 'filled');
colorbar;

figure;
plotprof36(nominal36, uncerts36(:,1), depths, MAP(1), MAP(2), MAP(3),...
    scaling_model);
