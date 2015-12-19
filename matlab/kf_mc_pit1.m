% Age profile calculator script for PV06 profile

clear

rng('default');
rng(69);

addpath('C:\research\cosmo\dal_geo_centre\DepthCodeSAMarrero\profilecalc\')
addpath('C:\research\cosmo\dal_geo_centre\DepthCodeSAMarrero\production\')

addpath('~/research/cosmo/CRONUS-calc-June-23-2015/common/production/')
addpath('~/research/cosmo/CRONUS-calc-June-23-2015/sa/production/')
addpath('~/research/cosmo/CRONUS-calc-June-23-2015/common/profilecalc/')


%addpath('~\research\cosmo\DepthCodeSAMarrero\profilecalc\')
%addpath('~\research\cosmo\DepthCodeSAMarrero\production\')

scaling_model = 'sa';

% Load data
load kf_pit1_chem.mat

% correct so that depth is in g/cm^2  (ONLY FOR KF PIT 1)
nominal36(:,37) = nominal36(:,37) .* nominal36(:,5);

% Calculate depth-to-middle of the samples
depths = nominal36(:,37) + 0.5*nominal36(:,6).*nominal36(:,5);
% ndepths=length(depths);

max_erosion_cm = 200;
max_erosion_g_cm2 = mean(nominal36(:,5)) * max_erosion_cm;

% Define the paramter space
er_min = 0;
er_max = 2;
er_minmax = [er_min er_max];

age_min = 100;
age_max = 300;
age_minmax = [age_min age_max];

inher_min = 0;
inher_max = 58e3;
inher_minmax = [inher_min inher_max];

n_iters = 100000;

[age_priors, er_priors, inher_priors] = sample_priors(age_minmax, ...
                      er_minmax, inher_minmax, n_iters, max_erosion_g_cm2);

[posterior_er, posterior_age, posterior_inher, MAP, rel_likes, ...
    likelihoods] = depth_profile_mc_36(nominal36, uncerts36(:,1), ...
                          depths, age_priors, er_priors, inher_priors, ...
                          scaling_model);

fprintf('MAP_erate = %.2f g/cm^2/kyr \n',MAP(1))
fprintf('MAP_age   = %.2f k years \n',MAP(2))
fprintf('MAP_inher = %.2f years \n \n',MAP(3))

% Create plots
%makeplots36(erates,ages,inhers,posterior_er,posterior_age,...
%    posterior_inher,nominal36,uncerts36,depths,MAP,jposterior);
save kf_p1_mc_results

figure;
scatter(posterior_age, posterior_er, [], rel_likes, 'filled');
colorbar;

figure;
plotprof36(nominal36, uncerts36(:,1), depths, MAP(1), MAP(2), MAP(3),...
    scaling_model);

figure;
subplot(131);
[age_bw, age_pdf, age_xs, age_cdf] = kde(posterior_age, 1000, age_min,...
                                         age_max);
plot(age_xs, age_pdf);
xlabel('Age (ka)');

subplot(132);
[er_bw, er_pdf, er_xs, er_cdf] = kde(posterior_er, 1000, er_min,...
                                         er_max);
plot(er_xs, er_pdf);
xlabel('Erosion rate (some units)');

subplot(133);
[in_bw, in_pdf, in_xs, in_cdf] = kde(posterior_inher, 1000, inher_min,...
                                     inher_max);
plot(in_xs, in_pdf);
xlabel('Inheritance (years)');


