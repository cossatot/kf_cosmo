% Age profile calculator script for PV06 profile

clear
tic

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

% Calculate depth-to-middle of the samples
depths = nominal36(:,37) + 0.5*nominal36(:,6).*nominal36(:,5);
% ndepths=length(depths);


% Define the paramter space
%erates=linspace(-1,5,10)';
%ages=linspace(10,500,5)';
%inhers=linspace(0,3.5e3,10)';

inhers = linspace(0,0.01,2)';

unif_rand = @(umax, umin, unum) (umin + (umax - umin) .* rand(unum, 1));

er_min = 0;
er_max = 1;
er_num = 2;

%erates = unif_rand(er_min, er_max, er_num);
erates = linspace(er_min, er_max, er_num)';

age_min = 150;
age_max = 300;
age_num = 100;

%ages = unif_rand(age_min, age_max, age_num);
ages = linspace(age_min, age_max, age_num)';

% Define the prior distribution
pr=prior(erates,ages,inhers,'uniform','uniform','uniform');


% Run ageprofilecalc1026
[posterior_er,posterior_age,posterior_inher,MAP,mu_bayes,chi2grid,lhsurf,...
    jposterior]=ageprofilecalc36(nominal36,uncerts36(:,1),depths,...
    erates,ages,inhers,pr, scaling_model);


% Display the MAP solution and Bayesian credible intervals
erCI=bayesianCI(posterior_er,erates,0.05);
ageCI=bayesianCI(posterior_age,ages,0.05);
inherCI=bayesianCI(posterior_inher,inhers,0.05);

toc

fprintf('MAP_erate = %.2f g/cm^2/kyr \n',MAP(1))
fprintf('MAP_age   = %.2f k years \n',MAP(2))
fprintf('MAP_inher = %.2f years of exposure \n \n',MAP(3))
%fprintf('mean erate = %.2f \n',bmm(1))
%fprintf('mean_age = %.2f \n',bmm(2))
%fprintf('mean_inher = %.2f \n \n',bmm(3))
fprintf('68%% credible intervals \n \n')
fprintf('P(%.2f < erate < %.2f)=0.68 \n',erCI)
fprintf('P(%.2f < age < %.2f)=0.68 \n',ageCI)
fprintf('P(%.2f < inher < %.2f)=0.68 \n \n',inherCI)


% Create plots
makeplots36(erates,ages,inhers,posterior_er,posterior_age,...
    posterior_inher,nominal36,uncerts36,depths,MAP,jposterior);
%save kf_p1_results

