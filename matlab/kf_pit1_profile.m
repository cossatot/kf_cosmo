% Age profile calculator script for PV06 profile

clear
tic


addpath('C:\research\cosmo\dal_geo_centre\DepthCodeSAMarrero\profilecalc\')
addpath('C:\research\cosmo\dal_geo_centre\DepthCodeSAMarrero\production\')


% Load data
load kf_pit1_chem.mat

% Calculate depth-to-middle of the samples
depths = nominal36(:,37) + 0.5*nominal36(:,6).*nominal36(:,5);
% ndepths=length(depths);


% Define the paramter space
erates=linspace(-1,5,15)';
ages=linspace(10,200,20)';
inhers=linspace(0e3,3.5e2,10)';


% Define the prior distribution
pr=prior(erates,ages,inhers,'uniform','uniform','uniform');


% Run ageprofilecalc1026
[posterior_er,posterior_age,posterior_inher,MAP,mu_bayes,chi2grid,lhsurf,...
    jposterior]=ageprofilecalc36(nominal36,uncerts36(:,1),depths,...
    erates,ages,inhers,pr);


% Display the MAP solution and Bayesian credible intervals
erCI=bayesianCI(posterior_er,erates,0.05);
ageCI=bayesianCI(posterior_age,ages,0.05);
inherCI=bayesianCI(posterior_inher,inhers,0.05);

fprintf('MAP_erate = %.2f g/cm^2/kyr \n',MAP(1))
fprintf('MAP_age   = %.2f k years \n',MAP(2))
fprintf('MAP_inher = %.2f years of exposure \n \n',MAP(3))
% fprintf('mean erate = %.2f \n',bmm(1))
% fprintf('mean_age = %.2f \n',bmm(2))
% fprintf('mean_inher = %.2f \n \n',bmm(3))
% fprintf('68%% credible intervals \n \n')
% fprintf('P(%.2f < erate < %.2f)=0.68 \n',erCI)
% fprintf('P(%.2f < age < %.2f)=0.68 \n',ageCI)
% fprintf('P(%.2f < inher < %.2f)=0.68 \n \n',inherCI)


% Create plots
makeplots36(erates,ages,inhers,posterior_er,posterior_age,...
    posterior_inher,nominal36,uncerts36,depths,MAP,jposterior);
save T6profileresults

toc