function [posterior_er, posterior_age, posterior_inher, MAP, rel_likes,...
    likelihoods] = depth_profile_mc_36(comp, sigma, depths, age_priors, er_priors, ...
                          inher_priors, scaling_model)
%
%   Monte Carlo-based Age Profile Calculator for Cl36
%
%   [posterior_er,posterior_age,posterior_inher,MAP,mu_bayes,chi2grid]=...
%       ageprofilecalc36(comp,concsig,depths,erates,ages,inhers,...
%       scaling_model,maxerosion)
%
%   Inputs:
%
%       comp                profile composistion
%       concsig             vector of concentration uncertanties
%       depths              vector of depths
%       erates              vector of erosion rates
%       age                 vector of ages
%       inhers              vector of inhers
%       prior               prior distribution
%  		scaling_model		The scaling model being used
%       maxerosion          maximum erosion allowed
%
%   Outputs:
%   
%       posterior_er        posterior erosion rate distribution
%       posterior_age       posterior age distribution
%       posterior_inher     posterior inheritance distribution
%       MAP                 maximum a posterior solution   
%       mu_bayes                 bayesian mean solution
%



% Set global variables
global PP;
global STOREDSP;
global STOREDSF;
global STOREDCP;
global maxdepth;


%
% Figure out the number of samples.
%
[nsamples,thirtyeight]=size(comp);

STOREDSP=cell(nsamples);
STOREDSF=cell(nsamples);
STOREDCP=cell(nsamples);

pred_36cl_conc = {};

%conc=comp(:,1);
%ndepths=length(depths);

measuredconc36=comp(:,1);
sigma36=sigma(:,1);



n_iters = length(age_priors);

% Define maxdepth
maxdepth=max(er_priors)*max(age_priors)*max(comp(:,6))+3000;


% Run getpars on the first sample(choice abitrary) to get and store
% physical parameters(pp) and scale factors(sf)

disp('storing phys. params')
tic;
for k=1:nsamples
  [pp,sp,sf,cp]=getpars36(comp(k,:),maxdepth);
  PP=pp;
  STOREDSP{k}=sp;
  STOREDSF{k}=sf;
  STOREDCP{k}=cp;
  pred_36cl_conc{k} = zeros(n_iters,1);
end
toc;

disp('doing 36Cl Monte Carlo')

ProdtotalCl=prodz36(0,pp,sf,cp);

tic;
parfor sample = 1:nsamples
    %ProdtotalCl=prodz36(0,pp,sf,cp);
    for iter = 1:n_iters
        epsilon = er_priors(iter)/1000.;
        pred_36cl_conc{sample}(iter) = predN36depth(pp, sp, sf, cp, ... 
            age_priors(iter), depths(sample), epsilon, scaling_model) ...
            + inher_priors(iter) * ProdtotalCl;
    end
end
toc

disp('calculating misfit')

tic;

chi_sq = zeros(n_iters, 1);

for sample = 1:nsamples
    chi2 = (( pred_36cl_conc{sample} - measuredconc36(sample) )/ ...
        sigma36(sample) ).^2;

    chi_sq = chi_sq + chi2;
end
toc; 

likelihoods = prod(1 ./ sqrt(2 * pi .* sigma36)) * exp(-chi_sq/2);

[max_like, max_idx] = max(likelihoods);

MAP(1) = er_priors(max_idx);
MAP(2) = age_priors(max_idx);
MAP(3) = inher_priors(max_idx);

rel_likelihoods = likelihoods / max_like;

rand_filter = rand(n_iters, 1);

keep_inds = (rand_filter < rel_likelihoods);

rel_likes = rel_likelihoods(keep_inds);

posterior_age = age_priors(keep_inds);
posterior_er = er_priors(keep_inds);
posterior_inher = inher_priors(keep_inds);
