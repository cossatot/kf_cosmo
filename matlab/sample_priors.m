function [age_prior, er_prior, inher_prior] = sample_priors(ages, erates, ...
                                                  inhers, n_iters, max_erosion)

    age_prior = sample_pdf(ages(1), ages(2), n_iters, get_prior_type(ages));
    er_prior = sample_pdf(erates(1), erates(2), n_iters, ...
                          get_prior_type(erates));
    inher_prior = sample_pdf(inhers(1), inhers(2), n_iters, ...
                             get_prior_type(inhers));
    if nargin == 4
        max_erosion = +inf;
    elseif nargin == 5
        [age_prior, er_prior, inher_prior] = remove_too_much_erosion(age_prior,...
            er_prior, inher_prior, max_erosion);
    end
end


function priortype = get_prior_type(prior)
    if length(prior) == 2
        priortype = 'unif_rnd';
    elseif length(prior) > 2
        priortype = 'arb_rnd';
    end
end


function pdf_samples = sample_pdf(vals, probs, n_samples, type)

    switch type

        case 'unif_rnd'
            pdf_samples = unif_rand(vals, probs, n_samples);

        case 'arb_rnd'
            pdf_samples = inverse_transform_sample(vals, probs, n_samples);
    end
end


function pdf_samples = inverse_transform_sample(vals, probs, n_samples)

    [pdf_range, pdf_probs] = make_pdf(vals, probs);
    [cdf_range, cdf_probs] = make_cdf(pdf_range, pdf_probs);

    samps = rand(n_samples, 1);
    pdf_samples = interp1(cdf_probs, cdf_range, samps);
end


function [pdf_range, pdf_probs] = make_pdf(vals, probs)

    n_interp = 1000;

    val_min = min(vals);
    val_max = max(vals);

    pdf_range = linspace(val_min, val_max, n_interp);

    pmf_samples = interp1(vals, probs, pdf_range);
    pdf_probs = pmf_samples / sum(pmf_samples);
end


function [cdf_range, cdf_vals] = make_cdf(pdf_range, pdf_vals)
    cdf_range = pdf_range;
    cdf_vals = cumtrapz(pdf_vals);
end


function unif_rand_samples = unif_rand(umin, umax, unum) 
    unif_rand_samples = umin + (umax - umin) .* rand(unum, 1);
end


function mono_increase(x)
    if all(diff(x)>0)
        disp('all mono');
    else
        disp('not all mono');
    end
end


function [ages, erates, inhers] = remove_too_much_erosion(ages, erates, ...
                                                          inhers, max_erosion)
    not_too_much_erosion = (ages .* erates < max_erosion);
    ages = ages(not_too_much_erosion);
    erates = erates(not_too_much_erosion);
    inhers = inhers(not_too_much_erosion);
end
                                                       
