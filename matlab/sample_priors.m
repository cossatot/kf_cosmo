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


function unif_rand_samples = unif_rand(umax, umin, unum) 
    unif_rand_samples = umin + (umax - umin) .* rand(unum, 1);
end


function mono_increase(x)
    if all(diff(x)>0)
        disp('all mono');
    else
        disp('not all mono');
    end
end
