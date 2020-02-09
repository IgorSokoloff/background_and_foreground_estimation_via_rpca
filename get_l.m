
function l = get_l(sur_kind, sigma_v, gamma)
    assert (isvector(sigma_v));
    if strcmp (sur_kind, "laplace")
        l = exp(-sigma_v./gamma) ./gamma;
    elseif strcmp (sur_kind, "geman")
        l = (1 + gamma)*exp(sigma_v) ./(gamma + exp(sigma_v)).^2;
    else
        error("wrong surrogate function");
    end 
end