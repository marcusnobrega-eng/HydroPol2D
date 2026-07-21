function weights = pf_normalize_logweights(log_weights)
%PF_NORMALIZE_LOGWEIGHTS Stable log-sum-exp normalization.

log_weights = double(log_weights(:));
finite_idx = isfinite(log_weights);
if ~any(finite_idx)
    weights = ones(size(log_weights)) ./ numel(log_weights);
    return;
end

shift = max(log_weights(finite_idx));
w = zeros(size(log_weights));
w(finite_idx) = exp(log_weights(finite_idx) - shift);
total = sum(w);
if total <= 0 || ~isfinite(total)
    weights = ones(size(log_weights)) ./ numel(log_weights);
else
    weights = w ./ total;
end
end
