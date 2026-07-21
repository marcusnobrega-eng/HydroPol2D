function idx = pf_systematic_resample(weights, u0)
%PF_SYSTEMATIC_RESAMPLE Deterministic-size systematic resampling.
%
% u0 is optional and must be in [0, 1/N). It is exposed for self-tests.

weights = double(weights(:));
N = numel(weights);
weights(~isfinite(weights) | weights < 0) = 0;
if sum(weights) <= 0
    weights(:) = 1 / N;
else
    weights = weights ./ sum(weights);
end

if nargin < 2 || isempty(u0)
    u0 = rand() / N;
end
if u0 < 0 || u0 >= 1 / N
    error('pf_systematic_resample:u0', 'u0 must be in [0, 1/N).');
end

positions = u0 + (0:(N-1))' ./ N;
edges = cumsum(weights);
edges(end) = 1;

idx = zeros(N, 1);
j = 1;
for i = 1:N
    while positions(i) > edges(j)
        j = j + 1;
    end
    idx(i) = j;
end
end
