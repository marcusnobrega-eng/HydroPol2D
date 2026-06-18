function filteredData = applyKernelFilter(data, kernelSize, filterType, varargin)
% APPLYKERNELFILTER
% -------------------------------------------------------------------------
% Robust NaN-aware kernel filtering for raster data.
%
% This version fixes two common problems:
%
%   1) Border dilution:
%      conv2(...,'same') assumes zeros outside the raster, which biases
%      values near rectangular edges.
%
%   2) Irregular catchment domains:
%      NaNs represent inactive cells outside the watershed. A normal kernel
%      can incorrectly mix NaNs/zeros into valid cells near the watershed
%      boundary.
%
% Default behavior:
%   - NaNs are ignored.
%   - Cells whose full kernel neighborhood is not valid are preserved from
%     the original raster.
%   - This means boundary pixels are NOT modified.
%   - Uniform rasters remain uniform.
%
% Usage:
%   filteredData = applyKernelFilter(data, kernelSize, filterType)
%
% Examples:
%   A_sm = applyKernelFilter(A, 3, 'gaussian');
%   A_sm = applyKernelFilter(A, 5, 'mean');
%   A_sm = applyKernelFilter(A, 3, 'median');
%
% Optional name-value arguments:
%
%   'EdgeMode'
%       'preserve'  : default. Preserve cells where the full kernel does not
%                     fit inside the valid domain.
%       'normalize' : filter using only available valid neighbors, with
%                     normalized weights. This smooths up to edges.
%
%   'Sigma'
%       Gaussian sigma. Default = kernelSize/3.
%
%   'ValidMask'
%       Logical mask of valid cells. Default = isfinite(data).
%
%   'PreserveUniform'
%       If true, uniform valid rasters are returned unchanged.
%       Default = true.
%
%   'UniformTolerance'
%       Tolerance used to detect uniform rasters. Default = 1e-12.
%
%   'MinValidFraction'
%       Minimum fraction of valid cells inside the kernel required to compute
%       a filtered value when EdgeMode='normalize'. Default = 0.5.
%
% Inputs:
%   data       - 2D numeric matrix
%   kernelSize - odd integer, for example 3 for a 3x3 kernel
%   filterType - 'mean', 'gaussian', or 'median'
%
% Output:
%   filteredData - filtered matrix, same size as data
% -------------------------------------------------------------------------

    %% ---------------- Parse inputs ----------------

    if nargin < 3
        error('applyKernelFilter requires data, kernelSize, and filterType.');
    end

    if ~isnumeric(data) || ndims(data) ~= 2
        error('data must be a 2D numeric matrix.');
    end

    if ~isscalar(kernelSize) || kernelSize < 1 || mod(kernelSize, 2) == 0
        error('kernelSize must be a positive odd integer.');
    end

    filterType = lower(string(filterType));

    p = inputParser;
    p.FunctionName = 'applyKernelFilter';

    addParameter(p, 'EdgeMode', 'preserve');
    addParameter(p, 'Sigma', kernelSize / 3);
    addParameter(p, 'ValidMask', []);
    addParameter(p, 'PreserveUniform', true);
    addParameter(p, 'UniformTolerance', 1e-12);
    addParameter(p, 'MinValidFraction', 0.5);

    parse(p, varargin{:});

    edgeMode = lower(string(p.Results.EdgeMode));
    sigma = p.Results.Sigma;
    validMask = p.Results.ValidMask;
    preserveUniform = p.Results.PreserveUniform;
    uniformTolerance = p.Results.UniformTolerance;
    minValidFraction = p.Results.MinValidFraction;

    if isempty(validMask)
        validMask = isfinite(data);
    else
        validMask = logical(validMask) & isfinite(data);
    end

    if ~isequal(size(validMask), size(data))
        error('ValidMask must have the same size as data.');
    end

    if ~(edgeMode == "preserve" || edgeMode == "normalize")
        error('EdgeMode must be either "preserve" or "normalize".');
    end

    minValidFraction = min(max(minValidFraction, 0), 1);

    data = double(data);

    filteredData = data;

    %% ---------------- Quick exits ----------------

    if ~any(validMask(:))
        filteredData(:) = NaN;
        return
    end

    validValues = data(validMask);

    if preserveUniform
        if max(validValues, [], 'omitnan') - min(validValues, [], 'omitnan') <= uniformTolerance
            % Preserve exactly uniform rasters.
            filteredData(~validMask) = NaN;
            return
        end
    end

    %% ---------------- Build kernel ----------------

    switch filterType

        case "mean"

            kernel = ones(kernelSize, kernelSize);
            kernel = kernel ./ sum(kernel(:));

        case "gaussian"

            if ~isscalar(sigma) || sigma <= 0
                error('Sigma must be a positive scalar.');
            end

            r = floor(kernelSize / 2);
            [xx, yy] = meshgrid(-r:r, -r:r);

            kernel = exp(-(xx.^2 + yy.^2) ./ (2 * sigma^2));
            kernel = kernel ./ sum(kernel(:));

        case "median"

            kernel = [];

        otherwise

            error('Invalid filter type. Choose "mean", "gaussian", or "median".');

    end

    %% ---------------- Filter ----------------

    switch filterType

        case {"mean", "gaussian"}

            % NaN-aware normalized convolution.
            %
            % Instead of:
            %   conv2(data, kernel, 'same')
            %
            % use:
            %   conv2(data*valid, kernel) / conv2(valid, kernel)
            %
            % This prevents zeros or NaNs outside the domain from diluting
            % valid cells.

            data0 = data;
            data0(~validMask) = 0;

            numerator = conv2(data0, kernel, 'same');
            denominator = conv2(double(validMask), kernel, 'same');

            candidate = numerator ./ max(denominator, eps);
            candidate(denominator <= eps) = NaN;

            switch edgeMode

                case "preserve"

                    % Preserve cells where the full kernel footprint is not
                    % completely valid. This avoids modifying rectangular
                    % borders and irregular watershed boundaries.
                    fullKernelCount = conv2(double(validMask), ones(kernelSize), 'same');
                    fullKernelMask = fullKernelCount == kernelSize^2;

                    replaceMask = validMask & fullKernelMask & isfinite(candidate);

                case "normalize"

                    % Smooth up to the boundary using only valid neighbors,
                    % but require a minimum fraction of valid kernel support.
                    validFraction = conv2(double(validMask), ones(kernelSize), 'same') ./ kernelSize^2;

                    replaceMask = validMask & ...
                                  validFraction >= minValidFraction & ...
                                  isfinite(candidate);

            end

            filteredData(replaceMask) = candidate(replaceMask);

        case "median"

            candidate = nanMedianFilter2D(data, validMask, kernelSize, edgeMode, minValidFraction);

            replaceMask = validMask & isfinite(candidate);

            filteredData(replaceMask) = candidate(replaceMask);

    end

    %% ---------------- Restore inactive cells ----------------

    filteredData(~validMask) = NaN;

end


function B = nanMedianFilter2D(A, validMask, kernelSize, edgeMode, minValidFraction)
% NANMEDIANFILTER2D
% NaN-aware median filter for irregular domains.
%
% This avoids medfilt2 padding artifacts and ignores NaNs.

    [ny, nx] = size(A);
    r = floor(kernelSize / 2);

    B = NaN(size(A));

    for i = 1:ny

        i1 = max(1, i - r);
        i2 = min(ny, i + r);

        for j = 1:nx

            if ~validMask(i,j)
                continue
            end

            j1 = max(1, j - r);
            j2 = min(nx, j + r);

            window = A(i1:i2, j1:j2);
            maskWindow = validMask(i1:i2, j1:j2);

            nValid = sum(maskWindow(:));
            nFull = kernelSize^2;

            switch edgeMode

                case "preserve"

                    % Only filter if the full kernel fits inside the valid
                    % domain. Otherwise preserve original value.
                    if nValid < nFull
                        B(i,j) = A(i,j);
                        continue
                    end

                case "normalize"

                    % Allow partial windows, but require enough valid cells.
                    if nValid / nFull < minValidFraction
                        B(i,j) = A(i,j);
                        continue
                    end

            end

            vals = window(maskWindow);
            vals = vals(isfinite(vals));

            if isempty(vals)
                B(i,j) = A(i,j);
            else
                B(i,j) = median(vals, 'omitnan');
            end

        end

    end

end