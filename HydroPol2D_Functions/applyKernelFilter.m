function filteredData = applyKernelFilter(data, kernelSize, filterType)
    % APPLYKERNELFILTER - Applies a kernel filter to reduce noise in data.
    %
    % Usage:
    %   filteredData = applyKernelFilter(data, kernelSize, filterType)
    %
    % Inputs:
    %   data       - 2D matrix representing the noisy data
    %   kernelSize - Size of the kernel (e.g., 3 for a 3x3 kernel)
    %   filterType - Type of filter ('mean', 'gaussian', 'median')
    %
    % Output:
    %   filteredData - Smoothed data after applying the kernel
    
    % Ensure kernelSize is an odd number
    if mod(kernelSize, 2) == 0
        error('Kernel size must be an odd number');
    end
    
    % Select the filter type
    switch lower(filterType)
        case 'mean'
            % Create a mean filter (box filter)
            kernel = ones(kernelSize) / (kernelSize^2);
            filteredData = conv2(data, kernel, 'same');
        
        case 'gaussian'
            % Create a Gaussian filter
            sigma = kernelSize / 3; % Standard deviation
            kernel = fspecial('gaussian', kernelSize, sigma);
            filteredData = conv2(data, kernel, 'same');
        
        case 'median'
            % Apply median filtering
            filteredData = medfilt2(data, [kernelSize kernelSize]);
        
        otherwise
            error('Invalid filter type. Choose "mean", "gaussian", or "median".');
    end
end