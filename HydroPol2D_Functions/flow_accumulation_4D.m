function [flow_accum, flow_dir] = flow_accumulation_4D(dem)
    % FLOW_ACCUMULATION_4D computes flow accumulation using a 4-direction flow model
    % that ensures diagonal flows are split into cardinal directions.
    %
    % INPUT:
    %   dem - Digital Elevation Model (2D matrix)
    % OUTPUT:
    %   flow_accum - Flow accumulation matrix (number of cells draining to each cell)
    %   flow_dir - Flow direction matrix (1=N, 2=E, 3=S, 4=W)
    
    [rows, cols] = size(dem);
    flow_accum = zeros(rows, cols);  % Accumulation starts at zero
    flow_dir = zeros(rows, cols);   % Flow direction map
    
    % Compute flow directions (4D: N, E, S, W)
    for i = 2:rows-1
        for j = 2:cols-1
            % Get elevation differences
            dN = dem(i, j) - dem(i-1, j);
            dE = dem(i, j) - dem(i, j+1);
            dS = dem(i, j) - dem(i+1, j);
            dW = dem(i, j) - dem(i, j-1);
            
            % Compute potential flow directions
            slopes = [dN, dE, dS, dW];
            [maxSlope, dir] = max(slopes);
            
            if maxSlope > 0
                flow_dir(i, j) = dir; % Assign primary flow direction
                
                % Update flow accumulation (number of contributing cells)
                if dir == 1  % N
                    flow_accum(i-1, j) = flow_accum(i-1, j) + 1;
                elseif dir == 2  % E
                    flow_accum(i, j+1) = flow_accum(i, j+1) + 1;
                elseif dir == 3  % S
                    flow_accum(i+1, j) = flow_accum(i+1, j) + 1;
                elseif dir == 4  % W
                    flow_accum(i, j-1) = flow_accum(i, j-1) + 1;
                end
                
                % Handle diagonal flow splitting
                if dir == 1  % N
                    if dE > 0, flow_accum(i-1, j) = flow_accum(i-1, j) + 1; flow_accum(i, j+1) = flow_accum(i, j+1) + 1; end
                    if dW > 0, flow_accum(i-1, j) = flow_accum(i-1, j) + 1; flow_accum(i, j-1) = flow_accum(i, j-1) + 1; end
                elseif dir == 2  % E
                    if dN > 0, flow_accum(i, j+1) = flow_accum(i, j+1) + 1; flow_accum(i-1, j) = flow_accum(i-1, j) + 1; end
                    if dS > 0, flow_accum(i, j+1) = flow_accum(i, j+1) + 1; flow_accum(i+1, j) = flow_accum(i+1, j) + 1; end
                elseif dir == 3  % S
                    if dE > 0, flow_accum(i+1, j) = flow_accum(i+1, j) + 1; flow_accum(i, j+1) = flow_accum(i, j+1) + 1; end
                    if dW > 0, flow_accum(i+1, j) = flow_accum(i+1, j) + 1; flow_accum(i, j-1) = flow_accum(i, j-1) + 1; end
                elseif dir == 4  % W
                    if dN > 0, flow_accum(i, j-1) = flow_accum(i, j-1) + 1; flow_accum(i-1, j) = flow_accum(i-1, j) + 1; end
                    if dS > 0, flow_accum(i, j-1) = flow_accum(i, j-1) + 1; flow_accum(i+1, j) = flow_accum(i+1, j) + 1; end
                end
            end
        end
    end
    
    % Ensure self-contributing cells count as 1
    flow_accum = flow_accum + 1;
end
