function flow_mask_4D = convert_8D_to_4D(flow_mask_8D)
    % CONVERT_8D_TO_4D Transforms an 8-directional flow mask to a 4-directional one
    % ensuring diagonal flows (NE, NW, SE, SW) are split into horizontal and vertical components.
    %
    % INPUT:
    %   flow_mask_8D - A matrix representing flow directions in 8D (1-8 encoding)
    % OUTPUT:
    %   flow_mask_4D - A modified flow direction matrix adhering to 4D (N, E, S, W)
    
    [rows, cols] = size(flow_mask_8D);
    flow_mask_4D = zeros(rows, cols);
    
    % Mapping from 8D to 4D
    for i = 2:rows-1
        for j = 2:cols-1
            switch flow_mask_8D(i, j)
                case 1 % N
                    flow_mask_4D(i, j) = 1; 
                case 2 % NE (split into N and E)
                    flow_mask_4D(i, j) = 1; % N
                    flow_mask_4D(i, j) = 2; % E
                case 3 % E
                    flow_mask_4D(i, j) = 2;
                case 4 % SE (split into S and E)
                    flow_mask_4D(i, j) = 3; % S
                    flow_mask_4D(i, j) = 2; % E
                case 5 % S
                    flow_mask_4D(i, j) = 3;
                case 6 % SW (split into S and W)
                    flow_mask_4D(i, j) = 3; % S
                    flow_mask_4D(i, j) = 4; % W
                case 7 % W
                    flow_mask_4D(i, j) = 4;
                case 8 % NW (split into N and W)
                    flow_mask_4D(i, j) = 1; % N
                    flow_mask_4D(i, j) = 4; % W
            end
        end
    end
end