function [flow_mask_4D, flow_width_4D, flow_depth_4D] = Transform8Dto4DFlow(flow_mask_8D, flow_width_8D, flow_depth_8D)
    % CONVERT_8D_TO_4D Transforms an 8-directional river mask to a 4-directional one
    % ensuring diagonal river cells are extended into horizontal and vertical components.
    % Also updates flow width and depth accordingly.
    %
    % INPUT:
    %   flow_mask_8D - A binary matrix where 1 indicates a river cell, 0 is non-river
    %   flow_width_8D - Matrix of flow widths corresponding to 8D river mask
    %   flow_depth_8D - Matrix of flow depths corresponding to 8D river mask
    % OUTPUT:
    %   flow_mask_4D - A modified binary mask adhering to 4D connectivity
    %   flow_width_4D - Adjusted flow width matrix for 4D connectivity
    %   flow_depth_4D - Adjusted flow depth matrix for 4D connectivity
    
    [rows, cols] = size(flow_mask_8D);
    flow_mask_4D = flow_mask_8D; % Start with the same river mask
    flow_width_4D = flow_width_8D; % Copy width matrix
    flow_depth_4D = flow_depth_8D; % Copy depth matrix
    
    % Process river cells to ensure 4D connectivity and adjust width and depth
    for i = 2:rows-1
        for j = 2:cols-1
            if flow_mask_8D(i, j) == 1
                % Check diagonal neighbors and distribute width and depth
                if flow_mask_8D(i-1, j+1) == 1 % NE
                    flow_mask_4D(i-1, j) = 1; % Ensure N
                    flow_mask_4D(i, j+1) = 1; % Ensure E
                    flow_width_4D(i-1, j) = max(flow_width_4D(i-1, j), flow_width_8D(i, j) );
                    flow_width_4D(i, j+1) = max(flow_width_4D(i, j+1), flow_width_8D(i, j) );
                    flow_depth_4D(i-1, j) = max(flow_depth_4D(i-1, j), flow_depth_8D(i, j));
                    flow_depth_4D(i, j+1) = max(flow_depth_4D(i, j+1), flow_depth_8D(i, j));
                end
                if flow_mask_8D(i+1, j+1) == 1 % SE
                    flow_mask_4D(i+1, j) = 1; % Ensure S
                    flow_mask_4D(i, j+1) = 1; % Ensure E
                    flow_width_4D(i+1, j) = max(flow_width_4D(i+1, j), flow_width_8D(i, j) );
                    flow_width_4D(i, j+1) = max(flow_width_4D(i, j+1), flow_width_8D(i, j) );
                    flow_depth_4D(i+1, j) = max(flow_depth_4D(i+1, j), flow_depth_8D(i, j));
                    flow_depth_4D(i, j+1) = max(flow_depth_4D(i, j+1), flow_depth_8D(i, j));
                end
                if flow_mask_8D(i+1, j-1) == 1 % SW
                    flow_mask_4D(i+1, j) = 1; % Ensure S
                    flow_mask_4D(i, j-1) = 1; % Ensure W
                    flow_width_4D(i+1, j) = max(flow_width_4D(i+1, j), flow_width_8D(i, j) );
                    flow_width_4D(i, j-1) = max(flow_width_4D(i, j-1), flow_width_8D(i, j) );
                    flow_depth_4D(i+1, j) = max(flow_depth_4D(i+1, j), flow_depth_8D(i, j));
                    flow_depth_4D(i, j-1) = max(flow_depth_4D(i, j-1), flow_depth_8D(i, j));
                end
                if flow_mask_8D(i-1, j-1) == 1 % NW
                    flow_mask_4D(i-1, j) = 1; % Ensure N
                    flow_mask_4D(i, j-1) = 1; % Ensure W
                    flow_width_4D(i-1, j) = max(flow_width_4D(i-1, j), flow_width_8D(i, j) );
                    flow_width_4D(i, j-1) = max(flow_width_4D(i, j-1), flow_width_8D(i, j) );
                    flow_depth_4D(i-1, j) = max(flow_depth_4D(i-1, j), flow_depth_8D(i, j));
                    flow_depth_4D(i, j-1) = max(flow_depth_4D(i, j-1), flow_depth_8D(i, j));
                end
            end
        end
    end
end
