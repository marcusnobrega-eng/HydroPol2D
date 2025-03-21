function [flowAcc, flowDir] = D4_FlowAccumulation(DEM)

% Fill sinks in a DEM using only 4-directional flow rules (N, S, E, W)

[m, n] = size(DEM);
DEM_filled = DEM;

% Create a border frame with large values to act as boundary
border_val = max(DEM(:)) + 1000;
DEM_filled([1 end], :) = border_val;
DEM_filled(:, [1 end]) = border_val;

% Define 4D neighbor offsets
offsets = [0 1; 1 0; 0 -1; -1 0];

changed = true;
while changed
    changed = false;
    for i = 2:m-1
        for j = 2:n-1
            cell_val = DEM_filled(i, j);
            % Skip boundary cells
            if cell_val >= border_val
                continue;
            end

            % Find minimum elevation from 4D neighbors
            neighbor_vals = zeros(1,4);
            for k = 1:4
                di = offsets(k,1);
                dj = offsets(k,2);
                neighbor_vals(k) = DEM_filled(i+di, j+dj);
            end

            min_neighbor = min(neighbor_vals);

            % If current cell is lower than its neighbors, fill it
            if cell_val < min_neighbor
                DEM_filled(i, j) = min_neighbor;
                changed = true;
            end
        end
    end
end

DEM = DEM_filled;

[m, n] = size(DEM);
flowDir = zeros(m, n);   % Flow direction matrix

% Define offsets for 8 directions (E, SE, S, SW, W, NW, N, NE)
offsets8 = [0 1; 1 1; 1 0; 1 -1; 0 -1; -1 -1; -1 0; -1 1];

% Step 1: Calculate approximated flow direction (projected from D8 to D4)
for i = 2:m-1
    for j = 2:n-1
        slopes = zeros(1,8);
        for d = 1:8
            ni = i + offsets8(d,1);
            nj = j + offsets8(d,2);
            slopes(d) = DEM(i,j) - DEM(ni,nj);
        end

        [maxSlope, idx] = max(slopes);
        
        if maxSlope <= 0
            flowDir(i,j) = 0;  % pit
        else
            % If direction is diagonal, project to nearest cardinal direction
            if mod(idx, 2) == 0 % diagonal index (2, 4, 6, 8)
                % Diagonal flow direction handling: First, handle horizontal (East/West), then vertical (North/South)
                
                if idx == 2  % North-East (NE)
                    cardinalDirs = [1, 3]; % First East (1), then North (3)
                elseif idx == 4  % South-East (SE)
                    cardinalDirs = [1, 7]; % First East (1), then South (7)
                elseif idx == 6  % South-West (SW)
                    cardinalDirs = [5, 7]; % First West (5), then South (7)
                elseif idx == 8  % North-West (NW)
                    cardinalDirs = [5, 3]; % First West (5), then North (3)
                end
                
                % Update flowDir to allow two-way flow for diagonal directions
                % Find the strongest slope in the chosen directions
                [~, cidx] = max(slopes(cardinalDirs));
                finalDir = cardinalDirs(cidx); % Assign the best direction based on the strongest slope

                % Map finalDir to D4 codes: E=1, S=2, W=3, N=4
                switch finalDir
                    case 1 % E
                        flowDir(i,j) = 1; % East
                        % Flow also to North
                        flowDir(i-1,j) = 4; % North (N)
                    case 3 % S
                        flowDir(i,j) = 2; % South
                        % Flow also to East
                        flowDir(i,j+1) = 1; % East (E)
                    case 5 % W
                        flowDir(i,j) = 3; % West
                        % Flow also to South
                        flowDir(i+1,j) = 2; % South (S)
                    case 7 % N
                        flowDir(i,j) = 4; % North
                        % Flow also to West
                        flowDir(i,j-1) = 3; % West (W)
                    otherwise
                        flowDir(i,j) = 0;
                end
            else
                finalDir = idx; % For non-diagonal directions, use the current index
                % Map finalDir to D4 codes: E=1, S=2, W=3, N=4
                switch finalDir
                    case 1 % E
                        flowDir(i,j) = 1;
                    case 3 % S
                        flowDir(i,j) = 2;
                    case 5 % W
                        flowDir(i,j) = 3;
                    case 7 % N
                        flowDir(i,j) = 4;
                    otherwise
                        flowDir(i,j) = 0;
                end
            end
        end
    end
end


% Step 2: Calculate flow accumulation (same topological approach)
flowAcc = zeros(m, n);
inflowCount = zeros(m, n);

for i = 2:m-1
    for j = 2:n-1
        d = flowDir(i, j);
        if d ~= 0
            [di, dj] = directionOffset(d);
            inflowCount(i+di, j+dj) = inflowCount(i+di, j+dj) + 1;
        end
    end
end

queue = zeros(m*n, 2);
head = 1; tail = 0;

for i = 2:m-1
    for j = 2:n-1
        if inflowCount(i, j) == 0
            tail = tail + 1;
            queue(tail, :) = [i, j];
            flowAcc(i, j) = 1;
        end
    end
end

while head <= tail
    i = queue(head, 1);
    j = queue(head, 2);
    head = head + 1;
    d = flowDir(i, j);
    if d ~= 0
        [di, dj] = directionOffset(d);
        ni = i + di;
        nj = j + dj;
        flowAcc(ni, nj) = flowAcc(ni, nj) + flowAcc(i, j);
        inflowCount(ni, nj) = inflowCount(ni, nj) - 1;
        if inflowCount(ni, nj) == 0
            tail = tail + 1;
            queue(tail, :) = [ni, nj];
        end
    end
end
end

function [di, dj] = directionOffset(d)
    switch d
        case 1
            di = 0; dj = 1;   % East
        case 2
            di = 1; dj = 0;   % South
        case 3
            di = 0; dj = -1;  % West
        case 4
            di = -1; dj = 0;  % North
        otherwise
            di = 0; dj = 0;   % Pit
    end
end
