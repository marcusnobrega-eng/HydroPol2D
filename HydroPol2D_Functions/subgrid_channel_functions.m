% Sub-grid approach for 1D Channels
% Developer: Marcus Nobrega, Ph.D.
% Goal: define subgrid approach for dealing with channels
%
% Input:
% h: flood depth matrix [m]
% w: river width matrix [m]
% zc: channel elevation matrix [m]
% dx: grid resolution [m]
% nc: channel Manning's roughness coefficient matrix [sm^-1/3]
% Qc_prev: previous outflow for channel [m3/s]
% g: gravity acceleration [m2/s]
% dt: time-step in [sec]
% 
% Output:
% q: outflow [m2/s] with first entry for x and second entry for y
% similarly:
% C_a = river surface area [m2]

function [q,C_a] = subgrid_channel_functions(h, w, zc,dx, nc, Qc_prev,Qci_prev, g, dt, idx_rivers, Subgrid_Properties)

    % --------- Channel Calculations ----------- %
    yc = h + zc; % Channel head [m]
    yc(isnan(zc)) = nan;

    % Cells that are not river change to nan
    % yc(~idx_rivers) = inf;
    % h(~idx_rivers) = inf;
    % w(~idx_rivers) = inf;
    % nc(~idx_rivers) = inf;
    
    nan_row = nan*ones(1,size(h,2));
    nan_col = nan*ones(size(h,1),1);

    % Water Surface Slope Calculation
    Sc = slope_function(yc,dx,nan_col,nan_row);
  
    % Effective Channel Flow Depth
    hfcflow = hf_function(yc,zc,nan_col,nan_row);

    % Cells in which channel flow occurs
    mask = hfcflow <= 0;    
        
    % Compute flow width for all cells with their respective h_targets
    h(isnan(h)) = 0;

    [Ac(:,:,1), Ac(:,:,2), hydraulic_radius(:,:,1), hydraulic_radius(:,:,2)] = Compute_Subgrid_Properties_Polynomialwise(Subgrid_Properties.A_east_spline, ...
                                                   Subgrid_Properties.A_north_spline, Subgrid_Properties.Rh_east_spline, Subgrid_Properties.Rh_north_spline, ...
                                                   h);


    % Channel Effective Flow Width
    wc_flow = width_function(w,nan_col,nan_row);

    % Channel Flow Area
    % Acflow = wc_flow.*hfcflow;
    Acflow = area_function(Ac,nan_col,nan_row);

    % Channel surface area
    % C_a = dx*w;

    C_a = dx*dx*ones(size(zc));
    C_a(isnan(h)) = nan;

    % Channel Hydraulic Radius
    Rc = Rh_function(hydraulic_radius,nan_col,nan_row);

    % Channel Hydraulic Radius
    % Rc = cellfun(@(func, h) func(h) * ~isempty(func), Subgrid_Properties.hydraulic_radius, num2cell(h));

    % Channel Flow
    Qc_prev(:,:,3:4) = Qci_prev;
    Qc = (Qc_prev - g.*Acflow.*dt.*Sc)./(1 + g*dt.*nc.^2.*abs(Qc_prev)./(Rc.^(4/3).*Acflow)); % m3/s
    Qc(mask) = 0;
    Qc(isnan(Qc)) = 0;
    Qc(isinf(Qc)) = 0;

    % Total Flow (Ortogonal Interfaces)
    Q = Qc(:,:,1:2); Q(:,:,3:4) = zeros(size(h,1),size(h,2),2);
    
    q = Q./wc_flow; %  flow per unit width
    q(:,:,3:4) = [];

end


function [S] = slope_function(M,dx,nan_col,nan_row)
    S(:,:,1) = [(M(:,2:end) - M(:,1:(end-1)))/dx,nan_col]; % East
    S(:,:,2) = [nan_row; M(1:(end-1),:) - M(2:end,:)]/dx; % North
    S(2:end,1:(end-1),3) = (M(1:(end-1),2:end) - M(2:end,1:(end-1)))/(sqrt(2)*dx); % NE
    S(1:(end-1),1:(end-1),4) = (M(2:end,2:(end)) - M(1:(end-1),1:(end-1)))/(sqrt(2)*dx); % SE
end

function [hf] = hf_function(y,z,nan_col,nan_row)
    hf(:,:,1) = [max(y(:,2:end), y(:,1:(end-1))) - max(z(:,2:end), z(:,1:(end-1))), nan_col]; % East
    hf(:,:,2) = [nan_row; max(y(1:(end-1),:),y(2:end,:)) - max(z(1:(end-1),:),z(2:end,:))]; % North
    hf(2:end,1:(end-1),3) = max(y(1:(end-1),2:end), y(2:end,1:(end-1))) - max(z(1:(end-1),2:end), z(2:end,1:(end-1))); % NE
    hf(1:(end-1),1:(end-1),4) = max(y(2:end,2:(end)), y(1:(end-1),1:(end-1))) - max(z(2:end,2:(end)), z(1:(end-1),1:(end-1))); % SE
end

function [wc] = width_function(w,nan_col,nan_row)
    wc(:,:,1) = [min(w(:,1:(end-1)),(w(:,2:(end)))), nan_col]; % East
    wc(:,:,2) = [nan_row; min(w(2:(end),:),w(1:(end-1),:))]; % North
    wc(2:end,1:(end-1),3) = [min(w(1:(end-1),2:end),w(2:end,1:(end-1)))]; % NE 
    wc(1:(end-1),1:(end-1),4) = [min(w(2:end,2:(end)),w(1:(end-1),1:(end-1)))]; % SE
    % river_width = repmat(w,[1 1 4]);
    % wc(wc == 0) = river_width(wc == 0);
    % wc(repmat(~idx_rivers,[1 1 4])) = 0;
end

function [wc] = area_function(w,nan_col,nan_row)
    wc(:,:,1) = [min(w(:,1:(end-1),1),(w(:,2:(end),1))), nan_col]; % East
    wc(:,:,2) = [nan_row; min(w(2:(end),:,2),w(1:(end-1),:,2))]; % North
    wc(2:end,1:(end-1),3) = [min(w(1:(end-1),2:end,1),w(2:end,1:(end-1),2))]; % NE 
    wc(1:(end-1),1:(end-1),4) = [min(w(2:end,2:(end),2),w(1:(end-1),1:(end-1),1))]; % SE
    % river_width = repmat(w,[1 1 4]);
    % wc(wc == 0) = river_width(wc == 0);
    % wc(repmat(~idx_rivers,[1 1 4])) = 0;
end

function [wc] = Rh_function(w,nan_col,nan_row)
    wc(:,:,1) = [min(w(:,1:(end-1),1),(w(:,2:(end),1))), nan_col]; % East
    wc(:,:,2) = [nan_row; min(w(2:(end),:,2),w(1:(end-1),:,2))]; % North
    wc(2:end,1:(end-1),3) = [min(w(1:(end-1),2:end,1),w(2:end,1:(end-1),2))]; % NE 
    wc(1:(end-1),1:(end-1),4) = [min(w(2:end,2:(end),2),w(1:(end-1),1:(end-1),1))]; % SE
    % river_width = repmat(w,[1 1 4]);
    % wc(wc == 0) = river_width(wc == 0);
    % wc(repmat(~idx_rivers,[1 1 4])) = 0;
end

function [A, W, Rh] = Compute_Subgrid_Properties_Splinewise(A_splines, W_splines, Rh_splines, H)
% Computes A, W, and Rh using loops to evaluate splines for each depth in H.
% A_splines, W_splines, Rh_splines: [nrows x ncols] cell arrays of spline structures.
% H: [nrows x ncols] depth matrix.

[nrows, ncols] = size(H);
A = NaN(nrows, ncols);
W = NaN(nrows, ncols);
Rh = NaN(nrows, ncols);

% Loop over each cell in the matrix
for rowc = 1:nrows
    for colc = 1:ncols
        % Get the corresponding depth for the current cell
        depth = H(rowc, colc);

        % Check if the spline exists and the depth is valid (not NaN)
        if ~isempty(A_splines{rowc, colc}) && ~isnan(depth)
            % Evaluate the splines for area, width, and hydraulic radius at the given depth
            A(rowc, colc) = feval(A_splines{rowc, colc}, depth);
            W(rowc, colc) = feval(W_splines{rowc, colc}, depth);
            Rh(rowc, colc) = feval(Rh_splines{rowc, colc}, depth);
        end
    end
end

% Ensure physical constraints (no values below the threshold)
A = max(A, 1e-12);
W = max(W, 1e-12);
Rh = max(Rh, 1e-12);

end

function [Area_east, Area_north, Rh_east, Rh_north] = Compute_Subgrid_Properties_Polynomialwise(A_east_coeffs, A_north_coeffs, Rh_east_coeffs, Rh_north_coeffs, Depths)
    % Computes Area_east, Area_north, Rh_east, Rh_north using stored polynomial coefficients
    % A_east_coeffs, A_north_coeffs, Rh_east_coeffs, Rh_north_coeffs: [nrows x ncols] matrices containing polynomial coefficients.
    % Depths: [nrows x ncols] matrix containing the flow depths at each grid cell.

    % Ensure that the polynomial coefficients are in the correct shape for matrix operations
    [nrows, ncols] = size(Depths);

    % Evaluate polynomials across the entire matrix using polyval
    Area_east = zeros(size(Depths)); Area_north = zeros(size(Depths)); Rh_east = zeros(size(Depths)); Rh_north = zeros(size(Depths));

    for i = 1:size(A_east_coeffs,3)
        Area_east = Area_east + A_east_coeffs(:,:,i).*Depths.^i;
        Area_north = Area_north + A_north_coeffs(:,:,i).*Depths.^i;
        Rh_east = Rh_east + Rh_east_coeffs(:,:,i).*Depths.^i;
        Rh_north = Rh_north + Rh_north_coeffs(:,:,i).*Depths.^i;
    end
    
    % Area_east = polyval(reshape(A_east_coeffs, nrows * ncols, []), repmat(Depths(:),2)); % Evaluates Area_east for all depths in a vectorized way
    % Area_north = polyval(reshape(A_north_coeffs, nrows * ncols, []), Depths(:));
    % Rh_east = polyval(reshape(Rh_east_coeffs, nrows * ncols, []), Depths(:));
    % Rh_north = polyval(reshape(Rh_north_coeffs, nrows * ncols, []), Depths(:));

    % Reshape back to the original matrix shape
    % Area_east = reshape(Area_east, nrows, ncols);
    % Area_north = reshape(Area_north, nrows, ncols);
    % Rh_east = reshape(Rh_east, nrows, ncols);
    % Rh_north = reshape(Rh_north, nrows, ncols);

    % Ensure physical constraints (no values below the threshold)
    Area_east = max(Area_east, 1e-12);
    Area_north = max(Area_north, 1e-12);
    Rh_east = max(Rh_east, 1e-12);
    Rh_north = max(Rh_north, 1e-12);
end
