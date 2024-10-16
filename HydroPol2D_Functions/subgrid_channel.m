% Sub-grid approach for 1D Channels
% Developer: Marcus Nobrega, Ph.D.
% Goal: define subgrid approach for dealing with channels
%
% Input:
% h: flood depth matrix [m]
% w: river width matrix [m]
% zf: floodplain elevation matrix [m]
% zc: channel elevation matrix [m]
% dx: grid resolution [m]
% nc: channel Manning's roughness coefficient matrix [sm^-1/3]
% nf: floodplain Manning's roughness coefficient matrix [sm^-1/3] 
% Qc_prev: previous outflow for channel [m3/s]
% Qf_prev = previous floodplain outflow [m3/s]
% Qci_prev = channel floodplain outflow for NE and SE directions [m3/s]
% Qfi_prev = floodplain outflow for NE a
% g: gravity acceleration [m2/s]
% dt: time-step in [sec]
% 
% Output:
% Q: outflow [m3/s] with first entry for x and second entry for y
% similarly:
% Qc: channel outflow [m3/s]
% Qc: previous channel outflow [m3/s]
% Qf = previous floodplain outflow [m3/s]
% Qci = channel outflow for NE and SE [m3/s]
% Qfi = floodplain outflow for NE and SE [m3/s]
% C_a = river surface area [m2]

function [Q,Qc,Qf,Qci,Qfi,C_a] = subgrid_channel(h, w, zf, zc,dx, nc, nf, Qc_prev, Qf_prev,Qci_prev,Qfi_prev,g,dt,idx_rivers)

    % -------- Rivers outside of river domain follow original roughness
    nc(~idx_rivers) = nf(~idx_rivers);
    % 
    %w(~idx_rivers) = nan;

    % --------- Channel Calculations ----------- %
    yc = h + zc; % Channel head [m]
    yc(isnan(zc)) = nan;    
    nan_row = nan*ones(1,size(h,2));
    nan_col = nan*ones(size(h,1),1);

    % Water Surface Slope Calculation
    Sc = slope_function(yc,dx,nan_col,nan_row); % m

    % Effective Channel Flow Depth
    hfcflow = hf_function(yc,zc,nan_col,nan_row); % m
    
    % Channel Effective Flow Width
    wc_flow = width_function(w,nan_col,nan_row); % m

    % Channel Flow Area
    Acflow = wc_flow.*hfcflow; % m2

    % Channel surface area
    C_a = dx^2*ones(size(nf));
    % inbank_rivers = idx_rivers & h < (zf - zc); % logical mask
    inbank_rivers = h <= (zf - zc); % logical mask
    C_a(inbank_rivers) = dx*w(inbank_rivers); % m2 of surface area    

    % Channel Hydraulic Radius
    Rc = Acflow./(w + 2*hfcflow);

    % Channel Flow
    Qc_prev(:,:,3:4) = Qci_prev; % 3-4 are inclined directions
    Qc = (Qc_prev - g.*Acflow.*dt.*Sc)./(1 + (g*dt.*nc.^2.*abs(Qc_prev))./(Rc.^(4/3).*Acflow)); % m3/s
    Qc(isnan(Qc)) = 0; % nan to 0

    % -------- Floodplain Flow ------------- %
    % Floodplain depth
    hflood = max(0,h + zc -zf);

    yf = hflood + zf; % Floodplain head [m]
    
    % Effective Floodplain Depth [m]
    [hfflood] = hf_function(yf,zf,nan_col,nan_row);

    % Cells in which floodplain flow occurs
    mask = hfflood == 0;

    % Floodplain Slope
    [Sf] = slope_function(yf,dx,nan_col,nan_row);

    % Floodplain Flow
    Qf_prev(:,:,3:4) = Qfi_prev;
    Qf = ((Qf_prev - g.*hfflood.*dt.*Sf)./(1 + g*dt.*nf.^2.*abs(Qf_prev)./(hfflood.^(7/3)))).*(dx - wc_flow); % m3/s
    Qf(mask) = 0;
    Qf(isnan(Qf)) = 0; % nan to 0

    % Total Flow (Ortogonal Interfaces)
    Q = Qf(:,:,1:2) + Qc(:,:,1:2);

    % Inclined Interfaces (NE, SE)
 

    % Returning to original size
    Qci = Qc(:,:,3:4);
    Qfi = Qf(:,:,3:4);
    
    Qc = Qc(:,:,1:2);
    Qf = Qf(:,:,1:2);

    % Masks
    Q(isnan(Q)) = 0;
    Q(repmat(isnan(zf),[1,1,2])) = 0;
    Qc(isnan(Qc)) = 0;
    Qf(isnan(Qf)) = 0;

    Q(isinf(Q)) = 0;
    Qc(isinf(Qc)) = 0;
    Qf(isinf(Qf)) = 0;    

end


function [S] = slope_function(M,dx,nan_col,nan_row)
    S(:,:,1) = [(M(:,2:end) - M(:,1:(end-1)))/dx,nan_col]; % East
    S(:,:,2) = [nan_row; M(1:(end-1),:) - M(2:end,:)]/dx; % North
    S(2:end,1:(end-1),3) = (M(1:(end-1),2:end) - M(2:end,1:(end-1)))/(sqrt(2)*dx);
    S(1:(end-1),1:(end-1),4) = (M(2:end,2:(end)) - M(1:(end-1),1:(end-1)))/(sqrt(2)*dx);
end

function [hf] = hf_function(y,z,nan_col,nan_row)
    hf(:,:,1) = [max(y(:,2:end), y(:,1:(end-1))) - max(z(:,2:end), z(:,1:(end-1))), nan_col]; % right
    hf(:,:,2) = [nan_row; max(y(1:(end-1),:),y(2:end,:)) - max(z(1:(end-1),:),z(2:end,:))]; % up
    hf(2:end,1:(end-1),3) = max(y(1:(end-1),2:end), y(2:end,1:(end-1))) - max(z(1:(end-1),2:end), z(2:end,1:(end-1)));
    hf(1:(end-1),1:(end-1),4) = max(y(2:end,2:(end)), y(1:(end-1),1:(end-1))) - max(z(2:end,2:(end)), z(1:(end-1),1:(end-1)));
end

function [wc] = width_function(w,nan_col,nan_row)
    wc(:,:,1) = [min(w(:,1:(end-1)),(w(:,2:(end)))), nan_col]; % Improve it later
    wc(:,:,2) = [nan_row; min(w(2:(end),:),w(1:(end-1),:))];
    wc(2:end,1:(end-1),3) = [min(w(1:(end-1),2:end),w(2:end,1:(end-1)))]; % up 
    wc(1:(end-1),1:(end-1),4) = [min(w(2:end,2:(end)),w(1:(end-1),1:(end-1)))]; % up
end