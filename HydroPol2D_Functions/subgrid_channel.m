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

function [Q,Qc,Qf,Qci,Qfi] = subgrid_channel(h, w, zf, zc,dx, nc, nf, Qc_prev, Qf_prev,Qci_prev,Qfi_prev,g,dt,idx_rivers)

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

    % IX = find(idx_rivers); % River Cells
    % [row_river, col_river] = find(idx_rivers);    
    % 
    % % Right, % Up, % Left, % Down % NE, NW, SW, SE
    % dir_row = [0 -1 0 1 -1 -1 1 1];
    % dir_col = [1 0 -1 0 1 -1 -1 1];
    % 
    % for i = 1:8
    %     dir_y = dir_row(i); % y direction
    %     dir_x = dir_col(i); % x direction
    %     IX_n(:,i) = (row_river + dir_y) + (col_river + dir_x)*ny; % Index of neighbor cells
    %     hfcflow(:,i) = max(yc(IX), yc(IX_n(:,i))) - max(zc(IX), zc(IX_n(:,i))); % Channel flow depth [m]
    %     wc_flow(:,i) = max(w(IX), w(IX_n(:,i))); % Channel width [m]
    %     hfflood(:,i) = max(0,h(IX) + zc(IX) - zf(IX)); % Floodplain depth [m]
    %     mask = hfflood == 0; % Mask for cells with no floodplain
    %     if i < 5
    %         Sc(:,i) = (yc(IX) - yc(IX_n(:,i)))/dx; % Channel slope
    %         Sf(:,i) = (y(IX) - y(IX_n(:,i)))/dx; % Floodplain slope
    % 
    %     else
    %         Sc(:,i) = (yc(IX) - yc(IX_n(:,i)))/(sqrt(2)*dx);
    %         Sf(:,i) = (y(IX) - y(IX_n(:,i)))/(sqrt(2)*dx); % Floodplain slope
    %     end
    % end    

    % IDXs
    % % Sc
    % Sc(:,:,1) = [(yc(:,2:end) - yc(:,1:(end-1)))/dx,nan_col]; % East
    % Sc(:,:,2) = [nan_row; yc(1:(end-1),:) - yc(2:end,:)]/dx; % North
    % Sc(2:end,1:(end-1),3) = [yc(1:(end-1),2:end) - yc(2:end,1:(end-1))]/(sqrt(2)*dx);
    % Sc(1:(end-1),1:(end-1),4) = [yc(2:end,2:(end)) - yc(1:(end-1),1:(end-1))]/(sqrt(2)*dx);

    % Water Surface Slope Calculation
    Sc = slope_function(yc,dx,nan_col,nan_row);

    % % hcflow
    % hfcflow(:,:,1) = [max(yc(:,2:nx), yc(:,1:(nx-1))) - max(zc(:,2:nx), zc(:,1:(nx-1))), nan_col]; % right
    % hfcflow(:,:,2) = [nan_row; max(yc(1:(end-1),:),yc(2:end,:)) - max(zc(1:(end-1),:),zc(2:end,:))]; % up 
    % hfcflow(2:end,1:(end-1),3) = [max(yc(1:(end-1),2:end), yc(2:end,1:(end-1))) - max(zc(1:(end-1),2:end), zc(2:end,1:(end-1)))];
    % hfcflow(1:(end-1),1:(end-1),4) = [max(yc(1:(end-1),2:end), yc(1:(end-1),1:(end-1))) - max(zc(1:(end-1),2:end), zc(2:end,1:(end-1)))];
    
    % Effective Channel Flow Depth
    hfcflow = hf_function(yc,zc,nan_col,nan_row);
    % 
    % % Width flow calculation
    % wc_flow(:,:,1) = [max(w(:,1:(end-1)),(w(:,2:(end)))), nan_col]; % Improve it later
    % wc_flow(:,:,2) = [nan_row; max(w(2:(end),:),w(1:(end-1),:))];
    % wc_flow(2:end,1:(end-1),3) = [max(w(1:(end-1),2:end),w(2:end,1:(end-1)))]; % up 
    % wc_flow(1:(end-1),1:(end-1),4) = [max(w(1:(end-1),2:end),w(1:(end-1),1:(end-1)))]; % up 
    
    % Channel Effective Flow Width
    wc_flow = width_function(w,nan_col,nan_row);

    % Channel Flow Area
    Acflow = wc_flow.*hfcflow;

    % Channel Hydraulic Radius
    Rc = Acflow./(wc_flow + 2*hfcflow);

    % Channel Flow
    Qc_prev(:,:,3:4) = Qci_prev;
    Qc = (Qc_prev - g.*Acflow.*dt.*Sc)./(1 + g*dt.*nc.^2.*abs(Qc_prev)./(Rc.^(4/3).*Acflow)); % m3/s
    Qc(isnan(Qc)) = 0;
    % Qc(repmat(~idx_rivers,[1 1 4])) = 0; % Attention    

    % -------- Floodplain Flow ------------- %
    % Floodplain depth
    hflood = max(0,h + zc -zf);

    yf = hflood + zf; % Floodplain head [m]
    
    % Effective Floodplain Depth [m]
    [hfflood] = hf_function(yf,zf,nan_col,nan_row);

    % hfflood(:,:,1) = [max(yf(:,2:nx), yf(:,1:(nx-1))) - max(zf(:,2:nx), zf(:,1:(nx-1))), nan_col]; % right
    % hfflood(:,:,2) = [nan_row; max(yf(1:(end-1),:),yf(2:end,:)) - max(zf(1:(end-1),:),zf(2:end,:))]; % up 
    % hfflood(2:end,1:(end-1),3) = [max(yf(1:(end-1),2:end), yf(2:end,1:(end-1))) - max(zf(1:(end-1),2:end), zf(2:end,1:(end-1)))];
    % hfflood(1:(end-1),1:(end-1),4) = [max(yf(1:(end-1),2:end), yf(1:(end-1),1:(end-1))) - max(zf(1:(end-1),2:end), zf(1:(end-1),1:(end-1)))];
    
    % Cells in which floodplain flow occurs
    mask = hfflood == 0;
    
    % Sf(:,:,1) = [(yf(:,2:end) - yf(:,1:(end-1)))/dx,nan_col]; % East
    % Sf(:,:,2) = [nan_row; yf(1:(end-1),:) - yf(2:end,:)]/dx; % North
    % Sf(2:end,1:(end-1),3) = [yf(1:(end-1),2:end) - yf(2:end,1:(end-1))]/(sqrt(2)*dx);
    % Sf(1:(end-1),1:(end-1),4) = [yf(2:end,1:(end-1)) - yf(1:(end-1),2:end)]/(sqrt(2)*dx);

    % Floodplain Slope
    [Sf] = slope_function(yf,dx,nan_col,nan_row);

    % Floodplain Flow
    Qf_prev(:,:,3:4) = Qfi_prev;
    Qf = ((Qf_prev - g.*hfflood.*dt.*Sf)./(1 + g*dt.*nf.^2.*abs(Qf_prev)./(hfflood.^(7/3)))).*(dx - wc_flow); % m3/s
    % Qf(repmat(~idx_rivers,[1 1 4])) = 0; % Attention
    Qf(mask) = 0;
    Qc(isnan(Qf)) = 0;

    % Total Flow (Ortogonal Interfaces)
    Q = Qf(:,:,1:2) + Qc(:,:,1:2);

    % Inclined Interfaces (NE, SE)
    % Q(:,:,1) = Q(:,:,1) + sqrt(2)/2*(Qc(:,:,3) + Qc(:,:,4));
    % Q(:,:,2) = Q(:,:,2) + sqrt(2)/2*(Qc(:,:,3) - Qc(:,:,4));

    % Returning to original size
    Qci = Qc(:,:,3:4);
    Qfi = Qf(:,:,3:4);
    
    Qc = Qc(:,:,1:2);
    Qf = Qf(:,:,1:2);

    % Masks
    Q(isnan(Q)) = 0;
    Qc(isnan(Qc)) = 0;
    Qf(isnan(Qf)) = 0;

    Q(isinf(Q)) = 0;
    Qc(isinf(Qc)) = 0;
    Qf(isinf(Qf)) = 0;    

    % Q(repmat(~idx_rivers,[1,1,2])) = 0;
    % Qc(repmat(~idx_rivers,[1,1,2])) = 0;
    % Qf(repmat(~idx_rivers,[1,1,2])) = 0;
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
    wc(:,:,1) = [max(w(:,1:(end-1)),(w(:,2:(end)))), nan_col]; % Improve it later
    wc(:,:,2) = [nan_row; max(w(2:(end),:),w(1:(end-1),:))];
    wc(2:end,1:(end-1),3) = [max(w(1:(end-1),2:end),w(2:end,1:(end-1)))]; % up 
    wc(1:(end-1),1:(end-1),4) = [max(w(2:end,2:(end)),w(1:(end-1),1:(end-1)))]; % up
end