% Given cell
x_ij = round(size(d_p,2)/2);
y_ij = round(size(d_p,1)/2);
z_cell = z(y_ij,x_ij);
y_cell = y(y_ij,x_ij);

z_right = z(y_ij,x_ij+1);
z_left = z(y_ij,x_ij-1);
z_up = z(y_ij-1,x_ij);
z_down = z(y_ij+1,x_ij);

y_right = y(y_ij,x_ij+1);
y_left = y(y_ij,x_ij-1);
y_up = y(y_ij-1,x_ij);
y_down = y(y_ij+1,x_ij);

% Slopes
Slope_left = (y_cell - y_left)/Resolution;
if Slope_left ~= matrix_store(y_ij,x_ij,1) && ~isnan(Slope_left)
    ttt = 1;
end
Slope_right = (y_right - y_cell)/Resolution; 
if Slope_right ~= matrix_store(y_ij,x_ij,2) && ~isnan(Slope_right)
    ttt = 1;
end
Slope_up = (y_up - y_cell)/Resolution;
if Slope_up ~= matrix_store(y_ij,x_ij,3) && ~isnan(Slope_up)
    ttt = 1;
end
Slope_down = (y_cell - y_down)/Resolution;
if Slope_down ~= matrix_store(y_ij,x_ij,4) && ~isnan(Slope_down)
    ttt = 1;
end
% Hf
Hf_left = max(y_cell,y_left) - max(z_cell,z_left);
if Hf_left ~= Hf(y_ij,x_ij,1)
    ttt = 1;
end
Hf_right = max(y_right,y_cell) - max(z_right,z_cell);
if Hf_right ~= Hf(y_ij,x_ij,2)
    ttt = 1;
end
Hf_up = max(y_up, y_cell) - max(z_up,z_cell);
if Hf_up ~= Hf(y_ij,x_ij,3)
    ttt = 1;
end
Hf_down = max(y_cell,y_down) - max(z_cell,z_down);
if Hf_down ~= Hf(y_ij,x_ij,4)
    ttt = 1;
end

%% Checking Fluxes
% x-x flux
x_flux = outflow(:,:,2) - outflow(:,:,1);
y_flux = outflow(:,:,3) - outflow(:,:,4);
outlet_flux = 0;

net_flux = sum(sum(x_flux + y_flux - outlet_flux));
net_volume = net_flux*Resolution*dt; % m3

