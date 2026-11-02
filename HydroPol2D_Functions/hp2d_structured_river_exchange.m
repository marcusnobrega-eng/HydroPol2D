function [head,surface_depth_mm,rate_m_s,net_volume_m3] = hp2d_structured_river_exchange( ...
        head,surface_depth_mm,bedrock,surface,soil,river_mask,dt_s,cell_area,idx_nan,this_step_recharge_m)
%HP2D_STRUCTURED_RIVER_EXCHANGE River-aquifer exchange of the legacy grid groundwater.
%   The exact two-store exchange (hp2d_river_groundwater_exchange_step) between
%   the surface stage and the groundwater head of every active river cell. The
%   groundwater the river can take excludes this step's recharge
%   (this_step_recharge_m, a depth), as in HydroPol2D-Python.
assert(isfield(soil,'riverbed_hydraulic_conductivity_m_s') && ...
    isfield(soil,'cell_riverbed_area_m2') && isfield(soil,'riverbed_thickness_m'), ...
    'HydroPol2D:MissingRiverbedParameters', ...
    'River-groundwater exchange requires riverbed conductivity, area, and thickness.');
kbed=double(soil.riverbed_hydraulic_conductivity_m_s);
area=max(double(soil.cell_riverbed_area_m2),0);
thickness=double(soil.riverbed_thickness_m);
assert(isfinite(thickness) && thickness>0 && all(kbed(~idx_nan)>=0), ...
    'HydroPol2D:InvalidRiverbedParameters', ...
    'Riverbed conductivity must be nonnegative and thickness positive.');
active=logical(river_mask) & area>0 & ~idx_nan;
stage=surface+max(surface_depth_mm,0)./1000;
sy=local_groundwater_specific_yield(soil);
if nargin < 10 || isempty(this_step_recharge_m), this_step_recharge_m = 0; end
this_step_recharge_m(~isfinite(this_step_recharge_m)) = 0;
groundwater_available=max(max(head-bedrock,0).*sy-max(this_step_recharge_m,0),0).*cell_area;
surface_available=max(surface_depth_mm,0)./1000.*cell_area;
[volume,~]=hp2d_river_groundwater_exchange_step(kbed,thickness,area,head,stage, ...
    groundwater_available,surface_available,sy.*cell_area,cell_area,dt_s,active);
head=head-volume./max(sy.*cell_area,eps);
surface_depth_mm=surface_depth_mm+volume./max(cell_area,eps).*1000;
surface_depth_mm=max(surface_depth_mm,0);
rate_m_s=volume./max(cell_area,eps)./dt_s;
rate_m_s(idx_nan)=nan;
net_volume_m3=nansum(volume,'all');
end

function Sy_local = local_groundwater_specific_yield(Soil_Properties)
if isfield(Soil_Properties, 'Sy') && ~isempty(Soil_Properties.Sy)
    Sy_local = Soil_Properties.Sy;
else
    Sy_local = Soil_Properties.theta_sat - Soil_Properties.theta_r;
end
Sy_local = max(Sy_local, 1e-6);
end
