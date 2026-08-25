function [surface_volume, hu, hv, edge_q, diagnostics, boundary_diagnostics] = ...
    Voronoi_Full_Momentum_Step(mesh, surface_volume, hu, hv, roughness, boundary, dt, options)
%VORONOI_FULL_MOMENTUM_STEP Explicit HLL shallow-water step on UGRID faces.
%
% This is a resolved-2D solver.  It deliberately has no Neal/subgrid
% channel representation: the conserved state is U=[h,hu,hv] per polygon.

arguments
    mesh struct
    surface_volume (:,1) double
    hu (:,1) double
    hv (:,1) double
    roughness double
    boundary struct
    dt (1,1) double {mustBePositive}
    options.gravity (1,1) double = 9.81
    options.dry_tolerance_m (1,1) double = 1e-6
    options.maximum_velocity_m_s (1,1) double {mustBePositive} = 10
end

n = mesh.n_cells;
if isscalar(roughness), roughness = repmat(roughness, n, 1); end
assert(numel(surface_volume)==n && numel(hu)==n && numel(hv)==n && numel(roughness)==n);
surface_volume = surface_volume(:); hu = hu(:); hv = hv(:); roughness = roughness(:);
roughness(~isfinite(roughness) | roughness <= 0) = 1e-6;
h = max(surface_volume ./ mesh.surface_area(:), 0);
[u, v] = velocity_components(h, hu, hv, options.dry_tolerance_m);

owner = mesh.edge_owner(:); neighbor = mesh.edge_neighbor(:);
internal = neighbor > 0;
edge_q = zeros(mesh.n_edges,1);
flux_hu_owner = zeros(mesh.n_edges,1); flux_hv_owner = zeros(mesh.n_edges,1);
flux_hu_neighbor = zeros(mesh.n_edges,1); flux_hv_neighbor = zeros(mesh.n_edges,1);

% Internal hydrostatic-reconstruction HLL fluxes.
o = owner(internal); d = neighbor(internal);
[nx, ny] = outward_normals(mesh, internal);
tx = -ny; ty = nx;
eta = mesh.surface_bed(:) + h;
bed_face = max(mesh.surface_bed(o), mesh.surface_bed(d));
hL = max(eta(o) - bed_face, 0); hR = max(eta(d) - bed_face, 0);
qnL = hL .* (u(o).*nx + v(o).*ny); qtL = hL .* (u(o).*tx + v(o).*ty);
qnR = hR .* (u(d).*nx + v(d).*ny); qtR = hR .* (u(d).*tx + v(d).*ty);
[Fh,Fqn,Fqt] = hll_normal_flux(hL,qnL,qtL,hR,qnR,qtR,options.gravity,options.dry_tolerance_m);
correctionL = 0.5*options.gravity .* (h(o).^2 - hL.^2);
correctionR = 0.5*options.gravity .* (h(d).^2 - hR.^2);
edge_q(internal) = Fh;
flux_hu_owner(internal) = (Fqn + correctionL).*nx + Fqt.*tx;
flux_hv_owner(internal) = (Fqn + correctionL).*ny + Fqt.*ty;
flux_hu_neighbor(internal) = (Fqn + correctionR).*nx + Fqt.*tx;
flux_hv_neighbor(internal) = (Fqn + correctionR).*ny + Fqt.*ty;

% Boundary ghost states and HLL fluxes.  Unspecified boundary faces are
% reflective walls, so all exterior faces are hydraulically defined.
boundary_edges = find(~internal);
if ~isempty(boundary_edges)
    [types, values] = boundary_values(mesh, boundary, boundary_edges);
    b = owner(boundary_edges);
    [bnx,bny] = outward_normals(mesh, ~internal);
    btx=-bny; bty=bnx;
    hL=h(b); qnL=hL.*(u(b).*bnx+v(b).*bny); qtL=hL.*(u(b).*btx+v(b).*bty);
    hR=hL; qnR=qnL; qtR=qtL;
    wall = types=="wall";
    qnR(wall)=-qnL(wall); % reflective normal velocity
    stage = types=="stage";
    hR(stage)=max(values(stage)-mesh.surface_bed(b(stage)),0);
    qnR(stage)=0; qtR(stage)=qtL(stage);
    normal = types=="normal_flow";
    qnR(normal)=hR(normal).^(5/3)./roughness(b(normal)).*sqrt(max(values(normal),0));
    critical = types=="critical_flow";
    qnR(critical)=hR(critical).*sqrt(options.gravity.*hR(critical));
    inflow = types=="inflow";
    qnR(inflow)=-values(inflow)./mesh.edge_length(boundary_edges(inflow));
    [Fh,Fqn,Fqt] = hll_normal_flux(hL,qnL,qtL,hR,qnR,qtR,options.gravity,options.dry_tolerance_m);
    edge_q(boundary_edges)=Fh;
    flux_hu_owner(boundary_edges)=Fqn.*bnx+Fqt.*btx;
    flux_hv_owner(boundary_edges)=Fqn.*bny+Fqt.*bty;
end

% One donor limiter covers internal and boundary mass fluxes before all
% conserved variables are updated. Positive Fh is outward from edge owner.
Q = edge_q .* mesh.edge_length(:);
out_cell = owner;
out_cell(internal & Q < 0) = neighbor(internal & Q < 0);
outgoing = accumarray(out_cell(Q>0 | (internal & Q<0)), abs(Q(Q>0 | (internal & Q<0))).*dt, [n 1], @sum, 0);
scale = min(1, surface_volume ./ max(outgoing,eps));
edge_scale=ones(mesh.n_edges,1);
edge_scale(Q>0 | (internal & Q<0))=scale(out_cell(Q>0 | (internal & Q<0)));
Q=Q.*edge_scale; edge_q=Q./mesh.edge_length(:);
flux_hu_owner=flux_hu_owner.*edge_scale; flux_hv_owner=flux_hv_owner.*edge_scale;
flux_hu_neighbor=flux_hu_neighbor.*edge_scale; flux_hv_neighbor=flux_hv_neighbor.*edge_scale;

volume_change=accumarray(owner,-Q.*dt,[n 1],@sum,0);
hu_change=accumarray(owner,-flux_hu_owner.*mesh.edge_length(:).*dt,[n 1],@sum,0);
hv_change=accumarray(owner,-flux_hv_owner.*mesh.edge_length(:).*dt,[n 1],@sum,0);
if any(internal)
    volume_change=volume_change+accumarray(neighbor(internal),Q(internal).*dt,[n 1],@sum,0);
    hu_change=hu_change+accumarray(neighbor(internal),flux_hu_neighbor(internal).*mesh.edge_length(internal).*dt,[n 1],@sum,0);
    hv_change=hv_change+accumarray(neighbor(internal),flux_hv_neighbor(internal).*mesh.edge_length(internal).*dt,[n 1],@sum,0);
end
surface_volume=surface_volume+volume_change;
if any(surface_volume < -1e-10)
    error('HydroPol2D:NegativeSurfaceVolume','Full-momentum draining limiter failed.');
end
surface_volume=max(surface_volume,0);
hu=hu+hu_change./mesh.surface_area(:);
hv=hv+hv_change./mesh.surface_area(:);
hnew=surface_volume./mesh.surface_area(:);

% Exact pointwise implicit Manning friction.  It is safe for CPU now and
% vectorized for the later gpuArray execution path.
[hu,hv]=apply_manning_friction(hu,hv,hnew,roughness,dt,options.gravity,options.dry_tolerance_m);
[hu,hv]=limit_momentum(hu,hv,hnew,options.dry_tolerance_m,options.maximum_velocity_m_s);
internal_q=Q(internal); boundary_q=Q(~internal);
speed=sqrt(hu.^2+hv.^2)./max(hnew,options.dry_tolerance_m);
diagnostics=struct('max_depth_m',max(hnew),'max_velocity_m_s',max(speed,[],'omitnan'), ...
    'internal_flux_volume_m3',sum(abs(internal_q))*dt,'mass_change_m3',sum(volume_change));
boundary_diagnostics=struct('net_inflow_volume_m3',-sum(boundary_q)*dt, ...
    'max_discharge_m3_s',max(abs(boundary_q),[],'omitnan'));
end

function [nx,ny] = outward_normals(mesh, selector)
edge_id=find(selector); nx=mesh.edge_normal_x(edge_id); ny=mesh.edge_normal_y(edge_id);
norm_n=sqrt(nx.^2+ny.^2); nx=nx./max(norm_n,eps); ny=ny./max(norm_n,eps);
internal=mesh.edge_neighbor(edge_id)>0;
if any(internal)
    o=mesh.edge_owner(edge_id(internal)); d=mesh.edge_neighbor(edge_id(internal));
    dx=mesh.cell_x(d)-mesh.cell_x(o); dy=mesh.cell_y(d)-mesh.cell_y(o);
    flip=(nx(internal).*dx+ny(internal).*dy)<0;
    index=find(internal); nx(index(flip))=-nx(index(flip)); ny(index(flip))=-ny(index(flip));
end
end

function [types,values] = boundary_values(mesh,boundary,edge_id)
types=repmat("wall",numel(edge_id),1); values=zeros(numel(edge_id),1);
if isempty(boundary) || ~isfield(boundary,'edge_id') || isempty(boundary.edge_id), return; end
specified=double(boundary.edge_id(:));
assert(all(specified>=1 & specified<=mesh.n_edges & mesh.edge_neighbor(specified)==0));
kind=string(boundary.type(:)); if isscalar(kind), kind=repmat(kind,numel(specified),1); end
value=zeros(numel(specified),1); if isfield(boundary,'value'), value=double(boundary.value(:)); if isscalar(value), value=repmat(value,numel(specified),1); end, end
assert(numel(kind)==numel(specified) && numel(value)==numel(specified) && all(isfinite(value)));
[matched,position]=ismember(specified,edge_id); assert(all(matched));
valid=["wall" "inflow" "stage" "normal_flow" "critical_flow"];
assert(all(ismember(kind,valid)),'HydroPol2D:UnknownBoundaryType','Unknown full-momentum boundary type.');
types(position)=kind; values(position)=value;
end

function [u,v] = velocity_components(h,hu,hv,dry)
u=zeros(size(h)); v=zeros(size(h)); wet=h>dry;
u(wet)=hu(wet)./h(wet); v(wet)=hv(wet)./h(wet);
u(~isfinite(u))=0; v(~isfinite(v))=0;
end

function [Fh,Fqn,Fqt] = hll_normal_flux(hL,qnL,qtL,hR,qnR,qtR,g,dry)
[uL,vtL]=velocity_components(hL,qnL,qtL,dry); [uR,vtR]=velocity_components(hR,qnR,qtR,dry);
cL=sqrt(g.*max(hL,0)); cR=sqrt(g.*max(hR,0));
FLh=qnL; FLn=qnL.*uL+0.5*g*hL.^2; FLt=qtL.*uL;
FRh=qnR; FRn=qnR.*uR+0.5*g*hR.^2; FRt=qtR.*uR;
sL=min(uL-cL,uR-cR); sR=max(uL+cL,uR+cR);
Fh=zeros(size(hL)); Fqn=Fh; Fqt=Fh; left=sL>=0; right=sR<=0; middle=~(left|right) & sR>sL;
Fh(left)=FLh(left); Fqn(left)=FLn(left); Fqt(left)=FLt(left);
Fh(right)=FRh(right); Fqn(right)=FRn(right); Fqt(right)=FRt(right);
den=max(sR(middle)-sL(middle),eps);
Fh(middle)=(sR(middle).*FLh(middle)-sL(middle).*FRh(middle)+sL(middle).*sR(middle).*(hR(middle)-hL(middle)))./den;
Fqn(middle)=(sR(middle).*FLn(middle)-sL(middle).*FRn(middle)+sL(middle).*sR(middle).*(qnR(middle)-qnL(middle)))./den;
Fqt(middle)=(sR(middle).*FLt(middle)-sL(middle).*FRt(middle)+sL(middle).*sR(middle).*(qtR(middle)-qtL(middle)))./den;
bad=(hL<=dry & hR<=dry) | ~isfinite(Fh) | ~isfinite(Fqn) | ~isfinite(Fqt);
Fh(bad)=0; Fqn(bad)=0; Fqt(bad)=0;
end

function [hu,hv] = apply_manning_friction(hu,hv,h,roughness,dt,g,dry)
wet=h>dry; magnitude=sqrt(hu.^2+hv.^2); coefficient=zeros(size(h));
coefficient(wet)=g.*roughness(wet).^2./max(h(wet),dry).^(7/3);
a=dt.*coefficient; updated=zeros(size(h)); active=wet & magnitude>0;
updated(active)=2.*magnitude(active)./(1+sqrt(1+4.*a(active).*magnitude(active)));
scale=zeros(size(h)); scale(active)=updated(active)./magnitude(active); hu=hu.*scale; hv=hv.*scale;
end

function [hu,hv] = limit_momentum(hu,hv,h,dry,max_velocity)
wet=h>dry; hu(~wet)=0; hv(~wet)=0; speed=zeros(size(h)); speed(wet)=sqrt(hu(wet).^2+hv(wet).^2)./h(wet);
limited=wet & speed>max_velocity; factor=max_velocity./speed(limited); hu(limited)=hu(limited).*factor; hv(limited)=hv(limited).*factor;
hu(~isfinite(hu))=0; hv(~isfinite(hv))=0;
end
