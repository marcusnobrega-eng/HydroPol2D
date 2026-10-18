function Summary=run_adaptive_coupling_validation()
%RUN_ADAPTIVE_COUPLING_VALIDATION Focused checks for the new coupling clock.
root=fullfile(fileparts(mfilename('fullpath')),'..','..');
addpath(fullfile(root,'HydroPol2D_Functions'));

storage_area_gw=20; storage_area_river=100;
gw_volume=100; river_volume=100;
gw_head=2; river_stage=1;
[single,~,single_residual]=hp2d_river_groundwater_exchange_step( ...
    1e-6,0.5,100,gw_head,river_stage,gw_volume,river_volume, ...
    storage_area_gw,storage_area_river,600,true);

split_total=0;
for i=1:10
    [moved,~,~]=hp2d_river_groundwater_exchange_step( ...
        1e-6,0.5,100,gw_head,river_stage,gw_volume,river_volume, ...
        storage_area_gw,storage_area_river,60,true);
    split_total=split_total+moved;
    gw_volume=gw_volume-moved;
    river_volume=river_volume+moved;
    gw_head=gw_head-moved/storage_area_gw;
    river_stage=river_stage+moved/storage_area_river;
end
exchange_difference=abs(single-split_total);

bedrock=zeros(3,5);
head=ones(3,5);
head(:,3)=5;
surface=20.*ones(3,5);
specific_yield=0.2.*ones(3,5);
conductivity=1e-2.*ones(3,5);
domain=true(3,5);
storage_before=sum(specific_yield.*(head-bedrock).*100.*100,'all');
[head_after,~,~,~,mass_residual,~,adaptive]=hp2d_boussinesq_adaptive( ...
    600,600,0.01,1e-8,1e-9,100,100,head,bedrock,specific_yield, ...
    zeros(3,5),conductivity,surface,surface-bedrock,domain,[],[],[],[], ...
    false(3,5),0.25);
storage_after=sum(specific_yield.*(head_after-bedrock).*100.*100,'all');
storage_difference=abs(storage_after-storage_before);

assert(exchange_difference<1e-10, ...
    'HydroPol2D:ExchangeStepDependence','Exact river exchange depends on step subdivision.');
assert(abs(single_residual)<1e-12, ...
    'HydroPol2D:ExchangeMass','River exchange does not conserve mass.');
assert(adaptive.accepted_steps>0 && adaptive.rejected_steps>0, ...
    'HydroPol2D:AdaptiveClock','Adaptive Boussinesq controller did not exercise rejection.');
assert(storage_difference<1e-8 && abs(mass_residual)<1e-8, ...
    'HydroPol2D:BoussinesqMass','Closed-domain lateral flow does not conserve mass.');

Summary=table(exchange_difference,storage_difference,mass_residual, ...
    adaptive.accepted_steps,adaptive.rejected_steps,adaptive.suggested_dt_s, ...
    'VariableNames',{'exchange_step_difference_m3','storage_difference_m3', ...
    'mass_residual_m3','accepted_steps','rejected_steps','suggested_dt_s'});
disp(Summary);
end
