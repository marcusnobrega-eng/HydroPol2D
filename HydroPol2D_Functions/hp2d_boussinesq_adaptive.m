function [head,u_x,u_y,q_exf,mass_error,seepage_discharge_m3_s,diagnostics] = ...
        hp2d_boussinesq_adaptive(interval_s,initial_step_s,min_step_s, ...
        relative_tolerance,head_tolerance_m,dx,dy,head,bedrock,specific_yield, ...
        recharge_m_s,hydraulic_conductivity_m_s,surface,soil_depth,domain, ...
        dirichlet_mask,dirichlet_head,seepage_mask,seepage_head,perimeter,courant)
%HP2D_BOUSSINESQ_ADAPTIVE Error-controlled lateral groundwater integration.

% Local recharge, capillary rise, and river exchange are applied before this
% function. This function advances only lateral saturated flow over interval_s.

assert(interval_s>0 && min_step_s>0 && relative_tolerance>0 && head_tolerance_m>0, ...
    'HydroPol2D:InvalidGroundwaterTolerance', ...
    'Adaptive groundwater interval, floor, and tolerances must be positive.');

elapsed=0;
candidate=min(initial_step_s,interval_s);
if ~isfinite(candidate) || candidate<=0, candidate=interval_s; end
accepted=0; rejected=0;
exfiltration_depth_m=zeros(size(head));
seepage_volume_m3=0; mass_error=0;
u_x=zeros(size(head)); u_y=zeros(size(head));

while elapsed<interval_s-10*eps(max(interval_s,1))
    remaining=interval_s-elapsed;
    stable=local_stable_step(head,bedrock,specific_yield, ...
        hydraulic_conductivity_m_s,dx,dy,courant);
    step_s=min([candidate,stable,remaining]);
    floor_s=min(min_step_s,remaining);
    step_s=max(step_s,floor_s);

    [coarse,~,~,~,~,~,~]=take_step(head,step_s);
    [half_head,~,~,q_exf_1,error_1,seepage_1,~]=take_step(head,0.5*step_s);
    [fine,u_x_2,u_y_2,q_exf_2,error_2,seepage_2,~]=take_step(half_head,0.5*step_s);

    scale=head_tolerance_m+relative_tolerance.*max( ...
        abs(fine-bedrock),abs(coarse-bedrock));
    normalized=abs(fine-coarse)./max(scale,realmin);
    normalized(~domain)=0;
    local_error=max(normalized(:),[],'omitnan');
    if isempty(local_error) || ~isfinite(local_error)
        error('HydroPol2D:GroundwaterAdaptiveError', ...
            'Adaptive Boussinesq error estimate is not finite.');
    end

    if local_error>1
        if step_s<=floor_s*(1+1e-12)
            error('HydroPol2D:GroundwaterMinimumStep', ...
                ['Adaptive Boussinesq solver reached its internal minimum ', ...
                 'step before meeting the requested tolerance.']);
        end
        rejected=rejected+1;
        candidate=max(floor_s,step_s*max(0.2,0.9*local_error^(-1/3)));
        continue
    end

    head=fine; u_x=u_x_2; u_y=u_y_2;
    exfiltration_depth_m=exfiltration_depth_m+ ...
        0.5*step_s.*(max(q_exf_1,0)+max(q_exf_2,0));
    seepage_volume_m3=seepage_volume_m3+0.5*step_s*(seepage_1+seepage_2);
    mass_error=mass_error+error_1+error_2;
    elapsed=elapsed+step_s; accepted=accepted+1;
    if local_error==0
        growth=2;
    else
        growth=min(2,max(0.5,0.9*local_error^(-1/3)));
    end
    candidate=max(min_step_s,step_s*growth);
end

q_exf=exfiltration_depth_m./interval_s;
q_exf(~domain)=nan;
seepage_discharge_m3_s=seepage_volume_m3./interval_s;
diagnostics=struct('accepted_steps',accepted,'rejected_steps',rejected, ...
    'suggested_dt_s',candidate);

    function [h1,ux1,uy1,q1,e1,s1,qr1]=take_step(h0,dt_s)
        [h1,ux1,uy1,q1,qr1,e1,s1]=Boussinesq_2D_explicit( ...
            dt_s,dx,dy,h0,bedrock,specific_yield,recharge_m_s, ...
            hydraulic_conductivity_m_s,false(size(h0)),zeros(size(h0)), ...
            surface,surface,courant,soil_depth,domain,dirichlet_mask, ...
            dirichlet_head,seepage_mask,seepage_head,perimeter);
    end
end

function dt_s=local_stable_step(head,bedrock,specific_yield,conductivity,dx,dy,courant)
diffusivity=max(conductivity,0).*max(head-bedrock,1e-4)./max(specific_yield,1e-6);
maximum=max(diffusivity(:),[],'omitnan');
if isempty(maximum) || ~isfinite(maximum) || maximum<=0
    dt_s=inf;
else
    dt_s=courant*min(dx,dy)^2/maximum/4;
end
end
