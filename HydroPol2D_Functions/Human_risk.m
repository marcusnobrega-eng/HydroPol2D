function [map] = human_risk(flag,v,h,ro,g,mu,cd,slope,m,y,w,d)   
    if flag == 3
        % Person's Geometry for h<y
        sub_temp_1 = h./cosd(slope) < y./2;
        V_s1 = sub_temp_1.*(((2*h.*d^2).*pi())./(cosd(slope).*4));
        A_s1 = sub_temp_1.*((2*h.*d^2)./cosd(slope));
        X_gs1 = sub_temp_1.*(d/2);
        Y_gs1 = sub_temp_1.*(h./(2*cosd(slope)));
         
        % Person's Geometry for h>y
        sub_temp_2 = h./cosd(slope) > y/2;
        V_s2 = sub_temp_2.*(((y*pi()*d^2)/4) + ((h./cosd(slope) - y/2).*pi*((w^2)/4)));
        A_s2 = sub_temp_2.*((y*d) + ((h./cosd(slope))-(y/2))*w);
        X_gs2 = sub_temp_2.*(((2*((pi()*d^2)/4)*(y/2)*(d/2)) + ((pi()*w^2)/4)*((h./cosd(slope))-(y/2))*(w/2))./(2*((pi()*d^2)/4)*(y/2) + ((pi()*w^2)/4)*((h./cosd(slope))-(y/2))));
        Y_gs2 = sub_temp_2.*((2*((pi()*d^2)/4)*((y^2)/8) + ((pi()*w^2)/4)*((h./cosd(slope))-(y/2)).*((0.5*(h./cosd(slope)-(y/2)))+(y/2)))./(2*((pi()*d^2)/4)*(y/2) + ((pi()*w^2)/4)*((h./cosd(slope))-(y/2))));
        
        % Drag forces by the flow
        L1 = 0.5*ro*cd*(sind(90-slope).^2).*cosd(90-slope).*(v.^2).*A_s1;
        L2 = 0.5*ro*cd*(sind(90-slope).^2).*cosd(90-slope).*(v.^2).*A_s2;
        L = L1+L2;
        D1 = 0.5*ro*cd*(sind(90-slope).^3).*(v.^2).*A_s1;
        D2 = 0.5*ro*cd*(sind(90-slope).^3).*(v.^2).*A_s2;
        Drag = D1+D2;

        % Person weight forces acting
        Wp = m*g*sind(slope);
        Wn = m*g*cosd(slope);

        % Buoyancy effect
        Bn1 = ro*g*V_s1.*cosd(slope);
        Bn2 = ro*g*V_s2.*cosd(slope);
        Bn = Bn1+Bn2;
        X_gs = X_gs1 + X_gs2;
        Y_gs = Y_gs1 + Y_gs2;

        % Friccion available
        T = mu*(Wn-Bn-L);

        % Momemtum acting
        M = Drag.*(h./2) + (Wp.*(7/12).*y.*cosd(slope)) + (Bn.*((X_gs./cosd(slope)) + (Y_gs.*sind(slope)))) + L.*((X_gs./cosd(slope)) + ((h./2).*tan(slope)));

        % Slipping stability
        slide = double(Drag + Wp > T)*1;

        % Toppling stability
        topple = double(M > Wn.*((((5/6)*d)./cosd(slope)) + (((7/12)*y)./sind(slope))))*2;

        % Drawing risk
        drawing = double(h > ((13/16)*y))*3;

        %final
        map = max(max(slide, topple), drawing);
    end
end