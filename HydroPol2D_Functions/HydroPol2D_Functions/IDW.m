%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INVERSE DISTANCE WEIGHT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Vint]=IDW(xc,yc,vc,x,y,e,r1,r2)
%%%%%%%%%%%%%%%%%
%%% INPUTS
%xc = stations x coordinates (columns) [vector]
%yc = stations y coordinates (rows) [vector]
%vc = variable values on the point [xc yc]
%x = interpolation points  x coordinates [vector]
%y = interpolation points y coordinates [vector]
%e = distance weight
%r1 --- 'fr' = fixed radius ;  'ng' = neighbours
%r2 --- radius lenght if r1 == 'fr' / number of neighbours if  r1 =='ng'
%%% OUTPUTS
%Vint --- Matrix [length(y),length(x)] with interpolated  variable values
%%% EXAMPLES
%%% --> V_spa=IDW(x1,y1,v1,x,y,-2,'ng',length(x1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simone Fatichi -- simonef@dicea.unifi.it
%   Copyright 2009
%   $Date: 2009/06/19 $ 
%   $Updated 2012/02/24 $ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vint=zeros(length(y),length(x));
xc=reshape(xc,1,length(xc));
yc=reshape(yc,1,length(yc));
vc=reshape(vc,1,length(vc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  strcmp(r1,'fr')
    if  (r2<=0)
        disp('Error: Radius must be positive')
        return
    end
    for i=1:length(x)
        for j=1:length(y)
            D=[]; V=[]; wV =[]; vcc=[];
            D= sqrt((x(i)-xc).^2 +(y(j)-yc).^2);
            if min(D)==0
                disp('Error: One or more stations have the coordinates of an interpolation point')
                return
            end
            vcc=vc(D<r2); D=D(D<r2); 
            V = vcc.*(D.^e);
            wV = D.^e;
            if isempty(D)
                V=NaN;
            else
                V=sum(V)/sum(wV);
            end
            Vint(j,i)=V;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if (r2 > length(vc)) || (r2<1)
        disp('Error: Number of neighbours not congruent with data')
        return
    end
    for i=1:length(x)
        for j=1:length(y)
            D=[]; V=[]; wV =[];vcc=[];
            D= sqrt((x(i)-xc).^2 +(y(j)-yc).^2);
            if min(D)==0
                disp('Error: One or more stations have the coordinates of an interpolation point')
                return
            end
            [D,I]=sort(D);
            vcc=vc(I);
            V = vcc(1:r2).*(D(1:r2).^e);
            wV = D(1:r2).^e;
            V=sum(V)/sum(wV);
            Vint(j,i)=V;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
