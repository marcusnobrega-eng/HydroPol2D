%% Process Potential Evapotranspiration
% Matheus Schroeder dos Santos
% 15/10/2022
% Objective:  Developed a 2-D array with the potential evapotranspiration for all cells in a given watershed
% Bibliography <Marco Antônio Fonseca Conceição. (2006). Roteiro de cálculo da evapotranspiração de referência pelo método de Penman-Monteith-FAO. Circular Técnica - EMBRAPA,(ISSN 1808-6810), 1–8>.

function [ETP] = Evapotranspiration(DEM,Temp,Temp_max,Temp_min,Day_year,lat,U_2_input,U_R_input,Krs,alfa_albedo_input,G_input)
%% Input data
% Change data
% Temp = 18; % Temperatura média do ar (C)
% lat = -23.62;% latitude (graus)
% Day_year = 10; % Dia do ano (1 a 365 ou 366)
% Temp_max = 23; % Temperatura máxima do dia (C)
% Temp_min = 16; % Temperatura mínima do dia (C)

% Constant data
% U_2_input = 1.6; % Velocidade do vento a 2m de altura (m/s) - (e.g fixo em 2.0 como default)
% U_R_input = 81.6; % Umidade relativa do ar (%)
% alfa_albedo_input = 0.23; % Cultura de referência (grama)
%a = 0.25;% Coeficiente local para determinação da Radiação solar incidente (e.g fixo em 0.25 como default)
%b = 0.50; % Coeficiente local para determinação da Radiação solar incidente (e.g fixo em 0.50 como default)
% Krs = 0.16; % Coeficiente local para determinação da Radiação solar incidente (e.g fixo em 0.16 - regiões continentais, e 0.19 - regiões costeiras, como default).
teta_S_B_input = 4.903*(10^-9); % Constante de Stefan-Boltzmann (MJm^-2/dia)



% G_input = 0.6; % Fluxo total diário de calor no solo (MJm^-2/dia)
n_rows = size(DEM,1);
n_cols = size(DEM,2);
latitude = lat;
G = G_input;
T = Temp;
U_2 = U_2_input;
U_R = U_R_input;
Tmax = Temp_max;
Tmin = Temp_min;
%% Processing data
%DEM
neg_DEM = DEM < -200;
high_DEM = DEM > 99999;
inf_nan_MAPS = isinf(DEM) + isnan(DEM) + neg_DEM + high_DEM; % Logical array
idx = inf_nan_MAPS > 0;

% Setting Up Mask in All Data
G(idx) = nan;
T(idx) = nan;
U_2(idx) = nan;
U_R(idx) = nan;
Tmax(idx) = nan;
Tmin(idx) = nan;

% Dia do ano
J = Day_year*ones(n_rows, n_cols);
J(idx) = nan;

% alfa_albedo
alfa_albedo = alfa_albedo_input*ones(n_rows, n_cols);
alfa_albedo(idx) = nan;

% teta_S_B
teta_S_B = teta_S_B_input*ones(n_rows, n_cols);
teta_S_B(idx) = nan;

%% Equations
% clearvars -except G teta_S_B alfa_albedo U_2 U_R J latitude Rs Tmin Tmax T n_cols n_rows DEM Krs a b

phi = latitude*pi/180; % latitude (radianos)
dec_sol = 0.409*sin((2*pi/365)*J-1.39); %Declinação solar (radianos)
if (1-((tan(phi)).^2).*((tan(dec_sol)).^2)) <= 0
    X = 0.00001;
else
    X = (1-((tan(phi)).^2).*((tan(dec_sol)).^2));
end
ws = (pi/2)-atan(((-tan(phi)).*tan(dec_sol))./(X.^0.5)); % Ângulo horário ao nascer do sol (radianos)
%N = (24/pi).*ws;
dr = 1+0.033*cos(((2*pi)/365)*J); % Distância relativa entre a Terra e o Sol (radianos)
Ra = (118.08/pi)*dr.*(ws.*sin(phi).*sin(dec_sol)+cos(phi).*cos(dec_sol).*sin(ws)); % Radiação solar no topo da atmosfera (MJm^-2/dia)
%Rs_input = (a+(b*(n/N)))*Ra; % Radiação solar incidente (MJm^-2/dia) calculado em função do número de horas com radiação solar
Rs_input = Krs*Ra.*sqrt((Tmax-Tmin)); % Radiação solar incidente (MJm^-2/dia) calculado em função da temperatura
Rso = (0.75+2*(10^-5)*DEM).*Ra; % Radiação solar incidente na ausência de nuvens (MJm^-2/dia)
Rns = (1-alfa_albedo).*Rs_input; %Saldo de radiação de ondas curtas (MJm^-2/dia)
e_s = 0.6108*exp((17.27*T)./(T+237.3)); % Pressão de saturação do vapor (kPa)
%e_a = e_s.*U_R/100; % Pressão atual de vapor (kPa)
e_a = 0.61*exp((17.27*Tmin)./(Tmin+237.3)); % Pressão atual de vapor (kPa) quando não se tem U_R
Rnl = teta_S_B.*((((Tmax+273.16).^4)+((Tmin+273.16).^4))./2).*(0.34-0.14.*sqrt(e_a)).*(1.35*(Rs_input./Rso)-0.35); % Saldo da radiação de ondas longas (MJm^-2/dia)
Rn = Rns-Rnl; % Saldo de radiação diária (MJm^-2/dia)
delta = (4098*(0.6108*exp((17.27*T)./(T+237.3))))./((T+237.3).^2); % Declividade da curva de pressão de vapor em relação a temperatura (kPa/C)
Patm = 101.3*((((293-0.0065*DEM)/293)).^5.26); % Pressão atmosférica local (kPa)
gama = 0.665*(10^-3)*Patm; % Coeficiente psicométrico (kPa/C)
ETP = (0.408.*delta.*(Rn-G)+(gama.*900.*U_2.*(e_s-e_a))./(T+273))./(delta+gama.*(1+0.34.*U_2));

% Mask
ETP(idx) = nan;

end
