function[fis_controller] = Fuzzy_control(d_lim,A_inf,A_sup)

%%Creation of fuzzy interferance

fis_controller(1:length(d_lim)) = mamfis;

for i = 1:length(d_lim)

dlimref = d_lim(i); % depth reference per reservior
Ainf = A_inf(i);
Asup = A_sup(i);

varref = 0.2*dlimref;%variation limit
percentageStep = 0.20;% Define the percentage step (20%)
stepSize = percentageStep * (Asup - Ainf);% Calculate the step size
parameterValues = Ainf:stepSize:Asup;% Generate the values

reservfis = mamfis('Name', 'ReservoirFuzzyController');

% Definição das variáveis de entrada
reservfis = addInput(reservfis, [0, dlimref], 'Name', 'WaterLevel');
reservfis = addInput(reservfis, [-varref, varref], 'Name', 'VarWaterLevel');
%fis = addInput(fis, [0, max(qbase)], 'Name', 'FlowDifference');

% Definição da variável de saída
reservfis = addOutput(reservfis, [Ainf, Asup], 'Name', 'ValveOpening');

% Definir funções de pertinência para a lâmina de água
reservfis = addMF(reservfis, 'WaterLevel', 'zmf', [0.2*dlimref 0.5*dlimref], 'Name', 'Low');
reservfis = addMF(reservfis, 'WaterLevel', 'gbellmf', [0.2*dlimref 0.8*dlimref 0.5*dlimref], 'Name', 'Medium');
reservfis = addMF(reservfis, 'WaterLevel', 'smf', [0.5*dlimref 0.8*dlimref], 'Name', 'High');


% Definir funções de pertinência para a diferença de fluxo
reservfis = addMF(reservfis, 'VarWaterLevel', 'linzmf', [-0.75*varref -0.5*varref], 'Name', 'NH');
reservfis = addMF(reservfis, 'VarWaterLevel', 'trapmf', [-0.75*varref -0.5*varref -0.25*varref -0.1*varref], 'Name', 'NL');
reservfis = addMF(reservfis, 'VarWaterLevel', 'trapmf', [-0.25*varref -0.1*varref 0.1*varref 0.25*varref], 'Name', 'Z');
reservfis = addMF(reservfis, 'VarWaterLevel', 'trapmf', [0.1*varref 0.25*varref 0.5*varref 0.75*varref], 'Name', 'PL');
reservfis = addMF(reservfis, 'VarWaterLevel', 'linsmf', [0.5*varref 0.75*varref], 'Name', 'PH');



% Definir funções de pertinência para a saída (abertura da válvula)
reservfis = addMF(reservfis, 'ValveOpening', 'linzmf', [parameterValues(2)-0.5*stepSize parameterValues(2)], 'Name', 'Op_20');
reservfis = addMF(reservfis, 'ValveOpening', 'trapmf', [parameterValues(2) parameterValues(2)+0.5*stepSize parameterValues(3)-0.5*stepSize parameterValues(4)-0.5*stepSize], 'Name', 'Op_40');
reservfis = addMF(reservfis, 'ValveOpening', 'trapmf', [parameterValues(3)-0.5*stepSize parameterValues(3)+0.5*stepSize parameterValues(4)-0.5*stepSize parameterValues(5)-0.5*stepSize], 'Name', 'Op_60');
reservfis = addMF(reservfis, 'ValveOpening', 'trapmf', [parameterValues(4)-0.5*stepSize parameterValues(4)+0.5*stepSize parameterValues(5)-0.5*stepSize parameterValues(5)], 'Name', 'Op_80');
reservfis = addMF(reservfis, 'ValveOpening', 'linsmf', [parameterValues(5) parameterValues(6)-0.5*stepSize], 'Name', 'Op_100');

% Definir as regras fuzzy
ruleList = [
    "If WaterLevel is Low and VarWaterLevel is NH then ValveOpening is Op_60"
    "If WaterLevel is Low and VarWaterLevel is NL then ValveOpening is Op_100"
    "If WaterLevel is Low and VarWaterLevel is Z then ValveOpening is Op_80"
    "If WaterLevel is Low and VarWaterLevel is PL then ValveOpening is Op_20"
    "If WaterLevel is Low and VarWaterLevel is PH then ValveOpening is Op_40"
    "If WaterLevel is Medium and VarWaterLevel is NH then ValveOpening is Op_40"
    "If WaterLevel is Medium and VarWaterLevel is NL then ValveOpening is Op_80"
    "If WaterLevel is Medium and VarWaterLevel is Z then ValveOpening is Op_40"
    "If WaterLevel is Medium and VarWaterLevel is PL then ValveOpening is Op_40"
    "If WaterLevel is Medium and VarWaterLevel is PH then ValveOpening is Op_60"
    "If WaterLevel is High and VarWaterLevel is NH then ValveOpening is Op_20"
    "If WaterLevel is High and VarWaterLevel is NL then ValveOpening is Op_60"
    "If WaterLevel is High and VarWaterLevel is Z then ValveOpening is Op_20"
    "If WaterLevel is High and VarWaterLevel is PL then ValveOpening is Op_80"
    "If WaterLevel is High and VarWaterLevel is PH then ValveOpening is Op_100"
    
];

for j = 1:length(ruleList)
    reservfis = addRule(reservfis, ruleList(j));
end

%%%%test
%{
plotmf(fis, 'input', 1);  % Mostra as funções de pertinência da WaterLevel
plotmf(fis, 'input', 2);  % Mostra as funções de pertinência da FlowDifference
plotmf(fis, 'output', 1); % Mostra as funções de pertinência da ValveOpening

test_k = evalfis(reservfis, [0.8*max(dlimref), 1]); % Teste com valores médios
disp(['Abertura da válvula recomendada: ', num2str(test_k)]);

%}

fis_controller(i) = reservfis;

end
end
