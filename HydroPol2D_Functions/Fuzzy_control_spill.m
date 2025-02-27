function[fis_controller] = Fuzzy_control_spill(d_lim,A_inf,A_sup,tgate,reserv_id)

%%Creation of fuzzy interferance

fis_controller(1:length(d_lim)) = mamfis;
for i = 1:length(d_lim)

dlimref = d_lim(i); % depth reference per reservior
Ainf = A_inf(i);
Asup = A_sup(i);

%{
dlimref = 3.02; % depth reference per reservior
kinf = 0;
ksup = 4.52;
%}

if tgate(i) == 0
    varref = 0.15*dlimref;%variation limit (15%)
    percentageStep = 0.25;% Define the percentage step (25%)
    stepSize = percentageStep * (Asup - Ainf);% Calculate the step size
    parameterValues = Ainf:stepSize:Asup;% Generate the values

    reservfis = mamfis('Name', 'ReservoirFuzzyController');

    % Definição das variáveis de entrada
    reservfis = addInput(reservfis, [0, dlimref], 'Name', 'WaterLevel');
    reservfis = addInput(reservfis, [-varref, varref], 'Name', 'VarWaterLevel');

    % Definição da variável de saída
    reservfis = addOutput(reservfis, [Ainf, Asup], 'Name', 'ValveOpening');

    % Definir funções de pertinência para a lâmina de água
    reservfis = addMF(reservfis, 'WaterLevel', 'linzmf', [0.2*dlimref 0.4*dlimref], 'Name', 'Low');
    reservfis = addMF(reservfis, 'WaterLevel', 'trapmf', [0.2*dlimref 0.4*dlimref 0.6*dlimref 0.8*dlimref], 'Name', 'Medium');
    reservfis = addMF(reservfis, 'WaterLevel', 'linsmf', [0.6*dlimref 0.8*dlimref], 'Name', 'High');

    % Definir funções de pertinência para a diferença de fluxo
    reservfis = addMF(reservfis, 'VarWaterLevel', 'linzmf', [-0.75*varref -0.5*varref], 'Name', 'NH');
    reservfis = addMF(reservfis, 'VarWaterLevel', 'trapmf', [-0.75*varref -0.5*varref -0.25*varref -0.1*varref], 'Name', 'NL');
    reservfis = addMF(reservfis, 'VarWaterLevel', 'trapmf', [-0.25*varref -0.1*varref 0.1*varref 0.25*varref], 'Name', 'Z');
    reservfis = addMF(reservfis, 'VarWaterLevel', 'trapmf', [0.1*varref 0.25*varref 0.5*varref 0.75*varref], 'Name', 'PL');
    reservfis = addMF(reservfis, 'VarWaterLevel', 'linsmf', [0.5*varref 0.75*varref], 'Name', 'PH');

    % Definir funções de pertinência para a saída (abertura da comporta)
    reservfis = addMF(reservfis, 'ValveOpening', 'trimf', [parameterValues(1), parameterValues(1), parameterValues(2)], 'Name', 'Op_0');
    reservfis = addMF(reservfis, 'ValveOpening', 'trimf', [parameterValues(1), parameterValues(2), parameterValues(3)], 'Name', 'Op_25');
    reservfis = addMF(reservfis, 'ValveOpening', 'trimf', [parameterValues(2), parameterValues(3), parameterValues(4)], 'Name', 'Op_50');
    reservfis = addMF(reservfis, 'ValveOpening', 'trimf', [parameterValues(3), parameterValues(4), parameterValues(5)], 'Name', 'Op_75');
    reservfis = addMF(reservfis, 'ValveOpening', 'trimf', [parameterValues(4), parameterValues(5), parameterValues(5)], 'Name', 'Op_100');

    % Definir as regras fuzzy para o culvert
    
    if reserv_id(i) == 3
        %AR3 rules culv
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is Z then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is PL then ValveOpening is Op_50"
            "If WaterLevel is Low and VarWaterLevel is PH then ValveOpening is Op_50"
            "If WaterLevel is Medium and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is Z then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is PL then ValveOpening is Op_75"
            "If WaterLevel is Medium and VarWaterLevel is PH then ValveOpening is Op_75"
            "If WaterLevel is High and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is NL then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is Z then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is PL then ValveOpening is Op_100"
            "If WaterLevel is High and VarWaterLevel is PH then ValveOpening is Op_100"
        ]; 
    elseif reserv_id(i) == 4  
        %Caguacu rules culv
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is Z then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is PL then ValveOpening is Op_50"
            "If WaterLevel is Low and VarWaterLevel is PH then ValveOpening is Op_50"
            "If WaterLevel is Medium and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is Z then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is PL then ValveOpening is Op_75"
            "If WaterLevel is Medium and VarWaterLevel is PH then ValveOpening is Op_75"
            "If WaterLevel is High and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is NL then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is Z then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is PL then ValveOpening is Op_100"
            "If WaterLevel is High and VarWaterLevel is PH then ValveOpening is Op_100"
        ];
    elseif reserv_id(i) == 1    
        % AR1 rules culvert
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is Z then ValveOpening is Op_100"
            "If WaterLevel is Low and VarWaterLevel is PL then ValveOpening is Op_50"
            "If WaterLevel is Low and VarWaterLevel is PH then ValveOpening is Op_50"
            "If WaterLevel is Medium and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is Z then ValveOpening is Op_100"
            "If WaterLevel is Medium and VarWaterLevel is PL then ValveOpening is Op_50"
            "If WaterLevel is Medium and VarWaterLevel is PH then ValveOpening is Op_75"
            "If WaterLevel is High and VarWaterLevel is NH then ValveOpening is Op_75"
            "If WaterLevel is High and VarWaterLevel is NL then ValveOpening is Op_75"
            "If WaterLevel is High and VarWaterLevel is Z then ValveOpening is Op_100"
            "If WaterLevel is High and VarWaterLevel is PL then ValveOpening is Op_100"
            "If WaterLevel is High and VarWaterLevel is PH then ValveOpening is Op_100"
        ];
    elseif reserv_id(i) == 5    
        % limo culv
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is Z then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is PL then ValveOpening is Op_50"
            "If WaterLevel is Low and VarWaterLevel is PH then ValveOpening is Op_50"
            "If WaterLevel is Medium and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is Z then ValveOpening is Op_50"
            "If WaterLevel is Medium and VarWaterLevel is PL then ValveOpening is Op_75"
            "If WaterLevel is Medium and VarWaterLevel is PH then ValveOpening is Op_75"
            "If WaterLevel is High and VarWaterLevel is NH then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is NL then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is Z then ValveOpening is Op_75"
            "If WaterLevel is High and VarWaterLevel is PL then ValveOpening is Op_100"
            "If WaterLevel is High and VarWaterLevel is PH then ValveOpening is Op_100"
        ];
    elseif reserv_id(i) == 2    
        % AR2 rules culvert
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is Z then ValveOpening is Op_25"
            "If WaterLevel is Low and VarWaterLevel is PL then ValveOpening is Op_50"
            "If WaterLevel is Low and VarWaterLevel is PH then ValveOpening is Op_50"
            "If WaterLevel is Medium and VarWaterLevel is NH then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is NL then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is Z then ValveOpening is Op_25"
            "If WaterLevel is Medium and VarWaterLevel is PL then ValveOpening is Op_75"
            "If WaterLevel is Medium and VarWaterLevel is PH then ValveOpening is Op_75"
            "If WaterLevel is High and VarWaterLevel is NH then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is NL then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is Z then ValveOpening is Op_50"
            "If WaterLevel is High and VarWaterLevel is PL then ValveOpening is Op_100"
            "If WaterLevel is High and VarWaterLevel is PH then ValveOpening is Op_100"
        ];
    end

    for j = 1:length(ruleList)
        reservfis = addRule(reservfis, ruleList(j));
    end

    %%%%test
    %{
    plotmf(reservfis, 'input', 1);  % Mostra as funções de pertinência da WaterLevel
    plotmf(reservfis, 'input', 2);  % Mostra as funções de pertinência da FlowDifference
    plotmf(reservfis, 'output', 1); % Mostra as funções de pertinência da ValveOpening

    % Test multiple cases
    test_cases = [  % [WaterLevel, VarWaterLevel]
        0.8*dlimref, -0.5*varref;
        0.4*dlimref,  0.2*varref;
        0.9*dlimref,  0.7*varref;
        0.1*dlimref, -0.6*varref;
    ];

    disp('Test Cases:');
    for i = 1:size(test_cases,1)
        test_k = evalfis(reservfis, test_cases(i, :));
        disp(['Test ', num2str(i), ': Water Level = ', num2str(test_cases(i,1)), ...
          ', VarWaterLevel = ', num2str(test_cases(i,2)), ...
          ' → Valve Opening: ', num2str(test_k)]);
    end
    %}
elseif tgate(i) == 1 %spill gate
    varref = 0.15*dlimref;%variation limit (15%)
    percentageStep = 0.25;% Define the percentage step (25%)
    stepSize = percentageStep * (Asup - Ainf);% Calculate the step size
    parameterValues = Ainf:stepSize:Asup;% Generate the values

    reservfis = mamfis('Name', 'ReservoirFuzzyController');

    % Definição das variáveis de entrada
    reservfis = addInput(reservfis, [0, dlimref], 'Name', 'WaterLevel');
    reservfis = addInput(reservfis, [-varref, varref], 'Name', 'VarWaterLevel');
    %fis = addInput(fis, [0, max(qbase)], 'Name', 'FlowDifference');

    % Definição da variável de saída
    reservfis = addOutput(reservfis, [Ainf, Asup], 'Name', 'GateOpening');

    % Definir funções de pertinência para a lâmina de água
    reservfis = addMF(reservfis, 'WaterLevel', 'linzmf', [0.25/2*dlimref 0.25*dlimref], 'Name', 'Low');
    reservfis = addMF(reservfis, 'WaterLevel', 'trimf', [0.25*dlimref 0.375*dlimref 0.5*dlimref], 'Name', 'MediumL');
    reservfis = addMF(reservfis, 'WaterLevel', 'trimf', [0.5*dlimref 0.625*dlimref 0.75*dlimref], 'Name', 'MediumH');
    reservfis = addMF(reservfis, 'WaterLevel', 'linsmf', [0.75*dlimref 0.875*dlimref], 'Name', 'High');

    % Definir funções de pertinência para a diferença de fluxo
    reservfis = addMF(reservfis, 'VarWaterLevel', 'zmf', [0 0.5*varref], 'Name', 'N');
    reservfis = addMF(reservfis, 'VarWaterLevel', 'smf', [0.1*varref 0.5*varref], 'Name', 'P');

    % Definir funções de pertinência para a saída (abertura da válvula)
    reservfis = addMF(reservfis, 'GateOpening', 'linzmf', [(parameterValues(2)-parameterValues(1))*0.5 parameterValues(2)], 'Name', 'Op_0');
    reservfis = addMF(reservfis, 'GateOpening', 'trapmf', [(parameterValues(2)-parameterValues(1))*0.5, parameterValues(2)-0.1*stepSize, parameterValues(2)+0.1*stepSize, parameterValues(2)+(parameterValues(3)-parameterValues(2))*0.5], 'Name', 'Op_25');
    reservfis = addMF(reservfis, 'GateOpening', 'trapmf', [parameterValues(2)+(parameterValues(3)-parameterValues(2))*0.5, parameterValues(3)-0.1*stepSize, parameterValues(3)+0.1*stepSize , parameterValues(3)+(parameterValues(3)-parameterValues(2))*0.5], 'Name', 'Op_50');
    reservfis = addMF(reservfis, 'GateOpening', 'trimf', [parameterValues(3)+(parameterValues(3)-parameterValues(2))*0.5, parameterValues(4), parameterValues(4)+(parameterValues(4)-parameterValues(3))*0.5], 'Name', 'Op_75');
    reservfis = addMF(reservfis, 'GateOpening', 'linsmf', [parameterValues(4) parameterValues(4)+(parameterValues(4)-parameterValues(3))*0.5 ], 'Name', 'Op_100');

    % Definir as regras fuzzy para o vertedouro
    
    if reserv_id(i) == 3
        %AR3 rules spill
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is Low and VarWaterLevel is P then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is MediumH and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is MediumH and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is P then GateOpening is Op_50"
        ];
    elseif reserv_id(i) == 4
        %Caguacu rules spill
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is Low and VarWaterLevel is P then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is MediumL and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is MediumH and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is MediumH and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is P then GateOpening is Op_50"
        ];
    elseif reserv_id(i) == 1
        %AR1 rules spill
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is Low and VarWaterLevel is P then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is P then GateOpening is Op_0"
            "If WaterLevel is MediumH and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is MediumH and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is P then GateOpening is Op_50"
        ];
    elseif reserv_id(i) == 5    
        %Limo rules spill
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is Low and VarWaterLevel is P then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is MediumH and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is MediumH and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is P then GateOpening is Op_50"
        ];
    elseif reserv_id(i) == 2
        %AR2 rules spill
        ruleList = [
            "If WaterLevel is Low and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is Low and VarWaterLevel is P then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is MediumL and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is MediumH and VarWaterLevel is N then GateOpening is Op_0"
            "If WaterLevel is MediumH and VarWaterLevel is P then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is N then GateOpening is Op_25"
            "If WaterLevel is High and VarWaterLevel is P then GateOpening is Op_50"
        ];
    end   

    for j = 1:length(ruleList)
        reservfis = addRule(reservfis, ruleList(j));
    end

    %%%%test
    %{
    plotmf(reservfis, 'input', 1);  % Mostra as funções de pertinência da WaterLevel
    plotmf(reservfis, 'input', 2);  % Mostra as funções de pertinência da FlowDifference
    plotmf(reservfis, 'output', 1); % Mostra as funções de pertinência da ValveOpening

    % Test multiple cases
    test_cases = [  % [WaterLevel, VarWaterLevel]
        0.8*dlimref, -0.5*varref;
        0.4*dlimref,  0.2*varref;
        0.9*dlimref,  0.7*varref;
        0.1*dlimref, -0.6*varref;
    ];

    disp('Test Cases:');
    for i = 1:size(test_cases,1)
        test_k = evalfis(reservfis, test_cases(i, :));
        disp(['Test ', num2str(i), ': Water Level = ', num2str(test_cases(i,1)), ...
          ', VarWaterLevel = ', num2str(test_cases(i,2)), ...
          ' → Valve Opening: ', num2str(test_k)]);
    end
    %}
end

fis_controller(i) = reservfis;

end
end
