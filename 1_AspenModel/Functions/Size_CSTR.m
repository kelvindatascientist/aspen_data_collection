function [V_Req,V_Reactor,D_Reactor,L_Reactor,N_Reactor,Type_A,Type_B]=Size_CSTR(simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve data from cell called 'simulation'
% simulation{1} = Number_corresponding
% simulation{2} = Total_flow
load('equipments.mat')
Number_corresponding = simulation{1};
Total_flow           = simulation{2};

% retrieve the name of the equipment
% Type_A = 'Reactor'
Type_A=equipments{Number_corresponding,1,1}
% Type_B = 'Reactor Jacketed Agitated'
Type_B=equipments{Number_corresponding,2,1}

% retrieve the min and max capacity
load('Turton.mat')
Unit_Row = (1:size(TableA1,1))*strcmp(TableA1(:,6),Type_B);
Cap_Min  = TableA1{Unit_Row,4};   
Cap_Max  = TableA1{Unit_Row,5};

% Define the ratio
LDr = 5;  % L to D ratio

% Calculate the volume required based on simulation
% from the article of Cui (2017), there was a 15h time to convert >80% of
% methanol to acetic acid, therefore in this study 1/15 h-1 was used as
% residence time.
tau = 1/15;
V_Req = Total_flow*tau

% Calculate the number of reactors	
N_Reactor = ceil(V_Req/Cap_Max);
V_Prelim  = V_Req/N_Reactor;                % Prelim reactor volume
V_Reactor = max(V_Prelim,Cap_Min );         % Reactor volume
D_Reactor = (4*V_Reactor/(pi*LDr))^(1/3);	% Reactor diameter
L_Reactor = LDr*D_Reactor;                  % Reactor length  

end