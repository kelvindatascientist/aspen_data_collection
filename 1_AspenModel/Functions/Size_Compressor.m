function [W_Comp,N_Comp,Type_A,Type_B]=Size_Compressor(simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve data from cell called 'simulation'
% simulation{1} = Number_corresponding
% simulation{2} = Total_flow
load('equipments.mat')
Number_corresponding = simulation{1};
W_Total              = simulation{2};
  
% retrieve the name of the equipment
% Type_A = 'Pump'
Type_A=equipments{Number_corresponding,1,1}
% Type_B = '..'
Type_B=equipments{Number_corresponding,2,1}

% retrieve the min and max capacity
load('Turton.mat')
Unit_Row = (1:size(TableA1,1))*strcmp(TableA1(:,6),Type_B);
Cap_Min  = TableA1{Unit_Row,4};   
Cap_Max  = TableA1{Unit_Row,5};

% Calculate the number of pumps

N_Comp = ceil(W_Total/Cap_Max);
W_Comp = max(min(W_Total/N_Comp,Cap_Max),Cap_Min);

end