negoTable = load("negoTable.mat"); 

C = cell2mat(negoTable.scores);
% C = [C2,cell2mat(negoTable.scores)];

n = size(C,1);
b = ones(n,1);
d0 = 1;
gamma = 0.9;
epsilon = 10;

[selection, trade, profit] = TACo(C, b, d0, gamma, epsilon);
selection
trade = trade(:,selection)
postProfit = profit(:,selection)
preProfit = -C(:,selection)