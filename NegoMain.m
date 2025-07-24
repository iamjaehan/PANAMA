taxParam = 10;

negoTable = load("negoTable.mat"); 
C = cell2mat(negoTable.scores);

n = size(C,1);
b = ones(n,1);
d0 = 1;
gamma = 0.9;
epsilon = 10;

assetReserve = ones(n,1)*100;
assetReserve = [1,20,20];
b = taxParam ./ (assetReserve + 0.01);

[selection, trade, profit] = TACo(C, b, d0, gamma, epsilon);
selection
trade = trade(:,selection)
postCost = -profit(:,selection)
preCost = C(:,selection)

% Compute shortfall
shortfalls = zeros(n,1);
for i = 1:n
    shortfalls(i) = max(0, -trade(i) - assetReserve(i));
end
disp("Shortfalls: "+num2str(shortfalls))
disp("Normalized Shortfalls: "+num2str(shortfalls/sum(shortfalls)))