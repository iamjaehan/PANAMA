% Load the results from the .mat file
% data = load("MC_test_results.mat");
% data = load("MC_test_results_randomSampling_500.mat");
% data = load("MC_test_results_randomSampling_1000.mat");
% data = load("MC_test_results_randomSampling_10000.mat");
data = load("MC_test_results_randomSampling_1000_w_baselines.mat");
% The structure fields:
% data.assetReserve, data.taxParam, data.repeat, data.rounds, data.shortFall

%% Pre processing
dataLen = length(data.assetReserve);
n = length(data.assetReserve{1});

% Fuck..
for i = 1:dataLen
    data.taxParam(i) = 1/data.taxParam(i);
end

choiceCostDispersionHistory = cell(dataLen,1);
choiceCostHistory = cell(dataLen,1);

choiceCost = zeros(dataLen, n);
choiceCostInit = zeros(dataLen, n);
finalProfit = zeros(dataLen, n);
systemCost = zeros(dataLen,1);
systemVar = zeros(dataLen,1);
systemStd = zeros(dataLen,1);
for i = 1:dataLen
    % Initial data
    tempNegoInitOut = data.negoOut_history{i}{1};
    tempInitOutcome = tempNegoInitOut.outcome;

    % Final data
    tempNegoOut = data.negoOut_history{i}{end};
    tempOutcome = tempNegoOut.outcome;

    % Process
    choiceCostInit(i,:) = tempNegoInitOut.C(:, tempInitOutcome);
    choiceCost(i,:) = tempNegoOut.C(:,tempOutcome);
    finalProfit(i,:) = tempNegoOut.profit(:,tempOutcome);
    
    systemCost(i) = sum(choiceCost(i,:)) / sum(choiceCostInit(i,:));
    % systemStd(i) = std(choiceCost(i,:));% / std(choiceCostInit(i,:));
    systemVar(i) = ComputeGini(choiceCost(i,:)) / ComputeGini(choiceCostInit(i,:));
    systemStd(i) = systemVar(i);

    % History data
    numRounds = data.rounds(i);
    for j = 1:numRounds
        choiceCostDispersionHistory{i}(j) = mean(arrayfun(@(x) ComputeGini(data.negoOut_history{i}{j}.C(:,x)), [1,2,3]));
        % choiceCostDispersionHistory{i}(j) = ComputeGini(std(data.negoOut_history{i}{j}.C));
        choiceCostHistory{i}(j) = mean(sum(data.negoOut_history{i}{j}.C));
        % temp = zeros(3,1);
        % for k = 1:3
        %     temp(k) = ComputeGini(data.negoOut_history{i}{j}.C(:,k));
        % end
        % choiceCostHistory{i}(j) = mean(temp);
    end
end

% Define krenel
assetReserve = zeros(dataLen,n);
reserveVar = zeros(dataLen,1);
bVar = zeros(dataLen,1);
meanb = zeros(dataLen,1);
reserveGini = zeros(dataLen,1);
for i = 1:dataLen
    assetReserve(i,:) = double(data.assetReserve{i}); 
    reserveVar(i) = std(assetReserve(i,:)/max(assetReserve(i,:)));
    % reserveVar(i) = std(assetReserve(i,:));
    reserveGini(i) = ComputeGini(assetReserve(i,:)/max(assetReserve(i,:)));
    bVar(i) = std(data.negoOut_history{i}{1}.b/max(data.negoOut_history{i}{1}.b));
    
    % meanb(i) = min(data.negoOut_history{i}{1}.b)/max(data.negoOut_history{i}{1}.b);
    % meanb(i) = min(data.negoOut_history{i}{1}.b)/data.taxParam(i);
    
    meanb(i) = mean(data.negoOut_history{i}{1}.b);
    % meanb(i) = mean(data.negoOut_history{i}{1}.b)/data.taxParam(i);
    % meanb(i) = std(data.negoOut_history{i}{1}.b)/mean(data.negoOut_history{i}{1}.b);
    % meanb(i) = ComputeGini(double(data.assetReserve{i}));
    % bVar(i) = ComputeGini(data.negoOut_history{i}{1}.b);
    
    % meanb(i) = mean(double(data.assetReserve{i})./data.negoOut_history{i}{1}.b);
end

Q1 = prctile(data.taxParam, 25);
Q2 = prctile(data.taxParam, 75);

K1 = find(data.taxParam < Q1);
K2 = find(data.taxParam >= Q1 & data.taxParam < Q2);
K3 = find(data.taxParam >= Q2);

% Q1 = prctile(reserveGini, 33);
% Q2 = prctile(reserveGini, 66);
% 
% Z1 = find(reserveGini < Q1);
% Z2 = find(reserveGini >= Q1 & reserveGini < Q2);
% Z3 = find(reserveGini >= Q2);

Q1 = prctile(meanb, 25);
Q2 = prctile(meanb, 75);

Z1 = find(meanb < Q1);
Z2 = find(meanb >= Q1 & meanb < Q2);
Z3 = find(meanb >= Q2);



colors = lines(10);
% Exp 1. Rounds per k
figure(1)
clf
semilogx(data.taxParam(Z1), data.rounds(Z1), 'o', 'Color',[colors(1,:)])
hold on
semilogx(data.taxParam(Z2), data.rounds(Z2), '*', 'Color',[colors(2,:)])
semilogx(data.taxParam(Z3), data.rounds(Z3), 's', 'Color',[colors(3,:)],'MarkerSize',10)
legend({"Low asset value","Med asset value","High asset value"})
xlabel('tax param k');
ylabel('Rounds');
title('Rounds vs k');
grid on
set(gca, 'FontSize', 15);

% Exp 2. Sys-cost per k
figure(2)
clf
semilogx(data.taxParam(Z1), systemCost(Z1), 'o', 'Color',[colors(1,:)])
hold on
semilogx(data.taxParam(Z2), systemCost(Z2), '*', 'Color',[colors(2,:)])
semilogx(data.taxParam(Z3), systemCost(Z3), 's', 'Color',[colors(3,:)],'MarkerSize',8)
legend({"Low asset value","Med asset value","High asset value"})
xlabel("tax Param k");
ylabel("System Cost")
title("System Cost vs k")
grid on
set(gca, 'FontSize', 15);

% % Exp 2-2. Sys-cost per k
% figure(22)
% clf
% lastValues = arrayfun(@(i) choiceCostDispersionHistory{i}(end), Z1);
% semilogx(data.taxParam(Z1), lastValues, 'o', 'Color',[colors(1,:)])
% hold on
% lastValues = arrayfun(@(i) choiceCostDispersionHistory{i}(end), Z2);
% semilogx(data.taxParam(Z2), lastValues, '*', 'Color',[colors(2,:)])
% lastValues = arrayfun(@(i) choiceCostDispersionHistory{i}(end), Z3);
% semilogx(data.taxParam(Z3), lastValues, 's', 'Color',[colors(3,:)],'MarkerSize',8)
% legend({"Low asset value","Med asset value","High asset value"})
% xlabel("tax Param k");
% ylabel("C Cost")
% title("C cost per k")
% grid on
% set(gca, 'FontSize', 15);

% Exp 3. Cost-var per k
figure(3)
clf
semilogx(data.taxParam(Z1), systemStd(Z1), 'o', 'Color',[colors(1,:)])
hold on
semilogx(data.taxParam(Z2), systemStd(Z2), '*', 'Color',[colors(2,:)])
semilogx(data.taxParam(Z3), systemStd(Z3), 's', 'Color',[colors(3,:)],'MarkerSize',8)
legend({"Low asset value","Med asset value","High asset value"})
xlabel("tax Param k");
ylabel("System Cost Gini Index")
title("System Cost Gini Index vs k")
grid on
set(gca, 'FontSize', 15);
% 
% % Exp 11. Average Purchase power to round
% figure(11)
% clf
% plot(meanb, data.rounds, 'o', 'Color',[colors(1,:)])
% legend({"Low asset value","Med asset value","High asset value"})
% xlabel('Purchase powers');
% ylabel('Rounds');
% title('Rounds vs Purchase power');
% grid on
% set(gca, 'FontSize', 15);
% 
% % Exp 12.
% uniqueRounds = unique(data.rounds);
% edges = 0:0.05:1;  % purchase power의 bin
% 
% Z = zeros(length(uniqueRounds), length(edges)-1);  % histogram 저장용
% 
% for i = 1:length(uniqueRounds)
%     r = uniqueRounds(i);
%     idx = data.rounds == r;
%     Z(i, :) = histcounts(meanb(idx), edges);
% end
% 
% % Plot
% figure(12)
% clf
% [X, Y] = meshgrid(edges(1:end-1) + 0.025, uniqueRounds);  % bin center
% surf(X, Y, Z);
% xlabel('Purchase power');
% ylabel('Rounds');
% zlabel('Frequency');
% title('Histogram of purchase power by round');
% view(45, 30);  % 보기 각도 조절
% colorbar;
% 
% % Exp 13 (log10-transformed histogram for better log visualization)
% uniqueRounds = unique(data.rounds);
% log_edges = linspace(-2, 1, 40);  % log10(0.01) to log10(10), 40 bins
% 
% Z = zeros(length(uniqueRounds), length(log_edges)-1);  % histogram 저장용
% 
% % log10 취한 값으로 histogram 계산
% log_taxParam = log10(data.taxParam);
% 
% for i = 1:length(uniqueRounds)
%     r = uniqueRounds(i);
%     idx = data.rounds == r;
%     Z(i, :) = histcounts(log_taxParam(idx), log_edges);
% end
% 
% % bin 위치 설정 (log10 scale 상에서)
% [X, Y] = meshgrid(log_edges(1:end-1), uniqueRounds);
% 
% % 시각화
% figure(13); clf;
% surf(X, Y, Z);
% xlabel('Tax param');
% ylabel('Rounds');
% zlabel('Frequency');
% title('Histogram of purchase power by round (log-transformed tax param)');
% view(45, 30);
% colorbar;
% 
% % x축 눈금 설정 (log tick처럼 보이게)
% set(gca, 'XTick', [-2 -1 0 1], ...
%          'XTickLabel', {'10^{-2}', '10^{-1}', '10^{0}', '10^{1}'});
% 
% %%
% figure(4)
% clf
% for i = 1:dataLen
%     semilogx(data.taxParam(i), choiceCostDispersionHistory{i}(end)/choiceCostDispersionHistory{i}(1), 'o', 'Color',[colors(1,:)])
%     hold on
% end
% xlabel("tax Param k");
% ylabel("Choice Cost Dispecsion ")
% title("Choice Cost Dispersion vs k")
% grid on
% set(gca, 'FontSize', 15);
% 
% 
% % Exp 5. Rounds per Asset Reserve variation
% figure(5)
% clf
% for i = 1:dataLen
%     semilogx(data.taxParam(i), choiceCostHistory{i}(end)/choiceCostHistory{i}(1), 'o', 'Color',[colors(1,:)])
%     hold on
% end
% xlabel("tax Param k");
% ylabel("Choice Cost Variation")
% title("Choice Cost vs k")
% grid on
% set(gca, 'FontSize', 15);
% 
% %%
% % Exp 7. Evolution of choice cost variation
% figure(7)
% clf
% for i = 1:dataLen
%     plot(1:data.rounds(i), choiceCostDispersionHistory{i}/choiceCostDispersionHistory{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     % plot(1:data.rounds(i), choiceCostHistory{i},'Color',[0,0,0,0.2],'LineWidth',3)
%     plot(data.rounds(i), choiceCostDispersionHistory{i}(end)/choiceCostDispersionHistory{i}(1), 'r*')
%     hold on
% end
% grid on
% title("Choice Cost Dispersion vs Rounds")
% xlabel("Rounds")
% ylabel("Choice Cost")
% set(gca,'FontSize',18)
% xlim([1 inf])
% 
% % Exp 8. Evolution of short fall
% figure(8)
% clf
% for i = Z1'
%     plot(1:data.rounds(i), choiceCostHistory{i}/choiceCostHistory{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), choiceCostHistory{i}(end)/choiceCostHistory{i}(1), 'o', 'Color',colors(1,:))
%     hold on
% end
% for i = Z2'
%     plot(1:data.rounds(i), choiceCostHistory{i}/choiceCostHistory{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), choiceCostHistory{i}(end)/choiceCostHistory{i}(1), '*','Color',colors(2,:))
%     hold on
% end
% for i = Z3'
%     plot(1:data.rounds(i), choiceCostHistory{i}/choiceCostHistory{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), choiceCostHistory{i}(end)/choiceCostHistory{i}(1), 's','Color',colors(3,:),'MarkerSize',8)
%     hold on
% end
% grid on
% title("Choice Cost vs Rounds")
% xlabel("Rounds")
% ylabel("Choice Cost")
% set(gca,'FontSize',18)
% xlim([1 inf])
% 
% % figure(8)
% % clf
% % for i = 1:dataLen
% %     localShortFallHistory = data.shortFall_history{i};
% %     histLen = length(localShortFallHistory);
% %     shortFallHist = zeros(histLen,n);
% %     for j = 1:histLen
% %         % localShortFall = data.shortFall_history{i}{j}./double(data.assetReserve{i});
% %         localShortFall = data.shortFall_history{i}{j};
% %         shortFallHist(j,:) = localShortFall;
% %     end
% %     for j = 1:n
% %         plot(1:data.rounds(i), shortFallHist(:,j),'Color',[0,0,0,0.2],'LineWidth',3)
% %     end
% %     hold on
% % end
% % grid on
% % title("Normalized ShortFall vs Rounds")
% % xlabel("Rounds")
% % ylabel("Choice Cost")
% % set(gca,'FontSize',18)
% % xlim([1 inf])
% 
% % Exp 9. System cost variation
% % figure(9)
% % clf
% % for i = 1:dataLen
% %     rounds = double(data.rounds(i));
% %     localSystemStd = zeros(rounds,1);
% %     localSystemCost = zeros(rounds,1);
% %     for j = 1:rounds
% %         localC = data.negoOut_history{i}{j}.C;
% %         localOutcome = data.negoOut_history{i}{j}.outcome;
% %         localSystemCostVec = localC(:,localOutcome);
% %         % localSystemStd(j) = std(localSystemCostVec);
% %         localSystemStd(j) = ComputeGini(localSystemCostVec);
% %         % localSystemCost(j) = sum(localSystemCostVec);
% %     end
% %     % plot(1:rounds,localSystemStd,'Color',[0,0,0,0.2])
% %     plot(1:rounds,localSystemStd/localSystemStd(1),'Color',[0,0,0,0.05],'LineWidth',3)
% %     % plot(rounds,localSystemStd(end),'r*')
% %     plot(rounds, localSystemStd(end)/localSystemStd(1),'r*')
% %     hold on
% % end
% % grid on
% % title("System Cost Gini vs Rounds")
% % xlabel("Rounds")
% % set(gca,'FontSize',18)
% % xlim([1 inf])
% %%
% figure(9)
% clf
% SystemStd = cell(dataLen,1);
% for i = 1:dataLen
%     rounds = double(data.rounds(i));
%     localSystemStd = zeros(rounds,1);
%     for j = 1:rounds
%         localC = data.negoOut_history{i}{j}.C;
%         localOutcome = data.negoOut_history{i}{j}.outcome;
%         localSystemCostVec = localC(:,localOutcome);
%         % localSystemStd(j) = std(localSystemCostVec);
%         localSystemStd(j) = ComputeGini(localSystemCostVec);
%         % localSystemCost(j) = sum(localSystemCostVec);
%     end
%     SystemStd{i} = localSystemStd;
% end
% for i = Z1'
%     plot(1:data.rounds(i), SystemStd{i}/SystemStd{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), SystemStd{i}(end)/SystemStd{i}(1), 'o', 'Color',colors(1,:))
%     hold on
% end
% for i = Z2'
%     plot(1:data.rounds(i), SystemStd{i}/SystemStd{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), SystemStd{i}(end)/SystemStd{i}(1), '*','Color',colors(2,:))
%     hold on
% end
% for i = Z3'
%     plot(1:data.rounds(i), SystemStd{i}/SystemStd{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), SystemStd{i}(end)/SystemStd{i}(1), 's','Color',colors(3,:),'MarkerSize',8)
%     hold on
% end
% grid on
% title("Relative System Cost Gini vs Rounds")
% xlabel("Rounds")
% ylabel("Relative System Cost Gini Index")
% set(gca,'FontSize',18)
% xlim([1 inf])
% 
% 
% % Exp 10. System cost
% % figure(10)
% % clf
% % for i = 1:dataLen
% %     rounds = double(data.rounds(i));
% %     localSystemCost = zeros(rounds,1);
% %     for j = 1:rounds
% %         localC = data.negoOut_history{i}{j}.C;
% %         localOutcome = data.negoOut_history{i}{j}.outcome;
% %         localSystemCostVec = localC(:,localOutcome);
% %         localSystemCost(j) = sum(localSystemCostVec);
% %     end
% %     plot(1:rounds,localSystemCost/localSystemCost(1),'Color',[0,0,0,0.05],'LineWidth',3)
% %     % plot(1:rounds,localSystemCost,'Color',[0,0,0,0.2],'LineWidth',3)
% %     plot(rounds, localSystemCost(end)/localSystemCost(1),'r*')
% %     hold on
% % end
% % grid on
% % title("System Cost vs Rounds")
% % xlabel("Rounds")
% % set(gca,'FontSize',18)
% % xlim([1 inf])
% 
% figure(10)
% clf
% SystemCost = cell(dataLen,1);
% for i = 1:dataLen
%     rounds = double(data.rounds(i));
%     localSystemCost = zeros(rounds,1);
%     for j = 1:rounds
%         localC = data.negoOut_history{i}{j}.C;
%         localOutcome = data.negoOut_history{i}{j}.outcome;
%         localSystemCostVec = localC(:,localOutcome);
%         localSystemCost(j) = sum(localSystemCostVec);
%     end
%     SystemCost{i} = localSystemCost;
% end
% for i = Z1'
%     plot(1:data.rounds(i), SystemCost{i}/SystemCost{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), SystemCost{i}(end)/SystemCost{i}(1), 'o', 'Color',colors(1,:))
%     hold on
% end
% for i = Z2'
%     plot(1:data.rounds(i), SystemCost{i}/SystemCost{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), SystemCost{i}(end)/SystemCost{i}(1), '*','Color',colors(2,:))
%     hold on
% end
% for i = Z3'
%     plot(1:data.rounds(i), SystemCost{i}/SystemCost{i}(1),'Color',[0,0,0,0.05],'LineWidth',3)
%     plot(data.rounds(i), SystemCost{i}(end)/SystemCost{i}(1), 's','Color',colors(3,:),'MarkerSize',8)
%     hold on
% end
% grid on
% title("Relative System Cost vs Rounds")
% xlabel("Rounds")
% ylabel("Relative System Cost")
% set(gca,'FontSize',18)
% xlim([1 inf])
% 
% %%
% % Exp 11. b value to preference
% % We have b value --- connect it to 기회비용?
% % 기회비용 = Cinit - SelectedCost
% 
% figure(23)
% clf
% b2OppCost = zeros(dataLen,1);
% for i = 1:dataLen
%     Cinit = data.negoOut_history{i}{1}.C; % Initial selfish cost
%     Cfinal = data.negoOut_history{i}{end}.C; % Final cost
%     proposalCost = diag(Cfinal);
%     % initialOutcome = data.negoOut_history{i}{1}.outcome;
%     % proposalCost = Cinit(:,initialOutcome);
% 
%     finalOutcome = data.negoOut_history{i}{end}.outcome;
%     % oppCost = mean(Cfinal(:,finalOutcome) - proposalCost);
%     oppCost = Cfinal(:,finalOutcome) - proposalCost;
%     % oppCost = mean(diag(Cfinal) - proposalCost);
% 
%     bVec = data.negoOut_history{i}{1}.b;
%     assetReserve = double(data.assetReserve{i});
% 
%     % localK = data.taxParam(i);
%     % localK = data.rounds(i);
%     % localK = sum(Cfinal(:,finalOutcome));
%     % localK = normalize(bVec,'norm');
%     localK = bVec/max(bVec);
% 
%     b2OppCost(i) = corr(localK,oppCost);
% 
%     [~,I] = sort(localK); 
%     plot(localK(I), oppCost(I),'.-','Color',[0 0 0 0.02])
%     hold on
% end
% grid on
% xlabel("Normalized asset value")
% ylabel("Opportunity Cost")
% set(gca,'FontSize',15)
% 
% figure(24)
% clf
% for i = 1:dataLen
%     Cfinal = data.negoOut_history{i}{end}.C; % Final cost
%     outcome = data.negoOut_history{i}{end}.outcome;
%     bVec = data.negoOut_history{i}{1}.b;
%     assetReserve = double(data.assetReserve{i});
% 
%     [~,I] = sort(bVec);
%     plot(normalize(bVec(I),'norm'), Cfinal(I,outcome),'.-','Color',[0 0 0 0.02])
%     % plot(bVec, Cfinal(:,outcome),'o','Color',[0 0 0 0.2])
%     hold on
% end
% grid on
% xlabel("Normalized asset value")
% ylabel("Individual Cost")
% set(gca,'FontSize',15)
% %%
% 
% figure(25)
% clf
% b2CostCorr = zeros(dataLen,1);
% for i = 1:dataLen
%     Cfinal = data.negoOut_history{i}{end}.C; % Final cost
%     outcome = data.negoOut_history{i}{end}.outcome;
%     finalCost = Cfinal(:,outcome);
%     bVec = data.negoOut_history{i}{1}.b;
%     assetReserve = double(data.assetReserve{i});
%     k = data.taxParam(i);
% 
%     b2CostCorr(i) = corr(finalCost, normalize(bVec,'norm'));
% end
% histogram(b2OppCost,'BinWidth',0.02)
% % histogram(b2CostCorr,'BinWidth',0.02)
% title("Opportunity Cost Correlation to asset value Histogram")
% xlabel("Correlation Coefficient")
% ylabel("Counts")
% grid on
% set(gca, 'fontsize', 15)
% 
% %% Comparison Study
% cent_SystemCost = zeros(dataLen,1);
% fcfs_SystemCost = zeros(dataLen,1);
% cent_SystemStd = zeros(dataLen,1);
% fcfs_SystemStd = zeros(dataLen,1);
% for i = 1:dataLen
%     cent_SystemCost(i) = sum(cell2mat(data.centralized_cost{i}));
%     cent_SystemStd(i) = ComputeGini(cell2mat(data.centralized_cost{i}));
%     fcfs_SystemCost(i) = sum(cell2mat(data.fcfs_cost{i}));
%     fcfs_SystemStd(i) = ComputeGini(cell2mat(data.fcfs_cost{i}));
% end
% 
% ours_SystemCost = zeros(dataLen,1);
% ours_SystemStd = zeros(dataLen,1);
% for i = 1:dataLen
%     ours_SystemCost(i) = SystemCost{i}(end);
%     ours_SystemStd(i) = SystemStd{i}(end);
% end
% 
% CostDataPack = [ours_SystemCost, cent_SystemCost, fcfs_SystemCost];
% StdDataPack = [ours_SystemStd, cent_SystemStd, fcfs_SystemStd];
% 
% test = [2.5399535140875003, 0.5702869889284113, 0.18135935248821328];
% naive_SystemStd = zeros(dataLen,1);
% naive_SystemCost = zeros(dataLen,1);
% for i = 1:dataLen
%     naive_SystemStd(i) = ComputeGini(test);
%     naive_SystemCost(i) = sum(test);
% end

%% Realtion Study
Bmax = zeros(dataLen,1);
for i = 1:dataLen
    Bmax(i) = ComputeBMax(data,i);
end
roundBound = zeros(dataLen,1);
for i = 1:dataLen
    % roundBound(i) = data.taxParam(i) * Bmax(i);
    roundBound(i) = ceil(data.taxParam(i) * Bmax(i));
end

figure(41)
clf
semilogx(roundBound(Z1), data.rounds(Z1), 'o', 'Color',[colors(1,:)])
hold on
semilogx(roundBound(Z2), data.rounds(Z2), '*', 'Color',[colors(2,:)])
semilogx(roundBound(Z3), data.rounds(Z3), 's', 'Color',[colors(3,:)],'MarkerSize',10)
semilogx(1:30, 1:30, 'k')
legend({"Low asset value","Med asset value","High asset value"})
xlabel('Round Bound');
ylabel('Rounds');
title('Rounds vs k');
grid on
set(gca, 'FontSize', 15);

effBound = zeros(dataLen,1);
for i = 1:dataLen
    % effBound(i) = n*Bmax(i)/roundBound(i) + n*Bmax(i)*2*(1-1/3);
    effBound(i) = n*Bmax(i)*(2*(1-1/n) + 1/data.rounds(i));
end
optCost = zeros(dataLen,1);
optGap = zeros(dataLen,1);
for i = 1:dataLen
    optCost(i) = min(sum(data.negoOut_history{i}{end}.C));
    optGap(i) = sum(choiceCost(i,:)) - optCost(i);
end
figure(42)
clf
semilogx(effBound(Z1), optGap(Z1), 'o', 'Color',[colors(1,:)])
hold on
semilogx(effBound(Z2), optGap(Z2), '*', 'Color',[colors(2,:)])
semilogx(effBound(Z3), optGap(Z3), 's', 'Color',[colors(3,:)],'MarkerSize',10)
semilogx(1:30, 1:30, 'k')
legend({"Low asset value","Med asset value","High asset value"})
xlabel('Efficiency Bound');
ylabel('Opt Gap');
title('k vs opt gap');
grid on
set(gca, 'FontSize', 15);

%% Comparison Study

figure(50)
clf
group1 = ours_SystemCost(K1);
group2 = ours_SystemCost(K2);
group3 = ours_SystemCost(K3);
group4 = cent_SystemCost;
group5 = fcfs_SystemCost;
group6 = naive_SystemCost;
localPack = [group1;group2;group3;group4;group5;group6];
group = [ones(length(group1),1); repmat(2,length(group2),1); repmat(3,length(group3),1); repmat(4,dataLen,1); repmat(5,dataLen,1); repmat(6,dataLen,1)];
% boxplot(CostDataPack,'Labels',{'Ours', 'Centralized-CTOP', 'FCFS-CTOP'})
boxplot(localPack,group,'Labels',{'Ours (Low k)', 'Ours (Med k)', 'Ours (High k)', 'Centralized-CTOP', 'FCFS-CTOP', 'Naive-CTOP'})
title("System Cost Comparison")
ylabel("Cost")
grid on
set(gca,'fontsize',15,'FontWeight','bold')

figure(51)
clf
group1 = ours_SystemStd(K1);
group2 = ours_SystemStd(K2);
group3 = ours_SystemStd(K3);
group4 = cent_SystemStd;
group5 = fcfs_SystemStd;
group6 = naive_SystemStd;
localPack = [group1;group2;group3;group4;group5;group6];
group = [ones(length(group1),1); repmat(2,length(group2),1); repmat(3,length(group3),1); repmat(4,dataLen,1); repmat(5,dataLen,1); repmat(6,dataLen,1)];
boxplot(localPack,group,'Labels',{'Ours (Low k)', 'Ours (Med k)', 'Ours (High k)', 'Centralized-CTOP', 'FCFS-CTOP', 'Naive-CTOP'})
title("System Gini Index Comparison")
ylabel("Gini Index")
grid on
set(gca,'fontsize',15,'FontWeight','bold')

%% Define function
function out = ComputeBMax(data, idx)
    C = data.negoOut_history{idx}{1}.C;
    agentMax = abs(max(C,[],2) - min(C,[],2));
    out = max(agentMax);
end