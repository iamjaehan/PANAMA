% Load the results from the .mat file
data = load("MC_test_results.mat");
% The structure fields:
% data.assetReserve, data.taxParam, data.repeat, data.rounds, data.shortFall

% Example: Display all rounds for each run
disp('Rounds for each run:');
disp(data.rounds);

% Plot the shortfall history
% n_runs = length(data.shortFall_history);
% figure(1)
% clf
% for i = 1:n_runs
%     plot(data.shortFall_history{i}, '-o');
%     hold on;
% end
% legend(arrayfun(@(i) sprintf('Run %d', i), 1:n_runs, 'UniformOutput', false));
% xlabel('Round');
% ylabel('Shortfall');
% title('Shortfall History');
% set(gca, 'FontSize', 15);

% Plot rounds for each case
figure(2)
clf
plot(data.rounds)
xlabel('Round');
ylabel('Rounds');
title('Rounds for each case');
set(gca, 'FontSize', 15);