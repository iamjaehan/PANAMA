% Load the results from the .mat file
data = load('MC_test_results_tight.mat');  % or 'MC_test_results.mat'
data = load("MC_test_results.mat");
data = load("MC_test_results_merged.mat"); data = data.mergedData;
% The structure fields:
% data.assetReserve, data.taxParam, data.repeat, data.rounds, data.shortFall

% Example: Display all rounds for each run
disp('Rounds for each run:');
disp(data.rounds);

% Example: Display assetReserve and taxParam for each run
for i = 1:length(data.rounds)
    fprintf('Run %d: assetReserve = [%s], taxParam = %d, rounds = %d\n', ...
        i, num2str(data.assetReserve{i}), data.taxParam(i), data.rounds(i));
end

% Example: Plot shortFall for each run
figure;
for i = 1:length(data.shortFall)
    plot(data.shortFall{i}, '-o');
    hold on;
end
xlabel('Index');
ylabel('ShortFall');
title('ShortFall for each run');
legend(arrayfun(@(i) sprintf('Run %d', i), 1:length(data.shortFall), 'UniformOutput', false));