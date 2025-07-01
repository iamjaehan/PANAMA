function outcome = TACo(C, b, d0, gamma, epsilon)
% Inputs:
%   C      : [n x m] cost matrix
%   b      : [n x 1] private valuations
%   d0     : initial trading unit (scalar)
%   gamma  : trading unit reduction factor (0 < gamma < 1)
%   epsilon: convergence threshold

[n, m] = size(C);
O = zeros(n, m);           % Offer matrix
P = zeros(n, m);           % Pay matrix
d = d0;                    % Trading unit
isConverged = false;       % Termination flag
selections = zeros(n, 1);  % Track each agent's latest choice
recordedStates = {};       % For cycle detection
agentOrder = 1:n;          % Fixed sequential agent order
step = 0;                  % Step counter

% fprintf(' Step | Agent |   O   |   P   |     J (profit matrix)     | Selections\n');
% fprintf('--------------------------------------------------------------------------\n');

while ~isConverged
    step = step + 1;
    i = agentOrder(mod(step-1, n) + 1);  % Sequential agent
    
    % --- Profit computation ---
    J = b(i) * (O(i,:) - P(i,:)) - C(i,:);   % Profit vector for agent i
    [~, j_star] = max(J);                   % Best choice for agent i
    
    % --- Update matrices ---
    P(i, j_star) = P(i, j_star) + n * d;     % Agent i pays
    O(:, j_star) = O(:, j_star) + d;        % Everyone receives offer
    selections(i) = j_star;                 % Update selection
    
    % --- Display formatted output ---
    % fprintf('  %2d   |   %d    | ', step, i);
    % fprintf('[%s] | ', strjoin(arrayfun(@(x) sprintf('%d', x), O(i,:), 'UniformOutput', false), ' '));
    % fprintf('[%s] | ', strjoin(arrayfun(@(x) sprintf('%d', x), P(i,:), 'UniformOutput', false), ' '));
    % fprintf('[%s] | ', strjoin(arrayfun(@(x) sprintf('%5.1f', x), J, 'UniformOutput', false), ' '));
    % fprintf('[%s]\n', strjoin(arrayfun(@(x) sprintf('%d', x), selections(~isnan(selections)), 'UniformOutput', false), ' '));
    
    % --- Cycle detection ---
    key = [reshape(O - P, 1, []), i];        % Flattened key
    found = any(cellfun(@(x) isequal(x, key), recordedStates));
    
    if found
        d = gamma * d;                      % Reduce trading unit
        recordedStates = {};               % Reset recorded states
        
        % Check epsilon-termination condition
        allSatisfied = true;
        for k = 1:n
            Jk = b(k) * (O(k,:) - P(k,:)) - C(k,:);
            if max(Jk) - min(Jk) <= epsilon
                allSatisfied = false;
                break;
            end
        end
        
        if allSatisfied
            isConverged = true;
        end
    end
    
    % --- Record new state ---
    recordedStates{end+1} = key;
end

% Return most frequently selected outcome
outcome = mode(selections);

end
