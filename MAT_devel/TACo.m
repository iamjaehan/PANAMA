function [outcome, trade, profit, step, maxCycleLen] = TACo(C, b, d0, gamma, epsilon)
% Inputs:
%   C      : [n x m] cost matrix
%   b      : [n x 1] private valuations
%   d0     : initial trading unit (scalar)
%   gamma  : trading unit reduction factor (0 < gamma < 1)
%   epsilon: convergence threshold
%   U: trade limit

U = 100; %#ok<NASGU>  % 아직 안 쓰고 있음

[n, m] = size(C);
O = zeros(n, m);           % Offer matrix
P = zeros(n, m);           % Pay matrix
d = d0;                    % Trading unit
isConverged = false;       % Termination flag
selections = zeros(n, 1);  % Track each agent's latest choice

agentOrder = 1:n;          % Fixed sequential agent order
step = 0;                  % Step counter
maxCycleLen = 0;           % 최대 cycle 길이

% 상태(key: char) -> 처음 등장한 step(double)을 저장하는 map
stateMap = containers.Map('KeyType','char','ValueType','double');

% --- cycle 동치 판단을 위한 양자화 해상도 ---
% tol보다 작은 차이는 같은 state로 본다.
% d0 ~ 1, epsilon ~ 1e-1 ~ 1e-2 정도라면 이 정도면 충분히 보수적이다.
tol = 1e-8;

while ~isConverged
    step = step + 1;
    i = agentOrder(mod(step-1, n) + 1);  % Sequential agent
    
    % --- Profit computation ---
    J = b(i) * (O(i,:) - P(i,:)) - C(i,:);   % Profit vector for agent i
    [~, j_star] = max(J);                    % Best choice for agent i
    
    % --- Update matrices ---
    P(i, j_star) = P(i, j_star) + n * d;     % Agent i pays
    O(:, j_star) = O(:, j_star) + d;         % Everyone receives offer
    selections(i) = j_star;                  % Update selection

    % --- Cycle detection key (현재 상태 + agent index) ---
    % 1) 연속 상태 벡터 구성
    vec = [reshape(O - P, 1, []), i];        % double 벡터

    % 2) tol 격자 위로 양자화: vec_q는 정수 벡터
    vec_q = round(vec / tol);

    % 3) 정수 벡터를 문자열로 변환해서 key 생성
    key = mat2str(vec_q);                    % full precision 대신 정수라 truncation 걱정 없음

    % 이전에 같은 (양자화된) 상태가 있었는지 확인
    found = isKey(stateMap, key);
    
    % 전원 같은 선택이면 종료 (원래 로직 그대로)
    if all(selections == selections(1))
        break;
    end
    
    if found
        % disp("found!")
        % cycle 길이: 현재 step - 처음 그 상태를 봤던 step
        firstStep = stateMap(key);
        cycleLen  = step - firstStep;

        if cycleLen > maxCycleLen
            maxCycleLen = cycleLen;
        end

        % trading unit 감소 및 state 기록 초기화 (원래 로직 유지)
        d = gamma * d;                 
        stateMap = containers.Map('KeyType','char','ValueType','double');
        
        % epsilon-termination condition 체크 (네가 써둔 형태 그대로)
        allSatisfied = true;
        for k = 1:n
            Jk = b(k) * (O(k,:) - P(k,:)) - C(k,:);
            if max(Jk) - min(Jk) > epsilon
                allSatisfied = false;
                break;
            end
        end
        
        if allSatisfied
            isConverged = true;
        end
    end
    
    % --- Record new state (현재 상태를 기록) ---
    stateMap(key) = step;
end

% Return most frequently selected outcome
outcome = mode(selections);
trade   = O - P;
profit  = diag(b) * (O - P) - C;

end
