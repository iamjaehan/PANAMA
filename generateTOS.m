function TOS = generateTOS(start_idx, end_idx)

%% Compute K-Shortest Paths using Yen’s Algorithm (DISTINCT PATHS)
G = load("US_waypoint_graph.mat");
G = G.G;
node_coords = G.Nodes{:,2:3};

K = 5; % Number of alternative paths
all_paths = cell(K,1); % Store all paths
all_path_distances = zeros(K,1);

% Compute the first shortest path using Dijkstra
% start_idx = 2887; % Origin node index
% start_idx = 31;
% end_idx = 4127; % Austin arrival

% end_idx = num_samples + 2; % Destination node index
[all_paths{1}, all_path_distances(1)] = shortestpath(G, start_idx, end_idx);

% Find alternative paths ensuring distinct routes
for k = 2:K
    max_attempts = 5; % Maximum retries to find a unique path
    attempt = 0;
    unique_path_found = false;
    
    while ~unique_path_found && attempt < max_attempts
        attempt = attempt + 1;
        
        % Create a temporary graph by increasing the weight of the previous paths
        G_temp = G;
        
        % Penalize edges used in the previous paths
        for prev_k = 1:k-1
            for i = 1:length(all_paths{prev_k})-1
                edge_idx = findedge(G_temp, all_paths{prev_k}(i), all_paths{prev_k}(i+1));
                if edge_idx > 0
                    G_temp.Edges.Weight(edge_idx) = G_temp.Edges.Weight(edge_idx) * (5 * 2^attempt); % Increase penalty progressively
                end
            end
        end
        
        % Compute the next shortest path
        [candidate_path, candidate_distance] = shortestpath(G_temp, start_idx, end_idx);
        
        % Check if the new path is distinct
        is_duplicate = any(cellfun(@(p) isequal(p, candidate_path), all_paths(1:k-1)));
        
        % Accept path if unique
        if ~is_duplicate
            all_paths{k} = candidate_path;
            all_path_distances(k) = candidate_distance;
            unique_path_found = true;
        end
    end
    
    % If no unique path found after max_attempts, break out
    if ~unique_path_found
        % warning('Could not find %d distinct paths. Stopping at %d.', K, k-1);
        K = k-1; % Adjust K to reflect actual number of paths found
        break;
    end
end



%% Plot results
if false
figure(2);
clf;
% geoplot(wplat, wplon, 'm.','MarkerSize',2); % Plot all waypoints
geoplot(node_coords(:,1), node_coords(:,2), 'b.', 'MarkerSize', 3); % Graph nodes
hold on
geoplot(node_coords(start_idx,1),node_coords(start_idx,2),'r*','MarkerSize',8,'LineWidth',5)
geoplot(node_coords(end_idx,1),node_coords(end_idx,2),'g*','MarkerSize',8,'LineWidth',5)

% Plot multiple paths with different colors
colors = {'r-', 'g-', 'c-','m-','y-'};
for k = 1:K
    path_lats = G.Nodes.Lat(all_paths{k});
    path_lons = G.Nodes.Lon(all_paths{k});
    % geoplot(path_lats, path_lons, colors{k}, 'LineWidth', 2);
    geoplot(path_lats, path_lons, '-', 'LineWidth', 4, 'Color', [rand(1,3) 0.5]);
end

geoplot(node_coords(start_idx,1),node_coords(start_idx,2),'r*','MarkerSize',8,'LineWidth',5)
geoplot(node_coords(end_idx,1),node_coords(end_idx,2),'g*','MarkerSize',8,'LineWidth',5)

legend('Waypoints','Origin','Destination');
title('RRG with Multiple Alternative Paths');
hold off;

set(gcf,"Position",[0 0 1200 800])
geolimits([13 57],[-130 -60])
end

%% Export TOS
if unique_path_found
    TOS.options = all_paths;
    RTK = all_path_distances / 9.96; %min
    TOS.RTK = RTK - RTK(1);
    % save("TOS","TOS");
else
    TOS = [];
end

end