data = readtable("Waypoint list.xlsx");
datalen = size(data,1);

%% Preprocessing
% dms to radians
for i = 1:datalen
    d = sscanf(data{i,2}{1}, '%d° %d'' %f"');
    d = dms2degrees(d');
    if contains(data{i,2}{1}, 'S')
        data{i,2}{1} = -d;
    else
        data{i,2}{1} = d;
    end
    d = sscanf(data{i,3}{1}, '%d° %d'' %f"');
    d = dms2degrees(d');
    if contains(data{i,3}{1}, 'W')
        data{i,3}{1} = -d;
    else
        data{i,3}{1} = d;
    end
end

% Cell to Mat
wplat = cell2mat(data{:,2});
wplon = cell2mat(data{:,3});

%% Draw US map and waypoints
load usapolygon.mat
figure(1)
clf
geoplot(uslat,uslon,"--")
hold on
geoscatter(wplat,wplon,[],"m")

%% RRG

% Define origin and destination
origin = [wplat(1), wplon(1)];  % Modify as needed
destination = [wplat(end), wplon(end)];

% Parameters for RRG
num_samples = datalen;   % Number of waypoints to sample for the graph
radius = 1.5;        % Connection radius (degrees)

% Construct RRG
G = graph(); % Initialize an empty graph
node_coords = [origin]; % Start with the origin node
G = addnode(G, table(1, origin(1), origin(2), 'VariableNames', {'ID', 'Lat', 'Lon'}));

% Sample waypoints to build the graph
for i = 1:num_samples
    idx = randi(datalen);  % Randomly pick a waypoint
    idx = i;
    new_wp = [wplat(idx), wplon(idx)];

    % Find nearby waypoints within a given radius
    % distances = vecnorm(node_coords - new_wp, 2, 2); % Euclidean distance
    distances = greatCircleDistance(node_coords, new_wp); % Great circle distance
    nearby_nodes = find(distances < radius * pi / 180 * 6371);

    % Add new waypoint to graph
    G = addnode(G, table(i+1, new_wp(1), new_wp(2), 'VariableNames', {'ID', 'Lat', 'Lon'}));
    node_coords = [node_coords; new_wp]; % Append to list

    % Connect to nearby nodes
    for j = 1:length(nearby_nodes)
        G = addedge(G, nearby_nodes(j), i+1, distances(nearby_nodes(j)));
    end
end

% Add the destination
G = addnode(G, table(num_samples+2, destination(1), destination(2), 'VariableNames', {'ID', 'Lat', 'Lon'}));
distances = vecnorm(node_coords - destination, 2, 2);
nearby_nodes = find(distances < radius);
for j = 1:length(nearby_nodes)
    G = addedge(G, nearby_nodes(j), num_samples+2, distances(nearby_nodes(j)));
end



%% Export TOS
TOS = all_paths;
RTK = all_path_distances;
RTK = RTK - all_path_distances(1);

TOS_len = length(all_paths);

for i = 1:TOS_len
    trajIdx = TOS{i};
    traj_len = length(trajIdx);

end
node_coords(start_idx,1);