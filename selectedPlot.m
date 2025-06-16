selected = load("[0]selected.mat");
tos = load("TOS.mat");

selected = selected.selected;
tos = tos.TOS_set;

G = load("US_waypoint_graph.mat");
G = G.G;
node_coords = G.Nodes{:,2:3};

figure(3);
clf;
% geoplot(wplat, wplon, 'm.','MarkerSize',2); % Plot all waypoints
geoplot(node_coords(:,1), node_coords(:,2), 'b.', 'MarkerSize', 3); % Graph nodes
hold on

% Plot multiple paths with different colors
selectIdx = 1;
selection = selected{selectIdx};
for k = 1:10
    localSelection = selection(k);
    localRoute = tos{k}.options{localSelection};
    path_lats = G.Nodes.Lat(localRoute);
    path_lons = G.Nodes.Lon(localRoute);
    % geoplot(path_lats, path_lons, colors{k}, 'LineWidth', 2);
    geoplot(path_lats, path_lons, '-', 'LineWidth', 4, 'Color', [0.1 0.1 0.1 0.4]);
end

% Draw FIR
fir = load("FIR_coord.mat");
fir = fir.data;
fir_names = fieldnames(fir);
fir_num = length(fir_names);
for i = 1:fir_num
    local_coord = fir.(fir_names{i});
    if i == 2 || i == 8 || i == 20
        geoplot(local_coord(:,2),local_coord(:,1),':','LineWidth',2,'Color',[1 0 0 0.9])
    else
        geoplot(local_coord(:,2),local_coord(:,1),':','LineWidth',2,'Color',[0 0 0 0.5])
    end
end

% legend('Waypoints','Origin','Destination');
title('TOS');
hold off;

set(gcf,"Position",[30 100 1200 800])
set(gca,"FontName","times","FontSize",20)
geolimits([13 57],[-130 -60])

saveas(gcf, "ctb"+num2str(selectIdx)+".png")