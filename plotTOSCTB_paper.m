tos = load("TOS.mat"); tos = tos.TOS_set;
G = load("US_waypoint_graph.mat");
G = G.G;
node_coords = G.Nodes{:,2:3};

figure(3);
clf;

%plot waypoints
geoplot(node_coords(:,1), node_coords(:,2), '.', 'MarkerSize', 3, 'Color', [0 0 0.9 0.1]); % Graph nodes
hold on

% plot ARTCC
fir = load("FIR_coord.mat");
fir = fir.data;
fir_names = fieldnames(fir);
fir_num = length(fir_names);
for i = 1:fir_num
    local_coord = fir.(fir_names{i});
    if i == 2 || i == 8 || i == 20 %12 3 13
        % geoplot(local_coord(:,2),local_coord(:,1),'.-','LineWidth',1,'Color',[1 0 0 0.9])
    else
        geoplot(local_coord(:,2),local_coord(:,1),'-','LineWidth',1,'Color',[0 0 0 1])
    end
end

for i = 1:fir_num
    local_coord = fir.(fir_names{i});
    if i == 2 || i == 8 || i == 20 %12 3 13
        geoplot(local_coord(:,2),local_coord(:,1),'-','LineWidth',3,'Color',[1 0 0 1])
    end
end

% Post processing

% For overall plot
title('ARTCC Boundary and Waypoint Map');
hold off;
set(gcf,"Position",[30 100 1000 700])
set(gca,"FontName","times","FontSize",23)
geolimits([13 57],[-130 -60])
exportgraphics(gca,'../Plots/giniPerk.pdf', 'Resolution',300);

% For tos plot