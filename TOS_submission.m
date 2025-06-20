n = 30; % Number of flights
airlineList = ["AA","UA","OZ","JB","SW"];

TOS_set = cell(n,1);

for i = 1:n
    % origin = randi(5000);
    % destination = randi(5000);
    % while origin == destination
    %     destination = randi(5000);
    % end
    origin = randi(5000);
    while ~((G.Nodes(origin,:).Lon < -93 && G.Nodes(origin,:).Lon > -100) && (G.Nodes(origin,:).Lat < 40 && G.Nodes(origin,:).Lat > 28))
        origin = randi(5000);
    end
    destination = randi(5000);
    while ~((G.Nodes(destination,:).Lon < -72 && G.Nodes(destination,:).Lon > -80) && (G.Nodes(destination,:).Lat < 44 && G.Nodes(destination,:).Lat > 36)) || destination == origin
        destination = randi(5000);
    end

    while true
        TOS = generateTOS(origin,destination,false);
        if ~isempty(TOS)
            break;
        else
            while origin == destination
                destination = randi(5000);
            end
        end
    end

    % Assign fake flight information
    flightNum = randsample(airlineList,1) + num2str(randsample(0:9,1)) + num2str(randsample(0:9,1)) ...
        + num2str(randsample(0:9,1)) + num2str(randsample(0:9,1));
    depTime = [2025 2 24 randsample(8:9,1) randsample(0:59,1)];
    altitude = randsample(250:10:360,1);

    % Add data
    TOS.flightNum = char(flightNum);
    TOS.depTime = depTime;
    TOS.altitude = altitude;
    TOS_set{i} = TOS;
end

save("TOS.mat","TOS_set");
