n = 100; % Number of flights
d = 5; % min - time resolution
fStep = 500; % final step length
tf = fStep * d; % min - final time
airlineList = ["AA","UA","OZ","JB","SW"];
v = 900; %km/h

TOS_set = cell(n,1);
G = load("US_waypoint_graph.mat");
G = G.G;

fir = load("FIR_coord.mat");
fir = fir.data;
fir_names = fieldnames(fir);
fir_num = length(fir_names);

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
    depTime = datetime(depTime(1), depTime(2), depTime(3), ...
               depTime(4), depTime(5), 0);
    depTime = posixtime(depTime);
    altitude = randsample(250:10:360,1);

    % Add data
    TOS.flightNum = char(flightNum);
    TOS.depTime = depTime;
    TOS.altitude = altitude;
    TOS_set{i} = TOS;
end

for i = 1:n
    TOS = TOS_set{i};
    dep = TOS.depTime;
    options = TOS.options;
    firTime = cell(length(options), 1);

    for j = 1:length(options)
        wpts = options{j};
        firLog = {}; % 결과 저장용
        
        curTime = dep; % 출발시간 offset (s 단위)
        prevFIR = "";
        inTime = NaN;

        for k = 1:length(wpts)-1
            lat1 = G.Nodes.Lat(wpts(k));
            lon1 = G.Nodes.Lon(wpts(k));
            lat2 = G.Nodes.Lat(wpts(k+1));
            lon2 = G.Nodes.Lon(wpts(k+1));
            dist = haversine(lat1, lon1, lat2, lon2);
            tseg = (dist / v) * 3600; % seconds

            % 샘플링 (예: 10개 점 생성)
            Nsample = 10;
            for s = 0:Nsample
                frac = s / Nsample;
                lat = lat1 + frac * (lat2 - lat1);
                lon = lon1 + frac * (lon2 - lon1);
                tNow = curTime + frac * tseg;

                % FIR 확인
                currentFIR = "";
                for f = 1:fir_num
                    fname = fir_names{f};
                    poly = fir.(fname);
                    [in,on] = inpolygon(lon, lat, poly(:,1), poly(:,2));
                    if in || on
                        currentFIR = fname;
                        break;
                    end
                end

                % FIR 변하면 기록
                if ~strcmp(currentFIR, prevFIR)
                    if ~(prevFIR == "")
                        firLog = [firLog; {prevFIR, inTime, tNow}];
                    end
                    prevFIR = currentFIR;
                    inTime = tNow;
                end
            end
            curTime = curTime + tseg;
        end

        % 마지막 FIR 퇴장
        if ~(prevFIR == "")
            firLog = [firLog; {prevFIR, inTime, curTime}];
        end

        firTime{j} = firLog;
    end

    TOS.firTime = firTime;
    TOS_set{i} = TOS;
end

save("TOS.mat","TOS_set");
