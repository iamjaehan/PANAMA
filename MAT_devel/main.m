% Load data
controlledFlights = load("controlledFlights.mat"); controlledFlights = controlledFlights.controlledFlights;
sectors = load("sectors.mat"); sectors = sectors.sectors;
flight_paths = load("flight_paths.mat"); flight_paths = flight_paths.flight_paths;
flight_sector_map = load("flight_sector_map.mat"); assigned_sector = flight_sector_map.flight_sector_map;
flightn = length(controlledFlights);

% CTOP initiated when 
%% CTOP


%% TACo