data = load('FIR_data.mat');  % Load FIR data
fir_names = fieldnames(data);  % Get all FIR names

% Loop through each FIR and convert to numeric matrix
for i = 1:length(fir_names)
    field_name = fir_names{i};
    if iscell(data.(field_name))  % Check if stored as cell
        data.(field_name) = cell2mat(data.(field_name));  % Convert to numeric matrix
    end
end

save("FIR_coord.mat","data");