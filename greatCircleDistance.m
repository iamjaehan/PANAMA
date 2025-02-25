function d = greatCircleDistance(origin,destination, r)
d2r = pi/180;
phi_s = origin(:,1) * d2r;
lambda_s = origin(:,2) * d2r;
phi_f = destination(1) * d2r;
lambda_f = destination(2) * d2r;

% If no radius supplied, assume the mean radius of the earth in km
if nargin < 3
    r = 6371.01; % km
end

% Compute Delta lambda (delta longitude)
Delta_lambda = lambda_f - lambda_s;

% Compute Delta sigma (central angle)
Delta_sigma = atan2(sqrt((cos(phi_f).*sin(Delta_lambda)).^2 + (cos(phi_s).*sin(phi_f) - sin(phi_s).*cos(phi_f).*cos(Delta_lambda)).^2), ...
    sin(phi_s).*sin(phi_f) + cos(phi_s).*cos(phi_f).*cos(Delta_lambda));

d = r*Delta_sigma;