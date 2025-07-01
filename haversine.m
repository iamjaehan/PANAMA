function out = haversine(lat1,lon1,lat2,lon2)

R = 6371; %Earth radius in km
d2r = pi/180;
dlat = d2r*(lat2 - lat1);
dlon = d2r*(lon2 - lon1);
a = sin(dlat / 2)^2 + cos(lat1*d2r) * cos(lat2*d2r) * sin(dlon / 2)^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));
out = R * c;

end