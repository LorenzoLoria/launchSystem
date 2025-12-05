function v_wind = wind_vs_altitude(h)
%
% INPUT
% - h       : altitude [km]
% 
% OUTPUT
% - v_wind  : wind speed [m/s]

if h >= 0 && h <= 9.6
    v_wind = (6.9288 * h + 9.144);  
elseif h > 9.6 && h <= 14
    v_wind = 76.2;
elseif h > 14 && h <= 20
    v_wind = (76.2 - 8.9474 * (h - 14));
elseif h > 20
    v_wind = 24.384;

end