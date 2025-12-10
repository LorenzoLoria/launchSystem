function [v_wind_annual_vec, variance_vec] = GRAM07_HWM07_annual(h_vec)
%
% INPUT
% - h_vec            : altitude 0–100 km (vector)
%
% OUTPUT
% - v_wind_annual_vec: annual wind profile (vector, same size as h_vec)
% - variance_vec     : variance of seasonal wind at each altitude (same size as h_vec)


% Profili stagionali
for i = 1 : length(h_vec)
    h = h_vec(i);
    v_wind_autumn_vec(i) = GRAM07_HWM07(h, 1);
    v_wind_winter_vec(i) = GRAM07_HWM07(h, 2);
    v_wind_spring_vec(i) = GRAM07_HWM07(h, 3);
    v_wind_summer_vec(i) = GRAM07_HWM07(h, 4);
end


% Media annuale (semplice media aritmetica delle 4 stagioni)
v_wind_annual_vec = 1/4 * ( ...
      v_wind_autumn_vec ...
    + v_wind_winter_vec ...
    + v_wind_spring_vec ...
    + v_wind_summer_vec );

% Scostamenti dalla media
diff_autumn = v_wind_autumn_vec - v_wind_annual_vec;
diff_winter = v_wind_winter_vec - v_wind_annual_vec;
diff_spring = v_wind_spring_vec - v_wind_annual_vec;
diff_summer = v_wind_summer_vec - v_wind_annual_vec;

% Varianza di popolazione sulle 4 stagioni:
% sigma^2(h) = (1/4) * Σ_s [ V_s(h) - mu(h) ]^2
variance_vec = 1/4 * ( ...
      diff_autumn.^2 ...
    + diff_winter.^2 ...
    + diff_spring.^2 ...
    + diff_summer.^2 );

end
