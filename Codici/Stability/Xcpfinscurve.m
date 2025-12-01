% close all;
% clc;
% 
% %% --------------------------------------------------------------
% % 1. Load and display the image
% %% --------------------------------------------------------------
% img = imread('Xcpfins.png');
% figure;
% imshow(img);
% title('Click curve points with getpts, then press ENTER');
% 
% %% --------------------------------------------------------------
% % 2. Digitize curve points manually
% %% --------------------------------------------------------------
% [x_pix, y_pix] = getpts;     % click points on the curve
% 
% %% --------------------------------------------------------------
% % 3. Convert pixels → data values
% % NOTE: adjust axis limits to match your plot
% %% --------------------------------------------------------------
% 
% % Get image size
% [H, W, ~] = size(img);
% 
% % True axis limits in your plot
% x_min = 0;      % Mach = 0
% x_max = 5;      % Mach = 5
% y_min = 0.2;    % bottom of the CP axis
% y_max = 0.5;    % top of the CP axis
% 
% % Convert X pixels → Mach numbers
% X = x_pix / W * (x_max - x_min) + x_min;
% 
% % Convert Y pixels → CP values (invert Y-axis)
% Y = (H - y_pix) / H * (y_max - y_min) + y_min;
% 
% %% --------------------------------------------------------------
% % 4. Select only the region 0.7 <= Mach <= 2
% %% --------------------------------------------------------------
% mask = (X >= 0.7) & (X <= 2);
% X_sub = X(mask);
% Y_sub = Y(mask);
% 
% %% --------------------------------------------------------------
% % 5. Fit a polynomial (order 2 recommended)
% %% --------------------------------------------------------------
% p = polyfit(X_sub, Y_sub, 2);
% 
% %% --------------------------------------------------------------
% % 6. Evaluate the polynomial on a fine grid
% %% --------------------------------------------------------------
% xx = linspace(0.7, 2, 300);
% yy = polyval(p, xx);
% 
% %% --------------------------------------------------------------
% % 7. Plot the results
% %% --------------------------------------------------------------
% figure;
% plot(X_sub, Y_sub, 'o', 'MarkerSize', 8, 'DisplayName', 'Digitized points'); hold on;
% plot(xx, yy, 'LineWidth', 2, 'DisplayName', 'Polynomial fit (0.7–2)');
% xlabel('Mach number');
% ylabel('X_{cp}/c_{mac}');
% title('Fit of X_{cp}/c_{mac} for Mach 0.7–2');
% legend('Location','best');
% grid on;
% 
% %% --------------------------------------------------------------
% % 8. Display the polynomial
% %% --------------------------------------------------------------
% disp('Polynomial coefficients p = [a b c] for:  a*M^2 + b*M + c');
% disp(p);

% close all;
% clc;
% 
% %% --------------------------------------------------------------
% % 1. Load and display the image
% %% --------------------------------------------------------------
% img = imread('Xcpfins.png');
% figure;
% imshow(img);
% title('Click curve points with getpts, then press ENTER');
% 
% %% --------------------------------------------------------------
% % 2. Digitize curve points manually
% %% --------------------------------------------------------------
% [x_pix, y_pix] = getpts;     % click points on the curve
% 
% %% --------------------------------------------------------------
% % 3. Convert pixels → data values
% %% --------------------------------------------------------------
% [H, W, ~] = size(img);
% 
% x_min = 0; x_max = 5;   % Mach limits
% y_min = 0.2; y_max = 0.5; % CP limits
% 
% X = x_pix / W * (x_max - x_min) + x_min;
% Y = (H - y_pix) / H * (y_max - y_min) + y_min;
% 
% %% --------------------------------------------------------------
% % 4. Select only the region 0.7 <= Mach <= 2
% %% --------------------------------------------------------------
% mask = (X >= 0.7) & (X <= 2);
% X_sub = X(mask);
% Y_sub = Y(mask);
% 
% %% --------------------------------------------------------------
% % 4b. Add exact endpoints by interpolation
% %% --------------------------------------------------------------
% y_start = interp1(X_sub, Y_sub, 0.7, 'linear', 'extrap');  % CP at Mach 0.7
% y_end   = interp1(X_sub, Y_sub, 2, 'linear', 'extrap');    % CP at Mach 2
% 
% X_sub = [0.7; X_sub(:); 2];
% Y_sub = [y_start; Y_sub(:); y_end];
% 
% %% --------------------------------------------------------------
% % 5. Fit a quadratic that passes exactly through endpoints + middle point
% %% --------------------------------------------------------------
% x1 = 0.7;  y1 = y_start;
% x2 = 2;    y2 = y_end;
% mid_idx = round(length(X_sub)/2);
% x3 = X_sub(mid_idx); y3 = Y_sub(mid_idx);
% 
% A = [x1^2 x1 1;
%      x2^2 x2 1;
%      x3^2 x3 1];
% b = [y1; y2; y3];
% 
% p = A\b;   % coefficients [a; b; c]
% 
% %% --------------------------------------------------------------
% % 6. Evaluate the quadratic on a fine grid
% %% --------------------------------------------------------------
% xx = linspace(0.7, 2, 300);
% yy = polyval(p, xx);
% 
% %% --------------------------------------------------------------
% % 7. Plot the results
% %% --------------------------------------------------------------
% figure;
% plot(X_sub, Y_sub, 'o', 'MarkerSize', 8, 'DisplayName', 'Digitized points'); hold on;
% plot(xx, yy, 'LineWidth', 2, 'DisplayName', 'Quadratic fit (constrained)');
% xlabel('Mach number');
% ylabel('X_{cp}/c_{mac}');
% title('Quadratic Fit of X_{cp}/c_{mac} for Mach 0.7–2');
% legend('Location','best');
% grid on;
% 
% %% --------------------------------------------------------------
% % 8. Display the polynomial
% %% --------------------------------------------------------------
% disp('Polynomial coefficients p = [a b c] for:  a*M^2 + b*M + c');
% disp(p);

function [p, fit_mid_fun, Mfull, Xcp] = Xcpfinscurve(imgFilename)
% Xcpfinscurve  Fit Xcp/cmac curve for A = 1 using a PARABOLIC fit
%
% Usage:
%   [p, fit_mid_fun, Mfull, Xcp] = Xcpfinscurve('Xcpfins.png');
%
% Interaction:
%   1) Click tick mark at Y = 0.5 (left axis)
%   2) Click tick mark at Y = 0.2
%   3) Click point at Mach = 0.7 on the curve
%   4) Click point at Mach = 2.0 on the curve
%   5) Click curve points ONLY between Mach 0.7 and 2.0
%      Press ENTER to finish

    if nargin < 1
        imgFilename = 'Xcpfins.png';
    end

%% ========================================================================
% 1) Load image & calibrate Y axis
%% ========================================================================
    img = imread(imgFilename);
    figure; imshow(img); title('Click Y=0.5 tick, then Y=0.2 tick');

    % Click ticks for Y = 0.5 and Y = 0.2
    [~, y05_pix] = ginput(1);   % top tick (0.5)
    [~, y02_pix] = ginput(1);   % bottom tick (0.2)

    y_min = 0.2;
    y_max = 0.5;

    % Linear mapping Y = aY * pixel + bY
    aY = (y_max - y_min) / (y02_pix - y05_pix);
    bY = y_max - aY * y05_pix;

%% ========================================================================
% 2) Mach axis calibration
%% ========================================================================
    figure; imshow(img);
    title('Click Mach=0.7 point, then Mach=2.0 point, then digitize curve');

    [x1_pix, y1_pix] = ginput(1);  % Mach = 0.7
    [x2_pix, y2_pix] = ginput(1);  % Mach = 2.0

    M1 = 0.7;
    M2 = 2.0;

    % Mach mapping M = aX * pixel + bX
    aX = (M2 - M1) / (x2_pix - x1_pix);
    bX = M1 - aX * x1_pix;

    % Known CP values at these Mach numbers
    X1 = aY * y1_pix + bY;
    X2 = aY * y2_pix + bY;

%% ========================================================================
% 3) Digitize curve points between Mach 0.7 and 2
%% ========================================================================
    [x_pix, y_pix] = ginput;  % digitized curve between 0.7–2

    M_points = aX * x_pix + bX;
    X_points = aY * y_pix + bY;

%% ========================================================================
% 4) Parabolic fit with derivative constraint at M2
%% ========================================================================
    % f(M) = a M^2 + b M + c
    % derivative constraint: f'(M2) = 0 → b = -2 a M2

    M_mat = [M_points(:).^2 - 2*M2*M_points(:), ones(length(M_points),1)];
    coeff = M_mat \ X_points(:);

    a = coeff(1);
    b = -2*a*M2;
    c = coeff(2);

    p = [a b c];
    fit_mid_fun = @(M) polyval(p, M);

%% ========================================================================
% 5) Assemble full piecewise curve from Mach 0 to 5
%% ========================================================================
    A = 1;
    X_supersonic = @(M) (A*sqrt(M.^2 - 1) - 0.67) ./ (2*A*sqrt(M.^2 - 1) - 1);

    Mfull = linspace(0, 5, 400);
    Xcp = zeros(size(Mfull));

    for k = 1:length(Mfull)
        Mk = Mfull(k);
        if Mk < M1
            Xcp(k) = X1;
        elseif Mk <= M2
            Xcp(k) = fit_mid_fun(Mk);
        else
            Xcp(k) = X_supersonic(Mk);
        end
    end

%% ========================================================================
% 6) Overlay final curve on image
%% ========================================================================
    % Convert curve to pixel coordinates
    x_full_pix = (Mfull - bX) / aX;
    y_full_pix = (Xcp   - bY) / aY;

    % Original digitized points for checking
    x_pix_points = x_pix;
    y_pix_points = y_pix;

    figure; imshow(img); hold on; axis ij;
    plot(x_full_pix, y_full_pix, 'r-', 'LineWidth', 2);
    plot(x_pix_points, y_pix_points, 'bo', 'MarkerFaceColor', 'b');

    title('Correct Overlay — Parabolic Fit');
    legend('Fitted curve', 'Digitized points');

end
