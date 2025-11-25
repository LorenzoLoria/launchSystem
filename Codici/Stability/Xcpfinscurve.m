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

close all;
clc;

%% --------------------------------------------------------------
% 1. Load and display the image
%% --------------------------------------------------------------
img = imread('Xcpfins.png');
figure;
imshow(img);
title('Click curve points with getpts, then press ENTER');

%% --------------------------------------------------------------
% 2. Digitize curve points manually
%% --------------------------------------------------------------
[x_pix, y_pix] = getpts;     % click points on the curve

%% --------------------------------------------------------------
% 3. Convert pixels → data values
%% --------------------------------------------------------------
[H, W, ~] = size(img);

x_min = 0; x_max = 5;   % Mach limits
y_min = 0.2; y_max = 0.5; % CP limits

X = x_pix / W * (x_max - x_min) + x_min;
Y = (H - y_pix) / H * (y_max - y_min) + y_min;

%% --------------------------------------------------------------
% 4. Select only the region 0.7 <= Mach <= 2
%% --------------------------------------------------------------
mask = (X >= 0.7) & (X <= 2);
X_sub = X(mask);
Y_sub = Y(mask);

%% --------------------------------------------------------------
% 4b. Add exact endpoints by interpolation
%% --------------------------------------------------------------
y_start = interp1(X_sub, Y_sub, 0.7, 'linear', 'extrap');  % CP at Mach 0.7
y_end   = interp1(X_sub, Y_sub, 2, 'linear', 'extrap');    % CP at Mach 2

X_sub = [0.7; X_sub(:); 2];
Y_sub = [y_start; Y_sub(:); y_end];

%% --------------------------------------------------------------
% 5. Fit a quadratic that passes exactly through endpoints + middle point
%% --------------------------------------------------------------
x1 = 0.7;  y1 = y_start;
x2 = 2;    y2 = y_end;
mid_idx = round(length(X_sub)/2);
x3 = X_sub(mid_idx); y3 = Y_sub(mid_idx);

A = [x1^2 x1 1;
     x2^2 x2 1;
     x3^2 x3 1];
b = [y1; y2; y3];

p = A\b;   % coefficients [a; b; c]

%% --------------------------------------------------------------
% 6. Evaluate the quadratic on a fine grid
%% --------------------------------------------------------------
xx = linspace(0.7, 2, 300);
yy = polyval(p, xx);

%% --------------------------------------------------------------
% 7. Plot the results
%% --------------------------------------------------------------
figure;
plot(X_sub, Y_sub, 'o', 'MarkerSize', 8, 'DisplayName', 'Digitized points'); hold on;
plot(xx, yy, 'LineWidth', 2, 'DisplayName', 'Quadratic fit (constrained)');
xlabel('Mach number');
ylabel('X_{cp}/c_{mac}');
title('Quadratic Fit of X_{cp}/c_{mac} for Mach 0.7–2');
legend('Location','best');
grid on;

%% --------------------------------------------------------------
% 8. Display the polynomial
%% --------------------------------------------------------------
disp('Polynomial coefficients p = [a b c] for:  a*M^2 + b*M + c');
disp(p);
