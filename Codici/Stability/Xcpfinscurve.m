close all; clear;

%% ============================================================
% 1) Load image
%% ============================================================
img = imread('Xcpfins.png');
imgH = size(img,1);
imgW = size(img,2);

figure; imshow(img); title('Click Mach=0 tick, then Mach=5 tick');
[x0, y0] = ginput(1);
[x5, y5] = ginput(1);

figure; imshow(img); title('Click Y=0.25 tick, then Y=0.50 tick');
[xY1, yY1] = ginput(1);
[xY2, yY2] = ginput(1);

%% ============================================================
% 2) Horizontal axis calibration (Mach)
%% ============================================================
M0  = 0;
M5  = 5;

aX = (M5 - M0) / (x5 - x0);
bX = M0 - aX * x0;

% Forward: Mach → pixel-x
pixX = @(M) (M - bX) / aX;


%% ============================================================
% 3) Vertical axis calibration (Ycp)
%% ============================================================
Y1 = 0.25;
Y2 = 0.50;

aY = (Y2 - Y1) / (yY2 - yY1);
bY = Y1 - aY * yY1;

% Determine if the axis is inverted
% If the user clicked Y=0.25 ABOVE Y=0.50, screen coordinates decrease upward
if yY1 < yY2
    % Y increases DOWNWARD → must invert
    pixY = @(Y) imgH - ((Y - bY) / aY);
else
    % Y increases UPWARD → do NOT invert
    pixY = @(Y) ((Y - bY) / aY);
end


%% ============================================================
% 4) Print calibration to verify
%% ============================================================
fprintf('\n=== Calibration check ===\n');
fprintf('Mach at clicked points:\n');
fprintf('  Mach(x0) = %.3f (expected 0)\n', aX*x0 + bX);
fprintf('  Mach(x5) = %.3f (expected 5)\n\n', aX*x5 + bX);

fprintf('Ycp at clicked points:\n');
fprintf('  Y(xY1) = %.3f (expected 0.25)\n', aY*yY1 + bY);
fprintf('  Y(xY2) = %.3f (expected 0.50)\n\n', aY*yY2 + bY);

%% ============================================================
% 5) Define curves
%% ============================================================
Xcp_supersonic = @(M,A) (A*sqrt(M.^2-1) - 0.67) ./ (2*A*sqrt(M.^2-1) - 1);

A_values = [1,2,3];
colors = [0 0.447 0.741;
          0.85 0.325 0.098;
          0.929 0.694 0.125];

M_low  = linspace(0,0.7,80);
M_mid  = linspace(0.7,2,120);
M_high = linspace(2,5,200);

%% ============================================================
% 6) Draw curves
%% ============================================================
figure; imshow(img); hold on;
title('Correctly Overlaid Curves');

for k = 1:length(A_values)

    A = A_values(k);
    col = colors(k,:);

    X2 = Xcp_supersonic(2,A);

    % Parabolic mid-section fit
    A_mat = [
        0.7^2,  0.7,  1;
        2^2,    2,    1;
        1.35^2, 1.35, 1 ];
    b_vec = [0.25; X2; 0.40];

    abc = A_mat \ b_vec;
    a = abc(1); b = abc(2); c = abc(3);

    X_low  = 0.25*ones(size(M_low));
    X_mid  = a*M_mid.^2 + b*M_mid + c;
    X_high = Xcp_supersonic(M_high, A);

    % Convert to pixels
    pxL = pixX(M_low);  pyL = pixY(X_low);
    pxM = pixX(M_mid);  pyM = pixY(X_mid);
    pxH = pixX(M_high); pyH = pixY(X_high);

    % Keep only points inside the image
    inside = @(x,y) x>=1 & x<=imgW & y>=1 & y<=imgH;

    valid = inside(pxL,pyL); plot(pxL(valid), pyL(valid), 'Color', col, 'LineWidth', 2);
    valid = inside(pxM,pyM); plot(pxM(valid), pyM(valid), 'Color', col, 'LineWidth', 2);
    valid = inside(pxH,pyH); plot(pxH(valid), pyH(valid), 'Color', col, 'LineWidth', 2);
end

hold off;
