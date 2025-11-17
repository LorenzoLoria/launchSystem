function [out] = splineGeneration(varargin)
% Spline matlab function works with not-a-knot type splines   
% Primo input: punti per spline
% Secondo input (optional): query points in cui restituire la spline
%
%
% Primi 3 output: spline completa con derivate 1° e 2° 
% 4° output (optional): valutazione spline nei query points
%
%

if nargin == 1

    points = cell2mat(varargin(1));
    flag = 0;

elseif nargin == 2

    points = cell2mat(varargin(1));
    queryPoints = cell2mat(varargin(2));
    flag = 1;

else

    error("Too many inputs")

end


% Spline coefficient retrieval
coeff = spline(points(:,1), points(:,2));
coeff = coeff.coefs;

% Spline intervals
xCoords = points(:,1)';
intervals = zeros(length(xCoords)-1,2);
for ii = 1:length(xCoords)-1
    intervals(ii,:) = xCoords(ii:ii+1);
end

% Unique function for the spline (all polynomials defined in the single interval are united)
correctInterval = @(x) find(([(x >= intervals(1,1)) .* (x <= intervals (1,2)) ; (x > intervals(2:end,1)) .* (x <= intervals (2:end,2))])~=0);
polynomial = @(x) sum( [coeff(correctInterval(x),:) * [(x-intervals(correctInterval(x),1)).^3 ; (x-intervals(correctInterval(x),1)).^2 ; (x-intervals(correctInterval(x),1)) ; 1] 0]) ;

% Unique functions for 1° and 2° spline derivatives 
firstDer = @(x) sum( [ coeff(correctInterval(x),1:end-1) * [3*(x-intervals(correctInterval(x),1)).^2 ; 2*(x-intervals(correctInterval(x),1)) ; 1] 0 ]);
secondDer = @(x) sum( [ coeff(correctInterval(x),1:end-2) * [6*(x-intervals(correctInterval(x),1)) ; 2 ] 0 ]);

if flag

    % Query points retrieval
    profile(:,1) = queryPoints ;
    dProfile(:,1)= queryPoints ;
    ddProfile(:,1) = queryPoints;

    for ii = 1:length(queryPoints)
        profile(ii,2) = polynomial(queryPoints(ii)) ;
        dProfile(ii,2) = firstDer(queryPoints(ii)) ;
        ddProfile(ii,2) = secondDer(queryPoints(ii));
    end

end


if nargin == 1

    out.polynomial = polynomial;
    out.firstDer = firstDer;
    out.secondDer = secondDer;



elseif nargin == 2

    out.polynomial =  polynomial;
    out.firstDer =  firstDer;
    out.secondDer =  secondDer;
    out.profile = profile;
    out.dProfile = dProfile;
    out.ddProfile = ddProfile;

else 

    error("Wrong number of inputs")

end


end