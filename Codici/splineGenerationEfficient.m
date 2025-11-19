function [splineStruct] = splineGenerationEfficient(points, discretizationLength)


% Spline coefficient retrieval
coeff = spline(points(1,:), points(2,:));
coeff = coeff.coefs;

% Spline intervals
xCoords = points(1,:);

intervalLength = abs(diff(xCoords));
nPoints = round(intervalLength./discretizationLength) ;
 

xQuery = zeros(1,sum(nPoints) - length(nPoints) + 1 ) ;
yQuery = zeros(1,sum(nPoints) - length(nPoints) + 1 ) ;
dyQuery = zeros(1,sum(nPoints) - length(nPoints) + 1 ) ;
ddyQuery = zeros(1,sum(nPoints) - length(nPoints) + 1 ) ;



xLocalQuery = linspace(xCoords(1), xCoords(2), nPoints(1)) ;
xQuery(1:nPoints(1)) = xLocalQuery;

yLocalQuery = coeff(1,:)*[(xLocalQuery-xCoords(1)).^3 ; (xLocalQuery-xCoords(1)).^2 ; (xLocalQuery-xCoords(1)) ; ones(1,length(xLocalQuery))];
yQuery(1:nPoints(1)) = yLocalQuery;

dyLocalQuery = coeff(1,:)*[3*(xLocalQuery-xCoords(1)).^2 ; 2*(xLocalQuery-xCoords(1)) ; ones(1,length(xLocalQuery)) ; zeros(1,length(xLocalQuery))];
dyQuery(1:nPoints(1)) = dyLocalQuery;

ddyLocalQuery = coeff(1,:)*[6*(xLocalQuery-xCoords(1)) ; 2*ones(1,length(xLocalQuery)) ; zeros(1,length(xLocalQuery)) ; zeros(1,length(xLocalQuery))];
ddyQuery(1:nPoints(1)) = ddyLocalQuery;



count = nPoints(1) + 1;

for ii = 2:length(intervalLength)

    idxStart = count ; 
    idxEnd = idxStart + nPoints(ii) - 2;
    count = idxEnd + 1 ;

    xLocalQuery = linspace(xCoords(ii), xCoords(ii+1), nPoints(ii)) ; 
    xQuery(idxStart:idxEnd) = xLocalQuery(2:end);

    yLocalQuery = coeff(ii,:)*[(xLocalQuery(2:end)-xCoords(ii)).^3 ; (xLocalQuery(2:end)-xCoords(ii)).^2 ; (xLocalQuery(2:end)-xCoords(ii)) ; ones(1,length(xLocalQuery)-1)];
    yQuery(idxStart:idxEnd) =  yLocalQuery;

    dyLocalQuery = coeff(ii,:)*[3*(xLocalQuery(2:end)-xCoords(ii)).^2 ; 2*(xLocalQuery(2:end)-xCoords(ii)) ;  ones(1,length(xLocalQuery)-1) ; zeros(1,length(xLocalQuery)-1)];
    dyQuery(idxStart:idxEnd) =  dyLocalQuery;

    ddyLocalQuery = coeff(ii,:)*[6*(xLocalQuery(2:end)-xCoords(ii)) ; 2* ones(1,length(xLocalQuery)-1) ;  zeros(1,length(xLocalQuery)-1) ; zeros(1,length(xLocalQuery)-1)];
    ddyQuery(idxStart:idxEnd) =  ddyLocalQuery;

end

splineStruct.profile(:,1) = xQuery ;
splineStruct.profile(:,2) = yQuery ;

splineStruct.dProfile(:,1) = xQuery ; 
splineStruct.dProfile(:,2) = dyQuery ;

splineStruct.ddProfile(:,1) = xQuery ;
splineStruct.ddProfile(:,2) = ddyQuery ;

end