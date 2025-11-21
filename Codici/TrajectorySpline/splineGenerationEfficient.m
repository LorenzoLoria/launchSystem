function [splineStruct] = splineGenerationEfficient(points, nPointsPerInterval)


% Spline coefficient retrieval
coeff = spline(points(1,:), points(2,:));
coeff = coeff.coefs;

% Spline intervals
xCoords = points(1,:); 

xQuery = zeros(1,(length(points(1,:))-1)*(nPointsPerInterval-1) + 1 ) ;
yQuery = zeros(1,(length(points(1,:))-1)*(nPointsPerInterval-1) + 1 ) ;
dyQuery = zeros(1,(length(points(1,:))-1)*(nPointsPerInterval-1) + 1 ) ;
ddyQuery = zeros(1,(length(points(1,:))-1)*(nPointsPerInterval-1) + 1 ) ;


xLocalQuery = linspace(xCoords(1), xCoords(2), nPointsPerInterval) ;
xQuery(1:nPointsPerInterval) = xLocalQuery;

yLocalQuery = coeff(1,:)*[(xLocalQuery-xCoords(1)).^3 ; (xLocalQuery-xCoords(1)).^2 ; (xLocalQuery-xCoords(1)) ; ones(1,length(xLocalQuery))];
yQuery(1:nPointsPerInterval) = yLocalQuery;

dyLocalQuery = coeff(1,:)*[3*(xLocalQuery-xCoords(1)).^2 ; 2*(xLocalQuery-xCoords(1)) ; ones(1,length(xLocalQuery)) ; zeros(1,length(xLocalQuery))];
dyQuery(1:nPointsPerInterval) = dyLocalQuery;

ddyLocalQuery = coeff(1,:)*[6*(xLocalQuery-xCoords(1)) ; 2*ones(1,length(xLocalQuery)) ; zeros(1,length(xLocalQuery)) ; zeros(1,length(xLocalQuery))];
ddyQuery(1:nPointsPerInterval) = ddyLocalQuery;



for ii = 2:length(points(1,:))-1

    xLocalQuery = linspace(xCoords(ii), xCoords(ii+1), nPointsPerInterval) ; 
    xQuery(1+nPointsPerInterval+(nPointsPerInterval-1)*(ii-2):nPointsPerInterval+(nPointsPerInterval-1)*(ii-1)) = xLocalQuery(2:end);

    yLocalQuery = coeff(ii,:)*[(xLocalQuery(2:end)-xCoords(ii)).^3 ; (xLocalQuery(2:end)-xCoords(ii)).^2 ; (xLocalQuery(2:end)-xCoords(ii)) ; ones(1,length(xLocalQuery)-1)];
    yQuery(1+nPointsPerInterval+(nPointsPerInterval-1)*(ii-2):nPointsPerInterval+(nPointsPerInterval-1)*(ii-1)) =  yLocalQuery;

    dyLocalQuery = coeff(ii,:)*[3*(xLocalQuery(2:end)-xCoords(ii)).^2 ; 2*(xLocalQuery(2:end)-xCoords(ii)) ;  ones(1,length(xLocalQuery)-1) ; zeros(1,length(xLocalQuery)-1)];
    dyQuery(1+nPointsPerInterval+(nPointsPerInterval-1)*(ii-2):nPointsPerInterval+(nPointsPerInterval-1)*(ii-1)) =  dyLocalQuery;

    ddyLocalQuery = coeff(ii,:)*[6*(xLocalQuery(2:end)-xCoords(ii)) ; 2* ones(1,length(xLocalQuery)-1) ;  zeros(1,length(xLocalQuery)-1) ; zeros(1,length(xLocalQuery)-1)];
    ddyQuery(1+nPointsPerInterval+(nPointsPerInterval-1)*(ii-2):nPointsPerInterval+(nPointsPerInterval-1)*(ii-1)) =  ddyLocalQuery;

end

splineStruct.profile(:,1) = xQuery ;
splineStruct.profile(:,2) = yQuery ;

splineStruct.dProfile(:,1) = xQuery ; 
splineStruct.dProfile(:,2) = dyQuery ;

splineStruct.ddProfile(:,1) = xQuery ;
splineStruct.ddProfile(:,2) = ddyQuery ;

end