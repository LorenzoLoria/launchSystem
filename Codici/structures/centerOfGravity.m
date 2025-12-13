function cg = centerOfGravity(massVec, positionVec)

% Calculate the sum of moments (m_i * x_i)
sumOfMoments = sum(massVec .* positionVec);

% Calculate the total mass (sum of m_i)
totalMass = sum(massVec);

% Calculate the center of gravity

cg = sumOfMoments / totalMass;


end