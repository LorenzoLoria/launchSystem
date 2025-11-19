clear
clc
close all

nOptPoints = 5 ;
rEarth = 6371 * 1e3 ;
x0 = [sqrt(2) * rEarth , 0, -45] ; 
xf = [sqrt(2) * rEarth , 180, 225] ;
target = xf ; 

initialGuessR = linspace(x0(1), x0(1) + 100 * 1e3, ceil(nOptPoints/2) + 1) ; 
initialGuessR = [initialGuessR, linspace( x0(1) + 100 * 1e3, x0(1), floor(nOptPoints/2) + 1) ]; 
initialGuessR = initialGuessY(2:end-1) ; 

initialGuessTheta = linspace(x0(2), xf(2), nOptPoints+2) ; 
initialGuessTheta = initialGuessTheta(2:end-1) ; 
initialGuessV = linspace(1, 1e2, nOptPoints) ; 
initialGuess = [initialGuessR, initialGuessTheta , initialGuessV] ; 
