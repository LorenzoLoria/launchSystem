
function [soundspeed] = soundSpeedFun(h)

[Temperature,Dens] =  atmosnrlmsise00(h,0,0,2025,1,1);

Temp   = Temperature(2);
n.HE = Dens(1);
n.O  = Dens(2);
n.N2 = Dens(3);
n.O2 = Dens(4);
n.Ar = Dens(5);
n.H  = Dens(7);
n.N  = Dens(8);

nAvo = 6.022e23;
mol.HE = n.HE / nAvo;
mol.O = n.O / nAvo;
mol.N2 = n.N2 / nAvo;
mol.O2 = n.O2 / nAvo;
mol.Ar = n.Ar / nAvo;
mol.H = n.H / nAvo;
mol.N = n.N / nAvo;

massMol.HE = 4;
massMol.O  = 16;
massMol.N2 = 28;
massMol.O2 = 32;
massMol.Ar = 40;
massMol.H  = 1;
massMol.N  = 14;

mass.HE = mol.HE * massMol.HE;
mass.O = mol.O * massMol.O;
mass.N2 = mol.N2 * massMol.N2;
mass.O2 = mol.O2 * massMol.O2;
mass.Ar = mol.Ar * massMol.Ar;
mass.H = mol.H * massMol.H;
mass.N = mol.N * massMol.N;

massToT = mass.HE + mass.O + mass.N2 + mass.O2 + mass.Ar + mass.H + mass.N;

massFrac.HE = mass.HE / massToT;
massFrac.O = mass.O / massToT;
massFrac.N2 = mass.N2 / massToT;
massFrac.O2 = mass.O2 / massToT;
massFrac.Ar = mass.Ar / massToT;
massFrac.H = mass.H / massToT;
massFrac.N = mass.N / massToT;


gamma.HE = 5/3;
gamma.O  = 5/3;
gamma.N2 = 1.4;
gamma.O2 = 1.4;
gamma.Ar = 5/3;
gamma.H  = 5/3;
gamma.N  = 5/3;

R.HE     = 8314/massMol.HE;
R.O      = 8314/massMol.O;
R.N2     = 8314/massMol.N2;
R.O2     = 8314/massMol.O2;
R.Ar     = 8314/massMol.Ar;
R.H      = 8314/massMol.H;
R.N      = 8314/massMol.N;

cv.HE    = R.HE/(gamma.HE-1);
cv.O     = R.O /(gamma.O- 1);
cv.N2    = R.N2/(gamma.N2-1);
cv.O2    = R.O2/(gamma.O2-1);
cv.Ar    = R.Ar/(gamma.Ar-1);
cv.H     = R.H /(gamma.H- 1);
cv.N     = R.N /(gamma.N- 1);

cp.HE    = R.HE/(gamma.HE-1) + R.HE;
cp.O     = R.O /(gamma.O- 1) + R.O;
cp.N2    = R.N2/(gamma.N2-1) + R.N2;
cp.O2    = R.O2/(gamma.O2-1) + R.O2;
cp.Ar    = R.Ar/(gamma.Ar-1) + R.Ar;
cp.H     = R.H /(gamma.H- 1) + R.H;
cp.N     = R.N /(gamma.N- 1) + R.N;


cpMix      = massFrac.HE * cp.HE + massFrac.O * cp.O + massFrac.N2 * cp.N2 + massFrac.O2 * cp.O2 + massFrac.Ar * cp.Ar + massFrac.H * cp.H + massFrac.N * cp.N;
cvMix      = massFrac.HE * cv.HE + massFrac.O * cv.O + massFrac.N2 * cv.N2 + massFrac.O2 * cv.O2 + massFrac.Ar * cv.Ar + massFrac.H * cv.H + massFrac.N * cv.N;
gammaMix   = cpMix/cvMix;
massMolMix = massFrac.HE * massMol.HE + massFrac.O * massMol.O + massFrac.N2 * massMol.N2 + massFrac.O2 * massMol.O2 + massFrac.Ar * massMol.Ar + massFrac.H * massMol.H + massFrac.N * massMol.N;
RMix       = 8314/massMolMix;

soundspeed= sqrt(Temp*RMix*gammaMix);
end