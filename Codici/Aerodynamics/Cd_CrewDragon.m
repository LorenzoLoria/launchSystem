function Cd = Cd_CrewDragon(M)
%
% INPUT
% - M     : mach (scalare)
%
% OUTPUT
% - Cd    : drag coefficient (scalare)


if M <= 0.722
    Cd = 0.68 * M^2 + 0.424;
else
    Cd = 1.73 * exp(-1.02 * (M + 0.25)) - 8.9 * exp(-1.95 * (M + 0.35)) + 1.23;
end



end