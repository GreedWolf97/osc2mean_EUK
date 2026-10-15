function [Ftable, dFtable] = KaulaFTable_EU(Lmax, inc_deg)
% Kaula inclination functions table

inc = deg2rad(inc_deg);
Ftable  = zeros(Lmax + 1, Lmax + 1, Lmax + 1);
dFtable = zeros(Lmax + 1, Lmax + 1, Lmax + 1);

for l = 2:Lmax
    for m = 0:l
        for p = 0:l
            [F, dF] = KaulaF_EU(l, m, p, inc);
            Ftable(l + 1, m + 1, p + 1)  = F;
            dFtable(l + 1, m + 1, p + 1) = dF;
        end
    end
end

end