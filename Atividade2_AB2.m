%% Initialization
clear ; close all; clc
syms T K z;
syms Fz(x);

% F(z) = z^3 + a2*z^2 + a1*z + a0 = 0
a0 = 1.395 * 10^(-4) * K * T^3 + 16.74*T - 111.6*T^2 - 1;
a1 = 3 - 33.48*T + 1.395 * 10^(-4) * K * T^3;
a2 = 111.6 * T^2 + 16.74 * T - 3;
a3 = 1;

Fz_coeff = [a3 a2 a1 a0];
F(z) = a3*z^3 + a2*z^2 + a1*z + a0;

conds(1) = F(1) > 1;
conds(2) = F(-1)*(-1)^(length(Fz_coeff)-1) > 1;

for i = 1:length(Fz_coeff)
    jTable(1,i) = Fz_coeff(length(Fz_coeff)+1-i);
end

size = 0;
for i = 2:2*length(Fz_coeff)-5
    sizeTable = length(Fz_coeff)-size;
    if mod(i,2) == 1
        for j = 1:sizeTable
            jTable(i,j) = jTable(i-2,1)*jTable(i-1,sizeTable-j+2) - jTable(i-2,sizeTable+1)*jTable(i-1,j);
        end
    else
        for j = 1:sizeTable
            jTable(i,j) = jTable(i-1,sizeTable+1-j);
        end
        size = size + 1;
    end
end

conds(3) = abs(jTable(1,1)) < abs(jTable(1,length(Fz_coeff)));
idx = 4;
shift = 1;
for i = 3:2:length(jTable)
    conds(idx) = abs(jTable(i,1)) > abs(jTable(i,length(jTable)-shift));
    shift = shift + 1;
    idx = idx + 1;
end

conds(idx+1) = T > 0;
conds(idx+2) = K > 0;
disp(vpa(jTable,3))

S = solve( conds , [K T], 'IgnoreAnalyticConstraints', true, 'Real', true, 'ReturnConditions', true);
