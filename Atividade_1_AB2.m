%% Initialization
clear ; close all; clc
syms T K;

a0 = 1.395 * 10^(-4) * K * T^3 + 16.74*T - 111.6*T^2 - 1;
a1 = 3 - 33.48*T + 1.395 * 10^(-4) * K * T^3;
a2 = 111.6 * T^2 + 16.74 * T - 3;

Fz_roots = roots([1 a2 a1 a0]);

%% Transform Bilinear
for i = 1:length(Fz_roots)
    r = Fz_roots(i);
    coeffVector(i) = (r+1)/(r-1);
end

%coeffVector = [(2.79*10^(-4)*K*T^3) (446.4*T^2 - 5.58*10^(-4)*K*T^3) (-446.4*T^2 + 66.96*T + 2.79*10^(-4)*K*T^3) (8 - 66.96*T)];

%% Compute Routh-Hurwitz
%coeffVector = [8.1 0.9 -0.9 -0.1]%input('input vector of your system coefficients: \n i.e. [an an-1 an-2 ... a0] = ');
ceoffLength = length(coeffVector);
rhTableColumn = round(ceoffLength/2);

%  Initialize Routh-Hurwitz table with empty zero array
rhTable = sym(zeros(ceoffLength,rhTableColumn));

%  Compute first row of the table
rhTable(1,:) = coeffVector(1,1:2:ceoffLength);

%  Check if length of coefficients vector is even or odd
if (rem(ceoffLength,2) ~= 0)
    % if odd, second row of table will be
    rhTable(2,1:rhTableColumn - 1) = coeffVector(1,2:2:ceoffLength);
else
    % if even, second row of table will be
    rhTable(2,:) = coeffVector(1,2:2:ceoffLength);
end

%% Calculate Routh-Hurwitz table's rows

%  Set epss as a small value
epss = 0.01;
%syms epss
%  Calculate other elements of the table
for i = 3:ceoffLength
   
    %  special case: row of all zeros
    if rhTable(i-1,:) == 0
        order = (ceoffLength - i);
        cnt1 = 0;
        cnt2 = 1;
        for j = 1:rhTableColumn - 1
            rhTable(i-1,j) = (order - cnt1) * rhTable(i-2,cnt2);
            cnt2 = cnt2 + 1;
            cnt1 = cnt1 + 2;
        end
    end
    
    for j = 1:rhTableColumn - 1
        %  first element of upper row
        firstElemUpperRow = rhTable(i-1,1);
        
        %  compute each element of the table
        rhTable(i,j) = ((rhTable(i-1,1) * rhTable(i-2,j+1)) - ....
            (rhTable(i-2,1) * rhTable(i-1,j+1))) / firstElemUpperRow;
    end
    
    
    %  special case: zero in the first column
    if rhTable(i,1) == 0
        rhTable(i,1) = epss;
    end
end

%%  Compute number of right hand side poles(unstable poles)
%   Initialize unstable poles with zero
unstablePoles = 0;

%   Check change in signs
for i = 1:ceoffLength - 1
    if sign(rhTable(i,1)) * sign(rhTable(i+1,1)) == -1
        unstablePoles = unstablePoles + 1;
    end
end

%   Print calculated data on screen
fprintf('\n Routh-Hurwitz Table:\n')

disp(vpa(rhTable,3))

for i = 1:length(rhTable)
    conds(i) = rhTable(i,1) > 0;
end
conds(length(rhTable)+1) = T > 0;
conds(length(rhTable)+2) = K > 0;

S = solve(conds, [K T],'IgnoreAnalyticConstraints', true, 'Real', true, 'ReturnConditions', true);
disp(vpa(S.conditions,3))
