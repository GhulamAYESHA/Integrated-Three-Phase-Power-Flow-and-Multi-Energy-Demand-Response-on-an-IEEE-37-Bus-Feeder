function [from, to, Z, loads, N] = read_ieee_excel_3ph(xlsx)
% === Reads IEEE 37-bus Excel (numeric-only, no headers) ===
% Robust against merged or non-English headers

fprintf('=== Reading IEEE Excel Data (numeric mode) ===\n');

% --- Step 1: Read numeric content only ---
[num, txt, raw] = xlsread(xlsx, 1); % Sheet 1 by default
fprintf('Read %d rows and %d columns (numeric only)\n\n', size(num,1), size(num,2));

% --- Step 2: Define column indices ---
% Adjust these if your Excel layout changes
iFROM = 2; 
iTO   = 3;

% Resistance ()
iRaa = 4; iRbb = 5; iRcc = 6;
iRab = 7; iRac = 8; iRbc = 9;

% Reactance ()
iXaa = 10; iXbb = 11; iXcc = 12;
iXab = 13; iXac = 14; iXbc = 15;

% Loads (kW, kvar)
iPa = 16; iPb = 17; iPc = 18;
iQa = 19; iQb = 20; iQc = 21;

% --- Step 3: Extract data safely ---
nCols = size(num, 2);

% Guard for short columns
col = @(i) (i <= nCols) * num(:,i) + (i > nCols) * 0;

from = col(iFROM);
to   = col(iTO);

Z.Raa = col(iRaa); Z.Rbb = col(iRbb); Z.Rcc = col(iRcc);
Z.Rab = col(iRab); Z.Rac = col(iRac); Z.Rbc = col(iRbc);

Z.Xaa = col(iXaa); Z.Xbb = col(iXbb); Z.Xcc = col(iXcc);
Z.Xab = col(iXab); Z.Xac = col(iXac); Z.Xbc = col(iXbc);

P0a = col(iPa); P0b = col(iPb); P0c = col(iPc);
Q0a = col(iQa); Q0b = col(iQb); Q0c = col(iQc);

% --- Step 4: Load packaging ---
loads.PQperPhase_kWkvar = {P0a, P0b, P0c, Q0a, Q0b, Q0c};

% --- Step 5: Compute N ---
N = max(max(from), max(to));

% --- Step 6: Summary printout ---
fprintf('=== Data Summary ===\n');
fprintf('Buses: %d\n', N);
fprintf('Branches: %d\n\n', length(from));

fprintf('Branch | From | To   | Raa      | Xaa      | Pa      | Qa\n');
fprintf('-------|------|------|----------|----------|---------|--------\n');
for i = 1:min(3, length(from))
    fprintf('  %2d   |  %2d  |  %2d  | %.4f | %.4f | %6.1f | %6.1f\n', ...
        i, from(i), to(i), Z.Raa(i), Z.Xaa(i), P0a(i), Q0a(i));
end

fprintf('\nTotal loads:\n');
fprintf('  Phase A: P = %.2f kW, Q = %.2f kvar\n', sum(P0a), sum(Q0a));
fprintf('  Phase B: P = %.2f kW, Q = %.2f kvar\n', sum(P0b), sum(Q0b));
fprintf('  Phase C: P = %.2f kW, Q = %.2f kvar\n', sum(P0c), sum(Q0c));
fprintf('  Total:   P = %.2f kW, Q = %.2f kvar\n\n', ...
    sum(P0a+P0b+P0c), sum(Q0a+Q0b+Q0c));

end