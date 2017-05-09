function num = precalc_states(n, m)
global states_id
k = 1;
for i = 0:(2 ^ n - 1)
    if ((sum(dec2bin(i) == '1') == 2) && (sum(dec2bin(bitand(i, bitshift(i, 1))) == '1') == 1))
        states_id.adjacent(i + 1) = 1;
    else
        states_id.adjacent(i + 1) = 0;
    end
    if (sum(dec2bin(i) == '1') == m)
        states_id.xtoid(i + 1) = k;
        states_id.idtox(k) = i;
        k = k + 1;
    else
        states_id.xtoid(i + 1) = 0;
    end
end
num = k - 1;
end