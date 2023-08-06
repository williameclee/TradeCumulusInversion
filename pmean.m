function pref = pmean(p, x)
    pref = sum((p(2:end) + p(1:end - 1)) .* (x(2:end) - x(1:end - 1)) / 2);
    pref = pref / 2;
end