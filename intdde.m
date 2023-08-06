function F = intdde(Fp, F0, Dx)
    F = zeros(size(Fp));
    hby24 = Dx / 24;
    F(1) = F0;
    F(2) = F0 + hby24 * dot([9, 19, -5, 1], Fp(1:4));
    for i = 3:length(Fp) - 1
        F(i) = F(i - 1) + hby24 * dot([-1, 13, 13, -1], Fp(i - 2:i + 1));
    end
    F(end) = F(end - 1) + hby24 * dot([1, -5, 19, 9], Fp(end - 3:end));
end
