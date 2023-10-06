for n = 1:100

    for kz = 2:3

        for niter = 1:4

            for k = kz:2:Kt - 1
                Jac = zeros(size(X,1));
                G = zeros([size(X, 1), 1]);
                % Compute residual and fill the tridiagonal Jacobian
                % M (even points)
                for j = 2:2:Jt - 1 % M
                    dsdY = X(j + 1, k) - X(j - 1, k);
                    dMdZZ = X(j, k - 1) - 2 * X(j, k) + X(j, k + 1);
                    Jac(j, j) = -2 * dsdY; % D
                    Jac(j, j - 1) = -dMdZZ; % L
                    Jac(j, j + 1) = dMdZZ; % U
                    G(j) = dsdY * dMdZZ ...
                        + Gamma(j, k) * Sigma(j, k) * fac2; 
                end

                % s (odd points)
                for j = 3:2:Jt - 2 % s
                    dMdY = X(j + 1, k) - X(j - 1, k);
                    Jac(j, j) = fac1(j) * (2 * X(j, k)) * (1 - Y(j) ^ 2) / (1 - X(j, k) ^ 2) ^ 2; % D
                    Jac(j, j - 1) = -1; % L
                    Jac(j, j + 1) = 1; % U
                    G(j) = dMdY - fac1(j) * (Y(j) ^ 2 - X(j, k) ^ 2) / (1 - X(j, k) ^ 2);
                end

                X(2:end - 1, k) = X(2:end - 1, k) - (Jac(2:end - 1, 2:end - 1) \ G(2:end - 1));
            end

        end

    end

    % updateBoundary

end
