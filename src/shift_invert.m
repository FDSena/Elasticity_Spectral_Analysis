function [eigenvalues, eigenvectors] = shift_invert(A, rho, n_modes, sigma, tol, maxit)
% SHIFT_INVERT_2D
% Calcule les premiers modes propres par itération inverse avec décalage.
% Résolution linéaire par GMRES préconditionné ILU.

    N = size(A, 1);
    As = A - sigma * speye(N);

    eigenvalues  = zeros(n_modes, 1);
    eigenvectors = zeros(N, n_modes);

    % Couper les warnings trop nombreux
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:illConditionedMatrix');

    % Préconditionneur ILU
    opts_ilu.type = 'ilutp';
    opts_ilu.droptol = 1e-3;

    [L_ilu, U_ilu] = ilu(As, opts_ilu);

    fprintf('Préconditionneur ILU calculé.\n');

    restart = min(100, N);
    max_outer = 500;

    for k = 1:n_modes

        v = randn(N, 1);

        for m = 1:k-1
            v = v - (eigenvectors(:,m)' * v) * eigenvectors(:,m);
        end

        v = v / norm(v);

        converged = false;
        theta_old = NaN;
        
        for it = 1:max_outer
        
            [w, flag] = gmres(As, v, restart, 1e-12, maxit, L_ilu, U_ilu);
        
            if flag ~= 0
                [w, flag] = gmres(As, v, restart, 1e-12, maxit*2, L_ilu, U_ilu);
            end
        
            for m = 1:k-1
                w = w - (eigenvectors(:,m)' * w) * eigenvectors(:,m);
            end
        
            if norm(w) < 1e-14
                w = randn(N,1);
            
                for m = 1:k-1
                    w = w - (eigenvectors(:,m)' * w) * eigenvectors(:,m);
                end
            end
            
            v = w / norm(w);        
            theta = v' * A * v;
            res = norm(A*v - theta*v) / norm(A*v);
        
            if ~isnan(theta_old)
                crit_theta = abs(theta - theta_old) / max(1, abs(theta));
        
                if crit_theta < tol && res < tol
                    converged = true;
                    break;
                end
            end
        
            theta_old = theta;
        end
        eigenvalues(k) = theta;
        eigenvectors(:,k) = v;

        omega = sqrt(abs(theta) / rho);
        freq = omega / (2*pi);

        if converged
            fprintf('Mode %d : omega^2 = %.4e  f = %.4f Hz  [convergé en %d iter, résidu = %.2e]\n', ...
                    k, theta, freq, it, res);
        else
            fprintf('Mode %d : omega^2 = %.4e  f = %.4f Hz  [NON convergé, résidu = %.2e]\n', ...
                    k, theta, freq, res);
        end
    end

    warning('on', 'MATLAB:nearlySingularMatrix');
    warning('on', 'MATLAB:singularMatrix');
    warning('on', 'MATLAB:illConditionedMatrix');

end