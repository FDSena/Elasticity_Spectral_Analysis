function verify_matrix(A, pen)

    diag_A   = full(diag(A));
    free_dof = find(diag_A < pen * 0.1);

    fprintf('DDL libres : %d / %d\n', length(free_dof), size(A,1));

    A_free = A(free_dof, free_dof);

    % Symétrie sur les DDL libres
    sym_err = norm(A_free - A_free', 'fro') / norm(A_free, 'fro');
    fprintf('Erreur de symétrie relative sur A_free : %.2e', sym_err);

    if sym_err < 1e-10
        fprintf('  -> OK\n');
    else
        fprintf('  -> ATTENTION : matrice libre non symétrique\n');
    end

    % Définie positive sur les DDL libres
    opts.tol   = 1e-6;
    opts.maxit = 1000;
    opts.isreal = true;
    opts.issym  = true;

    lambda_min = eigs(A_free, 1, 'smallestreal', opts);

    fprintf('Plus petite valeur propre (DDL libres) : %.4e', lambda_min);

    if lambda_min > 0
        fprintf('  -> définie positive OK\n');
    else
        fprintf('  -> ATTENTION : valeur propre <= 0\n');
    end

end