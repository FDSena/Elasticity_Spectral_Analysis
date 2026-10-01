function verify_matrix_3D(A, pen)

    fprintf('\n=== Vérification de la Matrice 3D ===\n');

    sym_err = norm(A - A', 'fro') / norm(A, 'fro');
    fprintf('Erreur de symétrie relative : %.2e\n', sym_err);

    diag_A = full(diag(A));
    fixed_indices = find(abs(diag_A - pen) < pen*1e-3);
    fprintf('Nombre de DDL fixés détectés : %d\n', length(fixed_indices));

    free_indices = setdiff(1:size(A,1), fixed_indices);
    A_free = A(free_indices, free_indices);

    opts.tol = 1e-8;
    opts.maxit = 1000;
    opts.isreal = true;

    min_eig = eigs(A_free, 1, 'smallestreal', opts);

    fprintf('Plus petite valeur propre (partie libre) : %.4e\n', min_eig);

    if min_eig > 0
        fprintf('OK : La matrice est définie positive.\n');
    else
        fprintf('ATTENTION : La matrice semble indéfinie.\n');
    end
end