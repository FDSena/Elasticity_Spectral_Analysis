function A = assemble_matrix_2D(nx, ny, h, lambda, mu, pen)
% ASSEMBLE_MATRIX_2D
% Assemble la matrice d'élasticité 2D pour les équations de Navier.
%
% Entrées :
%   nx, ny : nombre de points de grille
%   h      : pas de discrétisation
%   lambda : premier coefficient de Lamé
%   mu     : second coefficient de Lamé
%   pen    : coefficient de pénalisation
%
% Sortie :
%   A : matrice sparse d'élasticité

    if nargin < 6
        pen = 1e6;
    end

    N = 2 * nx * ny;

    idx_u = @(i,j) 2*((j-1)*nx + (i-1)) + 1;
    idx_v = @(i,j) 2*((j-1)*nx + (i-1)) + 2;

    c1 = lambda + 2*mu;
    c2 = mu;
    c3 = (lambda + mu) / (4*h^2);

    nnz_est = 30 * N;

    rows = zeros(nnz_est,1);
    cols = zeros(nnz_est,1);
    vals = zeros(nnz_est,1);

    ptr = 1;

    for j = 1:ny
        for i = 1:nx

            iu = idx_u(i,j);
            iv = idx_v(i,j);

            cross_ok = (i > 1) && (i < nx) && ...
                       (j > 1) && (j < ny);

            % ====================================
            % Equation sur u
            % ====================================

            diag_u = 0;

            % d²u/dx²
            if i+1 <= nx
                rows(ptr)=iu;
                cols(ptr)=idx_u(i+1,j);
                vals(ptr)=-c1/h^2;
                ptr=ptr+1;

                diag_u = diag_u + c1/h^2;
            end

            if i-1 >= 1
                rows(ptr)=iu;
                cols(ptr)=idx_u(i-1,j);
                vals(ptr)=-c1/h^2;
                ptr=ptr+1;

                diag_u = diag_u + c1/h^2;
            end

            % d²u/dy²
            if j+1 <= ny
                rows(ptr)=iu;
                cols(ptr)=idx_u(i,j+1);
                vals(ptr)=-c2/h^2;
                ptr=ptr+1;

                diag_u = diag_u + c2/h^2;
            end

            if j-1 >= 1
                rows(ptr)=iu;
                cols(ptr)=idx_u(i,j-1);
                vals(ptr)=-c2/h^2;
                ptr=ptr+1;

                diag_u = diag_u + c2/h^2;
            end

            % dérivée croisée
            if cross_ok

                rows(ptr)=iu;
                cols(ptr)=idx_v(i+1,j+1);
                vals(ptr)=-c3;
                ptr=ptr+1;

                rows(ptr)=iu;
                cols(ptr)=idx_v(i-1,j+1);
                vals(ptr)=c3;
                ptr=ptr+1;

                rows(ptr)=iu;
                cols(ptr)=idx_v(i+1,j-1);
                vals(ptr)=c3;
                ptr=ptr+1;

                rows(ptr)=iu;
                cols(ptr)=idx_v(i-1,j-1);
                vals(ptr)=-c3;
                ptr=ptr+1;
            end

            rows(ptr)=iu;
            cols(ptr)=iu;
            vals(ptr)=diag_u;
            ptr=ptr+1;

            % ====================================
            % Equation sur v
            % ====================================

            diag_v = 0;

            % d²v/dx²
            if i+1 <= nx
                rows(ptr)=iv;
                cols(ptr)=idx_v(i+1,j);
                vals(ptr)=-c2/h^2;
                ptr=ptr+1;

                diag_v = diag_v + c2/h^2;
            end

            if i-1 >= 1
                rows(ptr)=iv;
                cols(ptr)=idx_v(i-1,j);
                vals(ptr)=-c2/h^2;
                ptr=ptr+1;

                diag_v = diag_v + c2/h^2;
            end

            % d²v/dy²
            if j+1 <= ny
                rows(ptr)=iv;
                cols(ptr)=idx_v(i,j+1);
                vals(ptr)=-c1/h^2;
                ptr=ptr+1;

                diag_v = diag_v + c1/h^2;
            end

            if j-1 >= 1
                rows(ptr)=iv;
                cols(ptr)=idx_v(i,j-1);
                vals(ptr)=-c1/h^2;
                ptr=ptr+1;

                diag_v = diag_v + c1/h^2;
            end

            % dérivée croisée
            if cross_ok

                rows(ptr)=iv;
                cols(ptr)=idx_u(i+1,j+1);
                vals(ptr)=-c3;
                ptr=ptr+1;

                rows(ptr)=iv;
                cols(ptr)=idx_u(i-1,j+1);
                vals(ptr)=c3;
                ptr=ptr+1;

                rows(ptr)=iv;
                cols(ptr)=idx_u(i+1,j-1);
                vals(ptr)=c3;
                ptr=ptr+1;

                rows(ptr)=iv;
                cols(ptr)=idx_u(i-1,j-1);
                vals(ptr)=-c3;
                ptr=ptr+1;
            end

            rows(ptr)=iv;
            cols(ptr)=iv;
            vals(ptr)=diag_v;
            ptr=ptr+1;

        end
    end

    rows = rows(1:ptr-1);
    cols = cols(1:ptr-1);
    vals = vals(1:ptr-1);

    A = sparse(rows, cols, vals, N, N);

    % ====================================
    % Conditions de Dirichlet homogènes
    % sur les deux extremités du bord inférieur
    % ====================================

    fixed_dof = [idx_u(1,1), idx_v(1,1), idx_u(nx,1), idx_v(nx,1)];

    % pénalisation symétrique
    A(fixed_dof,:) = 0;
    A(:,fixed_dof) = 0;

    for k = 1:length(fixed_dof)
        A(fixed_dof(k), fixed_dof(k)) = pen;
    end

    % symétrisation finale
    A = 0.5 * (A + A');

end