function visualize_modes_3D(eigenvectors, eigenvalues, nx, ny, nz, L, H, Width, rho)

    n_modes = size(eigenvectors, 2);

    idx = @(i,j,k,c) 3 * ((k-1)*nx*ny + (j-1)*nx + (i-1)) + c;

    x = linspace(0, L, nx);
    y = linspace(0, H, ny);
    z = linspace(0, Width, nz);

    [X, Y, Z] = meshgrid(x, y, z);

    for m = 1:n_modes

        U = zeros(size(X));
        V = zeros(size(X));
        Wdisp = zeros(size(X));

        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    U(j,i,k)     = eigenvectors(idx(i,j,k,1), m);
                    V(j,i,k)     = eigenvectors(idx(i,j,k,2), m);
                    Wdisp(j,i,k) = eigenvectors(idx(i,j,k,3), m);
                end
            end
        end

        amp = max(sqrt(U.^2 + V.^2 + Wdisp.^2), [], 'all');
        alpha = 0.15 / amp;

        Xdef = X + alpha * U;
        Ydef = Y + alpha * V;
        Zdef = Z + alpha * Wdisp;

        freq = sqrt(abs(eigenvalues(m)) / rho) / (2*pi);

        figure('Position', [100 100 900 600]);
        hold on;

        % Arêtes selon X
        for k = 1:nz
            for j = 1:ny
                plot3(Xdef(j,:,k), Zdef(j,:,k), Ydef(j,:,k), 'b', 'LineWidth', 1.2);
            end
        end

        % Arêtes selon Z
        for i = 1:nx
            for j = 1:ny
                plot3(squeeze(Xdef(j,i,:)), squeeze(Zdef(j,i,:)), squeeze(Ydef(j,i,:)), ...
                      'r', 'LineWidth', 1.2);
            end
        end

        % Arêtes selon Y
        for k = 1:nz
            for i = 1:nx
                plot3(squeeze(Xdef(:,i,k)), squeeze(Zdef(:,i,k)), squeeze(Ydef(:,i,k)), ...
                      'k', 'LineWidth', 1.0);
            end
        end

        title(sprintf('Mode 3D %d - f = %.2f Hz', m, freq));

        xlabel('x');
        ylabel('z');
        zlabel('y');

        axis equal;
        grid on;
        view(3);
        rotate3d on;

    end
end