function visualize_modes(eigenvectors, eigenvalues, nx, ny, L, H, rho)

    n_modes = size(eigenvectors, 2);

    idx_u = @(i,j) 2*((j-1)*nx + (i-1)) + 1;
    idx_v = @(i,j) 2*((j-1)*nx + (i-1)) + 2;

    x = linspace(0, L, nx);
    y = linspace(0, H, ny);
    [X, Y] = meshgrid(x, y);

    figure('Units', 'normalized', 'OuterPosition', [0 0 0.7 1]);

    for k = 1:n_modes

        U_mode = zeros(ny, nx);
        V_mode = zeros(ny, nx);

        for j = 1:ny
            for i = 1:nx
                U_mode(j,i) = eigenvectors(idx_u(i,j), k);
                V_mode(j,i) = eigenvectors(idx_v(i,j), k);
            end
        end

        amp = max(sqrt(U_mode(:).^2 + V_mode(:).^2));

        if amp == 0
            alpha = 1;
        else
            alpha = 0.3 * H / amp;
        end

        Xdef = X + alpha * U_mode;
        Ydef = Y + alpha * V_mode;

        freq = sqrt(abs(eigenvalues(k)) / rho) / (2*pi);

        subplot(n_modes, 1, k);
        hold on;

        for j = 1:ny
            plot(X(j,:), Y(j,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
            plot(Xdef(j,:), Ydef(j,:), 'b-', 'LineWidth', 1.5);
        end

        for i = 1:nx
            plot(X(:,i), Y(:,i), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
            plot(Xdef(:,i), Ydef(:,i), 'b-', 'LineWidth', 1.5);
        end

        amp_map = sqrt(U_mode.^2 + V_mode.^2);

        scatter(Xdef(:), Ydef(:), 25, amp_map(:), 'filled');

        quiver(X, Y, U_mode, V_mode, 0.8, 'k');

        colormap(jet);
        colorbar;

        title(sprintf('Mode %d — f = %.1f Hz', k, freq), 'FontSize', 10);

        axis equal;
        axis tight;
        grid on;

        xlabel('x (m)');
        ylabel('y (m)');

        hold off;
    end

    sgtitle('Modes propres du pont : maillage déformé et champ de déplacement', ...
            'FontSize', 12);

end