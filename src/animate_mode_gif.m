function animate_mode_gif(eigenvector, eigenvalue, nx, ny, L, H, rho, filename, k)
    if nargin < 9
        k = 1;
    end
    if nargin < 8
        filename = 'mode_animation.gif';
    end

    idx_u = @(i,j) 2*((j-1)*nx + (i-1)) + 1;
    idx_v = @(i,j) 2*((j-1)*nx + (i-1)) + 2;

    x = linspace(0, L, nx);
    y = linspace(0, H, ny);
    [X, Y] = meshgrid(x, y);

    U = zeros(ny, nx);
    V = zeros(ny, nx);
    for j = 1:ny
        for i = 1:nx
            U(j,i) = eigenvector(idx_u(i,j));
            V(j,i) = eigenvector(idx_v(i,j));
        end
    end

    amp = max(sqrt(U(:).^2 + V(:).^2));
    alpha = 0.7 * H / amp;
    freq = sqrt(abs(eigenvalue) / rho) / (2*pi);

    fig = figure('Visible','off');
    n_frames = 40;

    for frame = 1:n_frames
        clf;
        hold on;
        phase = sin(2*pi*frame/n_frames);
        Xdef = X + alpha * phase * U;
        Ydef = Y + alpha * phase * V;

        for j = 1:ny
            plot(Xdef(j,:), Ydef(j,:), 'b-', 'LineWidth', 1.5);
        end
        for i = 1:nx
            plot(Xdef(:,i), Ydef(:,i), 'b-', 'LineWidth', 1.5);
        end

        title(sprintf('Mode %d — f = %.1f Hz — frame %d', k, freq, frame));
        xlabel('x (m)'); ylabel('y (m)');
        axis equal;
        axis([-0.5, L+0.5, -0.5, H+0.5]);
        grid on;

        drawnow;
        frame_data = getframe(fig);
        im = frame2im(frame_data);
        [imind, cm] = rgb2ind(im, 256);

        if frame == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
        end
    end

    close(fig);
    fprintf('GIF sauvegardé : %s\n', filename);
end