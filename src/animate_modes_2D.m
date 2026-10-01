function animate_modes_2D(eigvecs, eigvals, nx, ny, L, H, rho, modes, alpha_max)

    if nargin < 8
        modes = 1:size(eigvecs,2);
    end

    if nargin < 9
        alpha_max = 0.5;
    end

    x = linspace(0, L, nx);
    y = linspace(0, H, ny);
    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));

    for mode = modes

        v = eigvecs(:,mode);

        Ux = zeros(ny,nx);
        Uy = zeros(ny,nx);

        for j = 1:ny
            for i = 1:nx
                idu = 2*((j-1)*nx + (i-1)) + 1;
                idv = idu + 1;

                Ux(j,i) = v(idu);
                Uy(j,i) = v(idv);
            end
        end

        amp = max(sqrt(Ux.^2 + Uy.^2), [], 'all');
        alpha = alpha_max / amp;

        omega = sqrt(abs(eigvals(mode)) / rho);
        freq = omega / (2*pi);

        filename = sprintf('mode_%d_animation.gif', mode);

        figure('Position',[100 100 1000 400])

        for n = 1:120

            t = 2*pi*n/120;
            scale = alpha * sin(t);

            Xdef = X + scale*Ux;
            Ydef = Y + scale*Uy;

            clf

            mesh(X,Y,Z,'EdgeColor',[0.8 0.8 0.8])
            hold on
            mesh(Xdef,Ydef,Z,'EdgeColor','b','LineWidth',1.5)

            axis equal
            axis([-0.5 L+0.5 -0.8 H+0.8])
            view(2)

            title(sprintf('Mode %d - f = %.2f Hz', mode, freq))

            drawnow

            frame = getframe(gcf);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);

            if n == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.04);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.04);
            end
        end

        fprintf('Animation créée : %s\n', filename);

    end
end