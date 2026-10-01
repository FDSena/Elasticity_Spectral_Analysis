function animate_modes_3D(eigvecs, eigvals, nx, ny, nz, L, H, Width, rho, modes, alpha_max)

    if nargin < 10
        modes = 1:size(eigvecs,2);
    end

    if nargin < 11
        alpha_max = 0.4;
    end

    idx = @(i,j,k,c) 3*((k-1)*nx*ny + (j-1)*nx + (i-1)) + c;

    x = linspace(0,L,nx);
    y = linspace(0,H,ny);
    z = linspace(0,Width,nz);

    [X,Y,Z] = meshgrid(x,y,z);

    nframes = 80;

    for mode = modes

        U  = zeros(ny,nx,nz);
        V  = zeros(ny,nx,nz);
        Wd = zeros(ny,nx,nz);

        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    U(j,i,k)  = eigvecs(idx(i,j,k,1),mode);
                    V(j,i,k)  = eigvecs(idx(i,j,k,2),mode);
                    Wd(j,i,k) = eigvecs(idx(i,j,k,3),mode);
                end
            end
        end

        amp = max(sqrt(U.^2 + V.^2 + Wd.^2), [], 'all');
        alpha = alpha_max / amp;

        freq = sqrt(abs(eigvals(mode))/rho)/(2*pi);

        filename = sprintf('mode_%d_3D_animation.gif', mode);

        figure('Visible','off', ...
               'Color','w', ...
               'Position',[100 100 900 650]);
                for frame_id = 1:nframes

            t = 2*pi*frame_id/nframes;
            s = sin(t);

            Xd = X + alpha*s*U;
            Yd = Y + alpha*s*V;
            Zd = Z + alpha*s*Wd;

            clf
            hold on

            % Lignes selon x
            for k = 1:nz
                for j = 1:ny
                    plot3(squeeze(Xd(j,:,k)), ...
                          squeeze(Zd(j,:,k)), ...
                          squeeze(Yd(j,:,k)), ...
                          'b','LineWidth',1.4);
                end
            end

            % Lignes selon z
            for i = 1:nx
                for j = 1:ny
                    plot3(squeeze(Xd(j,i,:)), ...
                          squeeze(Zd(j,i,:)), ...
                          squeeze(Yd(j,i,:)), ...
                          'r','LineWidth',1.4);
                end
            end

            % Lignes selon y
            for i = 1:nx
                for k = 1:nz
                    plot3(squeeze(Xd(:,i,k)), ...
                          squeeze(Zd(:,i,k)), ...
                          squeeze(Yd(:,i,k)), ...
                          'k','LineWidth',1.0);
                end
            end

            axis equal
            grid on
            xlabel('x')
            ylabel('z')
            zlabel('y')

            title(sprintf('Mode 3D %d - f = %.2f Hz', mode, freq));

            view(35,25)

            xlim([-1 5])
            ylim([-0.5 Width+0.5])
            zlim([-0.8 H+0.8])

            drawnow

            img = frame2im(getframe(gcf));
            [A,map] = rgb2ind(img,256);

            if frame_id == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
            end
        end

        fprintf('Animation créée : %s\n', filename);

    end
end