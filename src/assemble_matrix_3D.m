function A = assemble_matrix_3D(nx, ny, nz, hx, hy, hz, lambda, mu, pen)

    if nargin < 9
        pen = 1e12;
    end

    N = 3 * nx * ny * nz;

    idx = @(i,j,k,c) 3*((k-1)*nx*ny + (j-1)*nx + (i-1)) + c;

    c1 = lambda + 2*mu;
    c2 = mu;

    cxy = (lambda + mu) / (4*hx*hy);
    cxz = (lambda + mu) / (4*hx*hz);
    cyz = (lambda + mu) / (4*hy*hz);

    nnz_est = 50 * N;

    rows = zeros(nnz_est,1);
    cols = zeros(nnz_est,1);
    vals = zeros(nnz_est,1);

    ptr = 1;

    for k = 1:nz
        for j = 1:ny
            for i = 1:nx

                iu = idx(i,j,k,1);
                iv = idx(i,j,k,2);
                iw = idx(i,j,k,3);

                % =========================
                % Equation sur u
                % =========================
                diag_u = 0;

                if i+1 <= nx
                    rows(ptr)=iu; cols(ptr)=idx(i+1,j,k,1); vals(ptr)=-c1/hx^2; ptr=ptr+1;
                    diag_u = diag_u + c1/hx^2;
                end
                if i-1 >= 1
                    rows(ptr)=iu; cols(ptr)=idx(i-1,j,k,1); vals(ptr)=-c1/hx^2; ptr=ptr+1;
                    diag_u = diag_u + c1/hx^2;
                end
                if j+1 <= ny
                    rows(ptr)=iu; cols(ptr)=idx(i,j+1,k,1); vals(ptr)=-c2/hy^2; ptr=ptr+1;
                    diag_u = diag_u + c2/hy^2;
                end
                if j-1 >= 1
                    rows(ptr)=iu; cols(ptr)=idx(i,j-1,k,1); vals(ptr)=-c2/hy^2; ptr=ptr+1;
                    diag_u = diag_u + c2/hy^2;
                end
                if k+1 <= nz
                    rows(ptr)=iu; cols(ptr)=idx(i,j,k+1,1); vals(ptr)=-c2/hz^2; ptr=ptr+1;
                    diag_u = diag_u + c2/hz^2;
                end
                if k-1 >= 1
                    rows(ptr)=iu; cols(ptr)=idx(i,j,k-1,1); vals(ptr)=-c2/hz^2; ptr=ptr+1;
                    diag_u = diag_u + c2/hz^2;
                end

                % - (lambda+mu) d²v/dxdy
                if i>1 && i<nx && j>1 && j<ny
                    rows(ptr)=iu; cols(ptr)=idx(i+1,j+1,k,2); vals(ptr)=-cxy; ptr=ptr+1;
                    rows(ptr)=iu; cols(ptr)=idx(i-1,j+1,k,2); vals(ptr)= cxy; ptr=ptr+1;
                    rows(ptr)=iu; cols(ptr)=idx(i+1,j-1,k,2); vals(ptr)= cxy; ptr=ptr+1;
                    rows(ptr)=iu; cols(ptr)=idx(i-1,j-1,k,2); vals(ptr)=-cxy; ptr=ptr+1;
                end

                % - (lambda+mu) d²w/dxdz
                if i>1 && i<nx && k>1 && k<nz
                    rows(ptr)=iu; cols(ptr)=idx(i+1,j,k+1,3); vals(ptr)=-cxz; ptr=ptr+1;
                    rows(ptr)=iu; cols(ptr)=idx(i-1,j,k+1,3); vals(ptr)= cxz; ptr=ptr+1;
                    rows(ptr)=iu; cols(ptr)=idx(i+1,j,k-1,3); vals(ptr)= cxz; ptr=ptr+1;
                    rows(ptr)=iu; cols(ptr)=idx(i-1,j,k-1,3); vals(ptr)=-cxz; ptr=ptr+1;
                end

                rows(ptr)=iu; cols(ptr)=iu; vals(ptr)=diag_u; ptr=ptr+1;

                % =========================
                % Equation sur v
                % =========================
                diag_v = 0;

                if i+1 <= nx
                    rows(ptr)=iv; cols(ptr)=idx(i+1,j,k,2); vals(ptr)=-c2/hx^2; ptr=ptr+1;
                    diag_v = diag_v + c2/hx^2;
                end
                if i-1 >= 1
                    rows(ptr)=iv; cols(ptr)=idx(i-1,j,k,2); vals(ptr)=-c2/hx^2; ptr=ptr+1;
                    diag_v = diag_v + c2/hx^2;
                end
                if j+1 <= ny
                    rows(ptr)=iv; cols(ptr)=idx(i,j+1,k,2); vals(ptr)=-c1/hy^2; ptr=ptr+1;
                    diag_v = diag_v + c1/hy^2;
                end
                if j-1 >= 1
                    rows(ptr)=iv; cols(ptr)=idx(i,j-1,k,2); vals(ptr)=-c1/hy^2; ptr=ptr+1;
                    diag_v = diag_v + c1/hy^2;
                end
                if k+1 <= nz
                    rows(ptr)=iv; cols(ptr)=idx(i,j,k+1,2); vals(ptr)=-c2/hz^2; ptr=ptr+1;
                    diag_v = diag_v + c2/hz^2;
                end
                if k-1 >= 1
                    rows(ptr)=iv; cols(ptr)=idx(i,j,k-1,2); vals(ptr)=-c2/hz^2; ptr=ptr+1;
                    diag_v = diag_v + c2/hz^2;
                end

                % - (lambda+mu) d²u/dydx
                if j>1 && j<ny && i>1 && i<nx
                    rows(ptr)=iv; cols(ptr)=idx(i+1,j+1,k,1); vals(ptr)=-cxy; ptr=ptr+1;
                    rows(ptr)=iv; cols(ptr)=idx(i-1,j+1,k,1); vals(ptr)= cxy; ptr=ptr+1;
                    rows(ptr)=iv; cols(ptr)=idx(i+1,j-1,k,1); vals(ptr)= cxy; ptr=ptr+1;
                    rows(ptr)=iv; cols(ptr)=idx(i-1,j-1,k,1); vals(ptr)=-cxy; ptr=ptr+1;
                end

                % - (lambda+mu) d²w/dydz
                if j>1 && j<ny && k>1 && k<nz
                    rows(ptr)=iv; cols(ptr)=idx(i,j+1,k+1,3); vals(ptr)=-cyz; ptr=ptr+1;
                    rows(ptr)=iv; cols(ptr)=idx(i,j-1,k+1,3); vals(ptr)= cyz; ptr=ptr+1;
                    rows(ptr)=iv; cols(ptr)=idx(i,j+1,k-1,3); vals(ptr)= cyz; ptr=ptr+1;
                    rows(ptr)=iv; cols(ptr)=idx(i,j-1,k-1,3); vals(ptr)=-cyz; ptr=ptr+1;
                end

                rows(ptr)=iv; cols(ptr)=iv; vals(ptr)=diag_v; ptr=ptr+1;

                % =========================
                % Equation sur w
                % =========================
                diag_w = 0;

                if i+1 <= nx
                    rows(ptr)=iw; cols(ptr)=idx(i+1,j,k,3); vals(ptr)=-c2/hx^2; ptr=ptr+1;
                    diag_w = diag_w + c2/hx^2;
                end
                if i-1 >= 1
                    rows(ptr)=iw; cols(ptr)=idx(i-1,j,k,3); vals(ptr)=-c2/hx^2; ptr=ptr+1;
                    diag_w = diag_w + c2/hx^2;
                end
                if j+1 <= ny
                    rows(ptr)=iw; cols(ptr)=idx(i,j+1,k,3); vals(ptr)=-c2/hy^2; ptr=ptr+1;
                    diag_w = diag_w + c2/hy^2;
                end
                if j-1 >= 1
                    rows(ptr)=iw; cols(ptr)=idx(i,j-1,k,3); vals(ptr)=-c2/hy^2; ptr=ptr+1;
                    diag_w = diag_w + c2/hy^2;
                end
                if k+1 <= nz
                    rows(ptr)=iw; cols(ptr)=idx(i,j,k+1,3); vals(ptr)=-c1/hz^2; ptr=ptr+1;
                    diag_w = diag_w + c1/hz^2;
                end
                if k-1 >= 1
                    rows(ptr)=iw; cols(ptr)=idx(i,j,k-1,3); vals(ptr)=-c1/hz^2; ptr=ptr+1;
                    diag_w = diag_w + c1/hz^2;
                end

                % - (lambda+mu) d²u/dzdx
                if k>1 && k<nz && i>1 && i<nx
                    rows(ptr)=iw; cols(ptr)=idx(i+1,j,k+1,1); vals(ptr)=-cxz; ptr=ptr+1;
                    rows(ptr)=iw; cols(ptr)=idx(i-1,j,k+1,1); vals(ptr)= cxz; ptr=ptr+1;
                    rows(ptr)=iw; cols(ptr)=idx(i+1,j,k-1,1); vals(ptr)= cxz; ptr=ptr+1;
                    rows(ptr)=iw; cols(ptr)=idx(i-1,j,k-1,1); vals(ptr)=-cxz; ptr=ptr+1;
                end

                % - (lambda+mu) d²v/dzdy
                if k>1 && k<nz && j>1 && j<ny
                    rows(ptr)=iw; cols(ptr)=idx(i,j+1,k+1,2); vals(ptr)=-cyz; ptr=ptr+1;
                    rows(ptr)=iw; cols(ptr)=idx(i,j-1,k+1,2); vals(ptr)= cyz; ptr=ptr+1;
                    rows(ptr)=iw; cols(ptr)=idx(i,j+1,k-1,2); vals(ptr)= cyz; ptr=ptr+1;
                    rows(ptr)=iw; cols(ptr)=idx(i,j-1,k-1,2); vals(ptr)=-cyz; ptr=ptr+1;
                end

                rows(ptr)=iw; cols(ptr)=iw; vals(ptr)=diag_w; ptr=ptr+1;

            end
        end
    end

    rows = rows(1:ptr-1);
    cols = cols(1:ptr-1);
    vals = vals(1:ptr-1);

    A = sparse(rows, cols, vals, N, N);

    % Appuis : quatre coins inférieurs du pont
    fixed_dof = [ ...
        idx(1,1,1,1),  idx(1,1,1,2),  idx(1,1,1,3), ...
        idx(nx,1,1,1), idx(nx,1,1,2), idx(nx,1,1,3), ...
        idx(1,1,nz,1), idx(1,1,nz,2), idx(1,1,nz,3), ...
        idx(nx,1,nz,1), idx(nx,1,nz,2), idx(nx,1,nz,3) ...
    ];

    A(fixed_dof,:) = 0;
    A(:,fixed_dof) = 0;

    for d = 1:length(fixed_dof)
        A(fixed_dof(d), fixed_dof(d)) = pen;
    end

    A = 0.5 * (A + A');

end