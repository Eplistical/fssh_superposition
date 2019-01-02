% exact solution for 1D model

function exact_1d(xI, kxI, sigmax, init_s, A, mass, fID)
    tic;

    c1 = sqrt(1 - init_s);
    c2 = sqrt(init_s);

    N = 2;
    L = 16;
    M = 256;                    

    dt = 0.1;
    Nstep = 3000; %ceil(7.5 / (k(1) / mass) / dt); % 5000;
    tgraph = 100;

    enable_plot = true;
    enable_plot_diab_surf = false;
    enable_plot_surf = false;
    % grids
    x0 = linspace(-L/2, L/2, M)';
    dx = x0(2) - x0(1);
    dkx = 2 * pi / M / dx;
    kx0 = (-M/2:M/2-1)' * dkx;
    % construct TU on k grid
    T = kx0.^2 / 2 / mass;
    TU = exp(-1i * dt * T);
    TU = fftshift(TU);
    % construct VU
    VU = zeros(M,N,N);
    Hs = zeros(M,N,N);
    for j=1:M
        B = 3.0;
        theta = pi / 2 * (erf(B * x0(j)) + 1);

        H = A * [-cos(theta), sin(theta); sin(theta), cos(theta)];
        VU(j,:,:) = expm(-1i * dt / 2 * H);
        Hs(j,:,:) = H;
    end
    % Initial wavefunction -- Gaussian wavepacket
    psi0 = zeros(M,M,N);
    psi0(:,1) = c1 * exp(1i*(kxI*x0)) .* exp(-(x0-xI).^2/sigmax^2);
    psi0(:,2) = c2 * exp(1i*(kxI*x0)) .* exp(-(x0-xI).^2/sigmax^2);
    psi0 = psi0 / sqrt(sum(sum(sum(abs(psi0).^2))));
    % psim & psi_k_m -- for plot
    psi0_k(:,1) = fftshift(fft(psi0(:,1)));
    psi0_k(:,2) = fftshift(fft(psi0(:,2)));
    psi0_k = psi0_k / sqrt(sum(sum(abs(psi0_k).^2)));
    psi_max = max(max(max(abs(psi0).^2)));
    psi_k_max = max(max(max(abs(psi0_k).^2)));
    % propagate WF
    psi = psi0;
    for t=0:Nstep-1
        % exp(-iVdt/2) * |Psi> in diab repr
        psi_k(:,1) = VU(:,1,1).*psi(:,1) + VU(:,1,2).*psi(:,2);
        psi_k(:,2) = VU(:,2,1).*psi(:,1) + VU(:,2,2).*psi(:,2);
        % exp(-iTdt) * psi
        psi_k(:,1) = TU .* fft(psi_k(:,1));
        psi_k(:,2) = TU .* fft(psi_k(:,2));
        % exp(-iVdt/2) * psi
        psi_k(:,1) = ifft(psi_k(:,1));
        psi_k(:,2) = ifft(psi_k(:,2));
        psi(:,1) = VU(:,1,1).*psi_k(:,1) + VU(:,1,2).*psi_k(:,2);
        psi(:,2) = VU(:,2,1).*psi_k(:,1) + VU(:,2,2).*psi_k(:,2);
        % analysis & report
        if mod(t,tgraph) == 0
            % analysis in diab
            psi_k(:,1) = fftshift(fft(psi(:,1)));
            psi_k(:,2) = fftshift(fft(psi(:,2)));
            psi_k = psi_k / sqrt(sum(sum(abs(psi_k).^2)));
            
            norm_k1 = sum(abs(psi_k(:,1)).^2) + 1e-16;
            norm_k2 = sum(abs(psi_k(:,2)).^2) + 1e-16;

            p1x = sum(abs(psi_k(:,1)).^2 .* kx0) / norm_k1;
            p2x = sum(abs(psi_k(:,2)).^2 .* kx0) / norm_k2;

            KE = sum((abs(psi_k(:,1)).^2 + abs(psi_k(:,2)).^2) .* (kx0.^2) / 2 / mass);
            PE = sum( conj(psi(:,1)) .* Hs(:,1,1) .* psi(:,1) + conj(psi(:,1)) .* Hs(:,1,2) .* psi(:,2) + conj(psi(:,2)) .* Hs(:,2,1) .* psi(:,1) + conj(psi(:,2)) .* Hs(:,2,2) .* psi(:,2) );

            % output
            if t == 0
                fprintf(fID, '# EXACT\n');
                fprintf(fID, '# xI = %8.4f kxI = %8.4f sigmax = %8.4f A = %8.4f init_s = %8.4f c1 = %8.4f c2 = %8.4f \n', ...
                                    xI, kxI, sigmax, A, init_s, c1, c2);
                fprintf(fID, '# L = %8.4f M = %8d dt = %8.4f Nstep = %8d tgraph = %8d\n', ...
                                    L, M, dt, Nstep, tgraph);
                fprintf(fID, '#%9s%16s%16s%16s%16s%16s\n', ...
                                't', 'n0', 'n1', 'px0', 'px1', 'Etot');
            end

            fprintf(fID, '#%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f\n', ...
                        t*dt, ...
                        sum(sum(abs(psi(:,1).^2))), ...
                        sum(sum(abs(psi(:,2).^2))), ...
                        p1x, p2x, ...
                        KE+PE ...
                        );
            % plot
            if enable_plot == true
                subplot(1,2,1);
                plot(x0,abs(psi(:,1)).^2, '-b', ...
                     x0,abs(psi(:,2)).^2, '--r');
                ylim([0, psi_max]);
                title('Real space Diab');

                subplot(1,2,2);
                plot(kx0,abs(psi_k(:,1)).^2, '-b', ...
                     kx0,abs(psi_k(:,2)).^2, '--r');
                ylim([0, psi_k_max]);
                title('Mom space Diab');

                drawnow;
            end

        end
    end
    fprintf(fID, '%16.10f%16.10f%16.10f%16.10f%16.10f\n', ...
            kxI, ...
            sum(sum(abs(psi(:,1).^2))), ...
            sum(sum(abs(psi(:,2).^2))), ...
            p1x, p2x ...
            );
    fprintf('# '); toc;
end
