classdef he_atom_class
    %HE_ATOM_CLASS does calculate eigenenegy for the 1s2 state of He atom using pseudospectral method
    % Written by Tsogbayar Tsednee (PhD), California State University Northridge;
    % Contact: tsog215@gmail.com
    % Reference: Ts.Tsogbayar & D. L. Yeager,  Chinese Physics B 26, 083101 (2017)
    % March 18, 2020

    properties
    end
    
    methods (Static)
        
        function [r, psi, Vee, En1, En_total] = he_1s2_state(Z,N,a,b,itermax,tol)
            %
            [r,wr,D] = he_atom_class.legDC2(N,a,b);
            %            
            wr = wr(2:N);
            r = r(2:N);
            D2 = (2/(b-a))^2*D^2;
            D2 = D2(2:N,2:N);
            %%%
            wf_1s = 2*Z^(3/2).*exp(-Z*r)/sqrt(4*pi); % 1s state
            %
            [Vee, Vpois] = he_atom_class.Poisson_solver(D2, wr, r, b, a, wf_1s, N) ;
            %
            [psi, En1,En_total] = he_atom_class.Hartree_Fock_solver(D2, wr, r, Z, Vee, Vpois, N);
            %            
            psi_old = psi;
            %
%            Iflag = 0.;  % converfence flag, 1 indicates convergence
            %
            for iter = 2:itermax
                
                iter;
                cormax = 0.;
                psi = psi_old;
                wf_1s = psi./sqrt(4.*pi); % radial wave function
                [Vee, Vpois] = he_atom_class.Poisson_solver(D2, wr, r, b, a, wf_1s, N) ;
                %
                [psi, En1, En_total] = he_atom_class.Hartree_Fock_solver(D2, wr, r, Z, Vee, Vpois, N);
                %
                psi_new = psi;
                %
                cor = rms(psi_new - psi_old);
                if (abs(cor) > cormax);
                    cormax = cor;
                end
                %
                psi_old = psi_new;
                %
                if (cormax < tol)
                    %    if (cor < tol)
                    Iflag = 1.
                    break
                end
                %
            end
        %   
        end
        %%%
        function [Vee, Vpois] = Poisson_solver(D2, wr, r, b, a, wf_1s, N)
            %
            f_rhs = -4.*pi*r.*conj(wf_1s).*wf_1s; % w/o 4*pi for radial wf
            u = D2\f_rhs;
            u = u + 1.*r/(b-a); % for 1s state or s-states
            Vpois = u./r;
            %
            sm = 0.;
            for i = 1:N-1
                sm = sm + wr(i)*wf_1s(i)*wf_1s(i)*Vpois(i)*4*pi*r(i)*r(i) ;    % (u(i)/r(i)) * r(i)^2 = u(i) * r(i)    J(1s^2|1s^2) = (5/8)*Z
            end
            %sm
            Vee = sm;
            %
            return
        end
        %%%
        %
        function [psi, En1, En_total] = Hartree_Fock_solver(D2, wr, r, Z, Vee, Vpois, N)
            %%% Hartee-Fock Hamiltonian
            H_ham = -0.5*D2 - diag((Z./r)) + diag(Vpois);
            [Vec,En] = eig(H_ham);                                     % Eigenvalue problem
            En = diag(En);
            [foo, ij] = sort(En);
            En = En(ij);
%            [En(1),En(2),En(3),En(4),En(5)];
            
            Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
            V1 = Vec(:,1);                         % The unnormalized eigenfunction for the ground state,
            % here 1 is for the ground state, 2,3... are corresponded
            % for the excited states
            %V1 = [0.,;V1,;0.];
            sm = 0.;
            for i = 1:N-1
                sm = sm + wr(i)*V1(i)*V1(i);
            end
            n_c = 1./sqrt(abs(sm));      % The normalization constant
            V1 = n_c*V1;
            n_wf = 0.;
            for i = 1:N-1
                n_wf = n_wf + wr(i)*V1(i)*V1(i);
            end
            psi = V1./r;
            %
            En1 = En(1);
            En_total = 2*En1 - Vee ;
            return
        end

        
        
        function [xi,w,D]=legDC2(N,a,b)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % legDc.m
            %
            % Computes the Legendre differentiation matrix with collocation at the
            % Legendre-Gauss-Lobatto nodes.
            %
            % Reference:
            %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
            %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
            %
            % Written by Greg von Winckel - 05/26/2004
            % Contact: gregvw@chtm.unm.edu
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % Truncation + 1
            N1=N+1;
            
            % CGL nodes
            xc=cos(pi*(0:N)/N)';
            
            % Uniform nodes
            xu=linspace(-1,1,N1)';
            
            % Make a close first guess to reduce iterations
            if N<3
                x=xc;
            else
                x=xc+sin(pi*xu)./(4*N);
            end
            
            % The Legendre Vandermonde Matrix
            P=zeros(N1,N1);
            
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            
            xold=2;
            while max(abs(x-xold))>eps
                
                xold=x;
                
                P(:,1)=1;    P(:,2)=x;
                
                for k=2:N
                    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
                end
                
                x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
            end
            
            X=repmat(x,1,N1);
            Xdiff=X-X'+eye(N1);
            
            L=repmat(P(:,N1),1,N1);
            L(1:(N1+1):N1*N1)=1;
            D=(L./(Xdiff.*L'));
            D(1:(N1+1):N1*N1)=0;
            D(1)=(N1*N)/4;
            D(N1*N1)=-(N1*N)/4;
            
            % Linear map from[-1,1] to [a,b]
            xi=(a*(1-x)+b*(1+x))/2;        % added by Tsogbayar Tsednee
            
            % Compute the weights
            w=(b-a)./(N*N1*P(:,N1).^2);    % added by Tsogbayar Tsednee
            
        end
        
    end
    
end

