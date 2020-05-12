function [] = DFT_Poisson_solver_H2_ground_state
%-----------------------------------------------------------------
% Solves the Hartree-Fock-Slater (Xalpha) equation for ground state of molecular
% hydrogen using pseudospectral method.
% The electron-electron interaction Vee is obtained as solution of the
% Poisson equation. 
% Coordinate system: Prolate spheroidal coordinates (mu, nu)
% Uses: function legDC2 - differentation matrix based on
% Legendre-Gauss-Lobatto nodes
%
% Written by Tsogbayar Tsednee (PhD), California State University Northridge 
% Contact: tsog215@gmail.com
% Reference: Tsogbayar Tsednee, PhD thesis, available at https://yorkspace.library.yorku.ca/xmlui/handle/10315/28160
% Date: May 1, 2013
%%%
% ----------------------------------------------------------------
clear;
clc;
format long
N = 64.; M = 12; a = 1.; b = 50.;  % Initial computaional parameters; you may change them
Rint = 1.400; % internuclear separation
itermax = 100; tol = 10^(-6);  % max number of iteration and tolerance; you may change them
%
%xalf = 2./3.
xalf = 0.6683 ; % x-alpha parameter in exchange correlation 
rle = -1.800; 
rleorb = -0.5;
%rleac = -0.890;
opts.p = 48.;

% mu part begins
[mu,wmu,D]=legDC2(N,a,b); mup = mu;
%[mup,wmup,Dp]=legDC2(N,a,b);
D2 = D*D;
Tmu = diag(mu.^2 - 1.)*(2/(b-a))^2*D2 + ...
      2.*diag(mu)*(2/(b-a))^1*D;
% mu part ends
mu = mu(2:N+1);

% nu part begins
[nu,wnu,Dn] = legDC2(M,-a,a); nup = nu;
D2n = Dn*Dn;
Tnu = diag(1. - nu.^2)*(2/(a-(-a)))^2*D2n - ...
      2.*diag(nu)*(2/(a-(-a)))*Dn;
% nu part ends
%%%
[mum,nun] = meshgrid(mu,nu); 
mum = mum(:);
nun = nun(:);
Smunu= sparse(diag(mum.^2 - nun.^2));
%  
Tmunu = (4./Rint^2)*(kron(Tmu(2:N+1,2:N+1),eye(M+1)) + ...
                     kron(eye(N),Tnu(1:M+1,1:M+1)));
%               
Vcmunu = -(4./Rint)*diag(mum);
%
Hmn2 = sparse(-0.5*Tmunu) + sparse(Vcmunu); % two-center hamiltonian
  
%%%% eigenvalue problem 
% Lanczos begins
[Vec,Ener] = eigs ( Hmn2,Smunu, 3, rle, opts );
En = diag(Ener);
[foo, ij] = sort(En);
Eni = En(ij);
[Eni(1),Eni(2),Eni(3)] ;
% Lanczos ends

Vec = Vec(:,ij);  % The unnormalized eigenfunctions
V1 = (Vec(:,1));    % The unnormalized eigenfunction for the ground state, 
%
V1 = reshape(V1,M+1,N)';
Vrow1 = zeros(M+1,1)' ;
V1 = [Vrow1; V1];  
%figure(1)
%mesh(V1)
sm = 0.;
for i = 1:N+1
    for j = 1:M+1
        sm = sm + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                  (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;        
    end
end
sm; 
n_c = 1./sqrt(sm);      % The normalization constant
V1 = n_c*V1;
%figure(2)
%mesh(mup,nu,V1')
%%%
%%%
n_wf = 0.;
for i = 1:N+1
    for j = 1:M+1
%        V1(i,j) = n_c*V1(i,j);
        n_wf = n_wf + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                      (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;
    end
end
n_wf;
%
rho1stwo = 2.0*conj(V1).*V1;    % rho = 2|psi|^2
%
%%% ----- <|Vex|> electron-electron interaction -------
%%%
smexch = 0.;
for i = 1:N+1
    for j = 1:M+1
        smexch = smexch + wmu(i)*wnu(j)*...
                      (3./2.)*xalf*...
                      (-3./4.)*(3./pi).^(1./3.)*rho1stwo(i,j)^(4./3.)*...
                      (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;
    end
end
Vexch = smexch;
%%% ----- <|Vex|> electron-electron interaction -------

%%%
%%% initial rho0 begins
rho0 = 2.*conj(V1).*V1;
rho0 = reshape(rho0',(N+1)*(M+1),1);
%%% initial rho0 ends

mup(1);
%[nunp,mump] = meshgrid(nup,mup);
[mump,nunp] = meshgrid(mup,nup); 
mump = mump(:);
nunp = nunp(:);
Smunup = diag(mump.^2 - nunp.^2);
%%%
% Impose boundary conditions by replacing appropraite rows:
bmul = find(abs(mump) == mup(1) );
%
Hmn_p = (4./Rint^2)*(kron(Tmu(1:N+1,1:N+1),eye(M+1)) + ...
                     kron(eye(N+1),Tnu(1:M+1,1:M+1))); 

Hmn_p(bmul,:) = zeros(size(bmul,1),(N+1)*(M+1));
Hmn_p(bmul,bmul) = eye(size(bmul,1));

%%%  BC implementation begins
rhs = zeros((N+1)*(M+1),1);
%rhs(bmul) = (mump(bmul) == mup(1)).*(2.*1./(Rint.*mump(bmul)));
rhs(bmul) = (mump(bmul) == mup(1)).*(1./((Rint/2.).*sqrt(mump(bmul).^2 +...
                                                         nunp(bmul).^2 -...
                                                         1. ) ));
rhsrsh = reshape(rhs,M+1,N+1);
%%%  BC implementation ends

%%% Poisson solver begins    
%%% right hand side calculation begins
f = -Smunup*4*pi*rho0;
frsh = reshape(f,M+1,N+1);
%%% right hand side calculation ends
%frsh(1:Nn+1,1) = 0.;
ff = frsh + rhsrsh;
f = reshape(ff,(N+1)*(M+1),1);
%figure(1), clf, spy(Hmn_p), drawnow
tic, u = Hmn_p\f; toc   
uu = reshape(u,M+1,N+1);
%mesh(reshape(u,Nn+1,N+1))
%uup = zeros(Nn+1,N+1);
%uup(1:Nn+1,2:N+1) = uu;
%u = uu(1:Nn+1,2:N+1);
%u = reshape(u,(Nn+1)*N,1);
uup = uu';
%%% Poisson solver ends

%%% Electron-electron interaction calculation begins
%uv = reshape(u,Nn+1,N+1)
%%% initial rho0 begins           
rho0p = rho0;
%%% initial rho0 ends
rho0rshp = reshape(rho0p,M+1,N+1)';
%%%
Vee = 0.;
for i = 1:N+1
    for j = 1:M+1
        Vee = Vee + wmu(i)*wnu(j)*rho0rshp(i,j)*uup(i,j)*...
                    (Rint^3/8)*(mup(i)^2-nu(j)^2)*2*pi;
    end
end
Vee = 0.5*Vee ;
%%% Electron-electron interaction calculation ends

%E_tot = T_rho + V_ext + Vee + Vexch

%%%
%%% right hand side calculation ends
%%% H2+ matrix element calculation begins
    Tmunuh2 = (4./Rint^2)*(kron(Tmu(1:N+1,1:N+1),eye(M+1)) + ...
                         kron(eye(N+1),Tnu(1:M+1,1:M+1)));             
    Vcmunuh2 = -(4./Rint)*diag(mump);
    Hmnh2 = -0.5*Tmunuh2 + Vcmunuh2;
%Smunu= diag(mum.^2 - nun.^2);
%%% H2+ matrix element calculation ends

    Hmnh2(bmul,:) = zeros(size(bmul,1),(N+1)*(M+1));
    Hmnh2(bmul,bmul) = eye(size(bmul,1));

    %%%  BC implementation begins
    rhswf = zeros((N+1)*(M+1),1);
    %rhswf(bmul) = (mump(bmul) == mup(1)).*(2.*1./(Rint.*mump(bmul)));
    rhswf(bmul) = (mump(bmul) == mup(1)).*exp(-sqrt(-2.*Eni(1))*...    %%%Eni(1)
                                              (Rint/2.).*mump(bmul) );
    %rhsrshwf = reshape(rhswf,Nn+1,N+1);
    rhsrshwf = diag(rhswf);
%%%  BC implementation ends

    %%% right hand side calculation begins
    fh2 = Smunup;
    %fh2(bmul,:) = zeros(size(bmul,1),(N+1)*(Nn+1));
    %fh2(bmul,bmul) = eye(size(bmul,1));
    Smunuh2 = sparse(fh2 - rhsrshwf) ;

%%%

for iter = 2:itermax
%    
    Vee_old = Vee;
%
    iter % 2    
    
%%% --- Schrodinger 1 electron eigenvalue solver ---
    Vex_p = sparse(-(3./2.)*xalf*(3./pi).^(1./3.)*(2.*rho0).^(1./3.));
%
%%%
    Hmn = sparse(1*Hmnh2) + sparse(1.*diag(u)*Smunup) + sparse(1.*diag(Vex_p))*Smunup;

% Lanczos begins
    [Vec,Ener] = eigs ( Hmn,Smunuh2, 3, rleorb, opts );
    En = diag(Ener);
    [foo, ij] = sort(En);
    En = En(ij);
%    [En(1),En(2),En(3)] ;
% Lanczos ends
%    En_tot_mat = 2*En(1) - Vee - 0.5*Vexch + 1/Rint ;
    
    Vec = Vec(:,ij);  % The unnormalized eigenfunctions
    V1 = (Vec(:,1));    % The unnormalized eigenfunction for the ground state, 
                      % here 1 is for the ground state, 2,3... are corresponded 
                      % for the excited states

    V1 = reshape(V1,M+1,N+1)'; 
    %figure(1)
    %mesh(V1)
    sm = 0.;
    for i = 1:N+1
        for j = 1:M+1
            sm = sm + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                      (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;        
        end
    end
    sm; 
    n_c = 1./sqrt(sm);      % The normalization constant
    V1 = n_c*V1;
    %figure(2)
    %mesh(mup,nu,V1')
%%%
    %%%
    n_wf = 0.;
    for i = 1:N+1
        for j = 1:M+1
    %        V1(i,j) = n_c*V1(i,j);
            n_wf = n_wf + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                          (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;
        end
    end
%    n_wf;

%    psi = V1;
%    rho1s = 2.0*conj(psi).*psi;
    rho1stwo = 2.0*conj(V1).*V1;    % rho = 2|psi|^2

    %%% ----- <|Vex|> electron-electron interaction -------
    %%%
    smexch = 0.;
    for i = 1:N+1
        for j = 1:M+1
            smexch = smexch + wmu(i)*wnu(j)*...
                          (3./2.)*xalf*...
                          (-3./4.)*(3./pi).^(1./3.)*rho1stwo(i,j)^(4./3.)*...
                          (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;
        end
    end
    Vexch = smexch ;
    %%% ----- <|Vex|> electron-electron interaction -------

    %%%
    %%% initial rho0 begins
    rho0 = 2.*conj(V1).*V1;
    rho0 = reshape(rho0',(N+1)*(M+1),1);
    %%% initial rho0 ends

    mup(1);
    %[nunp,mump] = meshgrid(nup,mup);
    [mump,nunp] = meshgrid(mup,nup); 
    mump = mump(:);
    nunp = nunp(:);
    Smunup = diag(mump.^2 - nunp.^2);
    %%%
    % Impose boundary conditions by replacing appropraite rows:
    bmul = find(abs(mump) == mup(1) );
    %
    Hmn_p = (4./Rint^2)*(kron(Tmu(1:N+1,1:N+1),eye(M+1)) + ...
                         kron(eye(N+1),Tnu(1:M+1,1:M+1))); 

    Hmn_p(bmul,:) = zeros(size(bmul,1),(N+1)*(M+1));
    Hmn_p(bmul,bmul) = eye(size(bmul,1));

    %%%  BC implementation begins
    rhs = zeros((N+1)*(M+1),1);
    %rhs(bmul) = (mump(bmul) == mup(1)).*(2.*1./(Rint.*mump(bmul)));
    rhs(bmul) = (mump(bmul) == mup(1)).*(1./((Rint/2.).*sqrt(mump(bmul).^2 +...
                                                             nunp(bmul).^2 -...
                                                             1. ) ));
    rhsrsh = reshape(rhs,M+1,N+1);
    %%%  BC implementation ends

    %%% Poisson solver begins    
    %%% right hand side calculation begins
    f = -Smunup*4*pi*rho0;
    frsh = reshape(f,M+1,N+1);
    %%% right hand side calculation ends
    %frsh(1:Nn+1,1) = 0.;
    ff = frsh + rhsrsh;
    f = reshape(ff,(N+1)*(M+1),1);
    %figure(1), clf, spy(Hmn_p), drawnow
    tic, u = Hmn_p\f; toc   
    uu = reshape(u,M+1,N+1);
    %mesh(reshape(u,Nn+1,N+1))
    %uup = zeros(Nn+1,N+1);
    %uup(1:Nn+1,2:N+1) = uu;
    %u = uu(1:Nn+1,2:N+1);
    %u = reshape(u,(Nn+1)*N,1);
    uup = uu';
    %%% Poisson solver ends

    %%% Electron-electron interaction calculation begins
    %uv = reshape(u,Nn+1,N+1)
    %%% initial rho0 begins           
    rho0p = rho0;
    %%% initial rho0 ends
    rho0rshp = reshape(rho0p,M+1,N+1)';
    %%%
    Vee = 0.;
    for i = 1:N+1
        for j = 1:M+1
            Vee = Vee + wmu(i)*wnu(j)*rho0rshp(i,j)*uup(i,j)*...
                        (Rint^3/8)*(mup(i)^2-nu(j)^2)*2*pi;
        end
    end
    Vee = 0.5*Vee;
    %%% Electron-electron interaction calculation ends

%    rho0ds = conj(V1).*V1;
%%% wave function calculation ends

    if (abs(Vee - Vee_old) < tol)
        Iflag = 1.
        break
    end
%   
%    rho0d = rho0ds;    
    
    
end
%%% 
%
En_tot = 2*En(1) - Vee - 0.5*Vexch + 1./Rint ; % total energy for ground state

%%%% some simple calculation for the ground state ---
% n_wf: check normalization of wave function
% smrsq: calculation of < |r^2| >  at given R to check wave function
% smQ: calculation of < |Q| >  at given R to check wave function 
n_wf = 0.;  
smrsq = 0.;
smQ = 0.;
for i = 1:N+1
    for j = 1:M+1
%        V1(i,j) = n_c*V1(i,j);
        n_wf = n_wf + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                  (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;
%              
        smrsq = smrsq + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                    (Rint^2/4)*(mup(i).^2 + nu(j).^2 - 1)*...
                    (Rint^3/8)*(mup(i).^2-nu(j).^2)*2*pi;              

%%% <Q>---        
        smQ = smQ + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                        0.5*(3.*(Rint^2/4)*mup(i).^2*nu(j).^2 - ...
                            (Rint^2/4)*(mup(i).^2 + nu(j).^2 - 1))*...
                            (Rint^3/8)*(mup(i).^2-nu(j).^2)*2*pi;
%%% <Q>---                  
                
    end
end
%n_wf;
%smrsq;   % < |r^2| >  at given R to check wave function 
%smQ;     % < |Q| >  at given R to check wave function 

% Output: 
[Rint, En(1), Vee, En_tot, smrsq, smQ] % 
% 1.400000000000000  -0.434191137784125   1.309266870261820  -1.174332586993145   2.508839933575087   0.267700737327968 % 

%%% En_tot = -1.17433 vs Exact = -1.17443 from W. Kolos and C.C. Roothaan, Rev. Mod. Phys. 32, 219 (1960)


return
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
xi=(a*(1-x)+b*(1+x))/2;      % added by Tsogbayar Tsednee

% Compute the weights
w=(b-a)./(N*N1*P(:,N1).^2); % added by Tsogbayar Tsednee
return
end

