% Solving At = b

%% INPUTS

% f:  Function handle for f(x,y) = RHS of PDE
% Nx: Number of x nodes with unknown temperatures
% Ny: Number of y nodes with unknown temperatures

%% OUTPUT
% T : Nx by Ny matrix of temperature value at the nodes
% storage: size rquirement of this method
%%
function [T, storage] = full_matrix_solver(f, Nx, Ny,P)
    %fprintf('FULL MATRIX SOLVER: ');
    % Initialisation
    
    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    N = Nx*Ny;
    t = zeros(N,1); % vector of unknown temperatures [At=b]
    A = zeros(N:N); % system matrix
    b = zeros(N,1); % vector for f(x,y) at interior nodes
    b_known_terms = zeros(N,1); % values to be added to f(x,y) in b
    x = linspace(hx,1-hx,Nx); % x coordinates of interior points
    y = linspace(hy,1-hy,Ny); % y coordinates of interior points
    T = zeros(Ny,Nx); % Matrix of node temperatures
    
    % populate the system matrix
    [A,b_known_terms] = matrixA(Nx,Ny);
    %spy(A);//to view the sparsity of the matrix- uncomment this
    
    % populate b vector
    for i = 1 : Nx
        for j = 1 : Ny
            iter = (i-1)*Ny + j;
            b(iter,1) = f(x(i),y(j));
        end
    end
    b = b + b_known_terms; % b_known terms will be all zeros because boundary condition is zero
    
    % Solve for t using matlab backslash operator
    t = A\b;
    
    % Populate T matrix from t vector (excluding boundary points)
    for i = 1 : Nx
        for j = 1 : Ny
            iter = (i-1)*Ny + j;
            T(j,i) = t(iter,1);
        end
    end
    
%% This section doesn't run when we call it thru the Timeit function.
    
    if(P==1) 
    % Calculating Storage Requirements
    
        storage = numel(A)+numel(t)+numel(b)+numel(b_known_terms)+numel(T) +numel(x)+numel(y)+ numel(Nx)+numel(Ny)+...
                  +numel(hx)+numel(hy)+numel(i)+numel(j);
    
    % Add zero-boundaries and plotting
        T_bound = [zeros(1,Nx); T; zeros(1,Nx)];
        T_bound = [zeros(Ny+2,1), T_bound, zeros(Ny+2,1)];
        title = strcat('Full matrix solver for Nx = Ny = ', num2str(Nx));
          
        plotter(T_bound,Nx,Ny,title);
    end
    
end
