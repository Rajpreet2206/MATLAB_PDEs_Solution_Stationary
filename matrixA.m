%% Function to Compute Matrix A and vector of known_terms for Ax=b system as below:
% The discretized equation for the Heat Equation is given by:
% [T(i+1,j) -2*T(i,j) + T(i-1,j)]/hx^2 + [T(i,j+1) -2*T(i,j) + T(i,j-1)]/hy^2 = f(xi,yi)
% where i,j denote the node at x=i.hx and y=j.hy
% This system can be written as a matrix system:
% Ax = b
% where x carries the unknowns, b carries the knowns
%%

% INPUTS: 
% Nx: Number of x nodes with unknown temperatures
% Ny: Number of y nodes with unknown temperatures
% OUTPUTS:
% A: Coefficient Matrix
% known_terms: vector of known temperatures (at boundary)
%%
function [A,known_terms] = matrixA(Nx,Ny)
    % Initialization
    N = Nx*Ny;
    A = zeros(N:N); % system matrix
    known_terms = zeros(N,1); % vector for knowns to be added to b in Ax=b
    hx = 1/(Nx+1); % assuming x domain is 0 to 1
    hy = 1/(Ny+1); % assuming y domain is 0 to 1
%     x = linspace(hx,1-hx,Nx); % x coordinates of interior points
%     y = linspace(hy,1-hy,Ny); % y coordinates of interior points
    
    % Create vector of all indexes [11;12;13;21;22;23...]
    i = 1:Nx; j = 1:Ny;
    [I,J] = meshgrid(i,j); c=cat(2,J',I'); 
    index=reshape(c,[],2);
    
    % Populate the first row
    for iter = 1:N  
        row = iter;
        col = iter;
        
        % DIAGONAL ELEMENT
        A(row,col) = -2/(hx^2) -2/(hy^2); 
        
        % ADJACENT LEFT (coeff. for i,j-1 node)
        if( index(iter,2) - 1 > 0)  % in ij, j-1 should be >0
            A(row, col-1) = 1/(hy^2);
        else
            known_terms(row) = known_terms(row)+0;
            % This section implements boundary conditions
        end
        
        % ADJACENT RIGHT COEFFICIENT (for i,j+1 node)
        if( index(iter,2) + 1 <= Ny)  % in ij, j+1 should be <Ny
            A(row, col+1) = 1/(hy^2);
        else
            known_terms(row) = known_terms(row)+0;
            % This section implements boundary conditions
        end
                
        % Ny ELEMENTS TO THE RIGHT (Coeff. for i+1,j node)
        if (col+Ny) <= N
            A(row,col+Ny) = 1/(hx^2);
        else
            known_terms(row) = known_terms(row)+0;
            
        end
        
        % Ny ELEMENTS TO THE LEFT (Coeff. for i-1,j node
        if (col-Ny) > 0
            A(row,col-Ny) = 1/(hx^2);
        else
            known_terms(row) = known_terms(row)+0;

        end
    end  
end