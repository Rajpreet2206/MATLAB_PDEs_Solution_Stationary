%% Function to Compute Matrix A
% The discretized equation for the Heat Equation is given by:
% [T(i+1,j) -2*T(i,j) + T(i-1,j)]/hx^2 + [T(i,j+1) -2*T(i,j) + T(i,j-1)]/hy^2 = f(xi,yi)
% where i,j denote the x and y nodes respectively
%%

% INPUTS: 
% Nx: Number of x nodes with unknown temperatures
% Ny: Number of y nodes with unknown temperatures
% OUTPUTS:
% T: Coefficient Matrix (symbolic)
%%
function [T,x,b] = matrixA_sym(Nx,Ny)
    % Initialization
    N = Nx*Ny;
    T = sym(zeros(N:N));
%     b = sym(zeros(N,1));
%     x = sym(zeros(size(b)));
    
    % Create vector of all indexes [11;12;13;21;22;23...]
    i = 1:Nx; j = 1:Ny;
    [I,J] = meshgrid(i,j); c=cat(2,J',I'); 
    index=reshape(c,[],2);
    
    % Heat Equation Loop
    %% A Matrix
    for iter = 1:N
        row = iter;
        col = iter;
        
        % DIAGONAL ELEMENT
        T(row,col) = -2/sym('hx')^2 -2/sym('hy')^2; 
        
        % ADJACENT LEFT (coeff. for i,j-1 node)
        if( index(iter,2) - 1 > 0)  % in ij, j-1 should be >0
            T(row, col-1) = 1/sym('hy')^2;
        end
        
        % ADJACENT RIGHT COEFFICIENT (for i,j+1 node)
        if( index(iter,2) + 1 <= Ny)  % in ij, j+1 should be <Ny
            T(row, col+1) = 1/sym('hy')^2;
        end
                
        % Ny ELEMENTS TO THE RIGHT (Coeff. for i+1,j node)
        if (col+Ny)<=N
            T(row,col+Ny) = 1/sym('hx')^2;
        end
        
        % Ny ELEMENTS TO THE RIGHT (Coeff. for i+1,j node)
        if (col-Ny)>0
            T(row,col-Ny) = 1/sym('hx')^2;
        end      
    end
    
    %% b && x Matrices
    b = sym("f_", [Nx,Ny]);
    x = sym("T_", [Nx,Ny]);
    b = reshape(b,[],1);
    x = reshape(x,[],1);
    
end