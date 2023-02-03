%% INPUTS
% T_noB : temperature values on NxbyNy nodes (excluding boundaries)
% Nx : no of grid points along x axis
% Ny : no of grid points along y axis
% f_anal : analytical solution

%% OUTPUTS
% error: error value of T_noB against f_anal

%%
function error_ = error_numeric(T_noB, Nx, Ny, f_anal)
    E = zeros(Nx,Ny);
    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    x = linspace(hx,1-hx,Nx); % x coordinates of interior points
    y = linspace(hy,1-hy,Ny); % y coordinates of interior points
    % Calculating Error
    for i = 1:Nx
        for j = 1:Ny
            E(i,j) = T_noB(i,j) - f_anal(x(i),y(j));
        end
    end

    error_ = sqrt(sum(sum((E.^2))/(Nx*Ny)));
    % error = sqrt((1/(Nx*Ny))*sumsqr(E(i,j)));
end