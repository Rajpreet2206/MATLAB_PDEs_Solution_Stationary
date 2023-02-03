%% Worksheet 4: Solving Stationary Partial Differential Equations

%% General Instructions 
% 1. After computing results for one grid size, the program waits for user prompt. At that 
%point, press any key to continue. 

%% Defining Constants and Functions

tol = 1e-4;
f = @(x,y) -2*pi^2*sin(pi*x)*sin(pi*y);
f_anal = @(x,y) sin(pi.*x).*sin(pi.*y);

%% Creating the Full matrix A for the linear system

[A_sym, x_sym, b_sym] = matrixA_sym(5,5);
disp('Coefficient Matrix A');
disp(A_sym);
disp('Approximate Temperature at Grid Points T');
disp(x_sym);
disp('Resultant Temperature Function x');
disp(b_sym);
spy(A_sym); % to visualise the sparsity of A

%% Solving the Analytical Solution
% T-matrix from analytical solution

T_anal = analytical_solver(f_anal,15,15);

%% Solving Full Matrix, Sparse Matrix and Gauss Seidel for various Grid Sizes

iter = 0;
N_array = [3,7,15,31,63,127];
for Grid_size = N_array
     
     pause;
     disp('Calculating results for grid size below: ');
     disp(Grid_size);
     disp('Computations in Smaller Grid Sizes take more time, press Any key after you see countour and surface plots');
     
      iter = iter+1;
      Nx = Grid_size; Ny = Grid_size; 
      
     if (Grid_size < 127) % we want matrix solvers to run only till 63x63 nodes
         
          % Solving using NxN full system matrix
          [T_d1, storage(1,iter)] = full_matrix_solver(f, Nx, Ny, 1);
          time(1,iter)= timeit(@()full_matrix_solver(f, Nx, Ny, 0));

          % Solving using NxN sparse system matrix
          [T_d2, storage(2,iter)] = sparse_matrix_solver(f, Nx, Ny, 1);
          time(2,iter)= timeit(@()sparse_matrix_solver(f, Nx, Ny, 0));
      
     end

     % Solving matrix using Gauss-Seidel
      [T_d3, storage(3,iter)] = gaussSeidal(f, Nx, Ny, tol, 1); % T_d3 is temperature matrix excluding the boundaries
      time(3,iter)= timeit(@()gaussSeidal(f, Nx, Ny, tol, 0));
      Err(iter) = error_numeric(T_d3, Nx, Ny, f_anal); % error calculation
      
      err_red(1) = 0;
      if(iter>1)
          err_red(iter) = Err(iter-1)/Err(iter);
      end        
end


%% Table for Runtime and Storage for Full Matrix Solver

RS_Table(time(1,(2:5)), storage(1,(2:5)), N_array(2:5), "Grid Size", "Runtime and Storage for Full Matrix Solver");

%% Table for Runtime and Storage for Sparse Matrix Solver

RS_Table(time(2,(2:5)), storage(2,(2:5)), N_array(2:5), "Grid Size", "Runtime and Storage for Sparse Matrix Solver");

%% Table for Runtime and Storage for Gauss Seidal Solver

RS_Table(time(3,(2:5)), storage(3,(2:5)), N_array(2:5), "Grid Size", "Runtime and Storage for Gauss Seidel Solver");

%% Table for Errors and Reduced Error in Gauss Seidal

Error_table(Err(2:6), err_red(2:6), N_array(2:6), "Grid Size", "Error & Error Reduced for Gauss Seidel");






