function [x] = OpenMPI_kz( A,b,lambd, iter)
% regularized Kaczmarz
% As published here: http://iopscience.iop.org/article/10.1088/0031-9155/55/6/003/meta on page 1582.
% Other references : Saunders M 1995: Solution of sparse rectangular
% systems using LSQR and CRAIG
% or
% From Dax A 1993: On row relaxation methods for large constrained 
% least squares problems

voxel = 19;

shuff = 0;
enforceReal = 1;
enforcePositive = 1;

iterations = iter;

% initialization of the variable
[N, M] = size(A);  % N:column  M:row 
x = complex(zeros(N,1)); 
residual = complex(zeros(M,1));
rowIndexCycle = 1:M;

x_prev=zeros(N,1);
energy = rowEnergy(A);

if shuff
    rowIndexCycle = randperm(M);
end

lambdZero = sum(energy.^2)/N;
l = 1;
lambdIter = lambd*lambdZero;

while l <= iterations
%     fprintf([num2str(l),'\n'])
    for m = 1:M
        k = rowIndexCycle(m); 
        
        if energy(k) > 0
            tmp = A(:,k).'*x;
             beta = (b(k) - tmp - sqrt(lambdIter)*residual(k)) / (energy(k)^2 + lambdIter);
            x = x + beta*conj(A(:,k));
            residual(k) = residual(k) + beta*sqrt(lambdIter);
        end
    end
    
    if enforceReal && ~isreal(x)
        x = complex(real(x),0);
    end

    if enforcePositive
        x(real(x) < 0) = 0;
    end
       
    l = l + 1;

end
x=reshape(x, voxel, voxel,voxel);
end

