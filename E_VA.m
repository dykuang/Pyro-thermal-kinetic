function [E, obj] = E_VA(csv_file,index)
% This function calculates the activation energy during the
% reaction at  each conversion rate by using VA's method.
%
% inputs:
% csv_file: name of file
% index: subgroup of data used, row vector 
%
% Outputs:
% E: a col vector containg calculated activation energy.
% obj: value of obj function used in VA, can be used to access
% roughly how "well" the result is. A col vector

%% parameter setting
par = set_parameter();
R = par.R;
beta = par.beta; % heating rates
n = 5; % power for I
batches = size(index,2);

%% importdata
data = importdata(csv_file); 
% data = data.data;
sample_density = size(data,1);
E = zeros(sample_density, 1);
obj = E;

Delta_T = [data(1, :); data(2:end, :) - data(1:end-1, :)];

[E(1), obj(1)] = fminbnd(@(x) ...
    ComputeOBJ0(x, Delta_T(1,index), beta(index)), 10, 1e6);

for i = 2:sample_density    
    [E(i), obj(i)] = fminbnd(@(x) ComputeOBJn(x, ...
        data(i,index), Delta_T(i,index), beta(index)),10, 1e6);
end

% E = E/1000;


%% visualization of the result, comment them out
% figure(1)
% plot(data(:,1), E, '-^');
% 
% figure(2)
% plot(data(:,1), obj);

    
%% compute object functions for minimization purpose.
function OBJ0 = ComputeOBJ0(E, dT, beta)
    OBJ0=0;

    I = zeros(1,batches);
    for s = 1 : batches
        I(s) = integral(@(T)exp(-E./T/R),0,dT(s));
    end
    J = (I./beta).^n;
%     J = log(I./beta);
%     J = I./beta;
    for s = 1: batches-1
        for j = s+1:batches
            OBJ0 = OBJ0 + J(s)/J(j) + J(j)/J(s);
%             OBJ0 = OBJ0 + (J(s) - J(j))^2;
        end
    end    

end


function OBJn = ComputeOBJn(E, T, dT, beta)
    OBJn=0;

    I = zeros(1,batches);
    for s = 1 : batches
        I(s) = integral(@(T)exp(-E./T/R),T(s)-dT(s),T(s));
    end
    J = (I./beta).^n;
%     J = log(I./beta);
%     J = I./beta;
    for s = 1: batches-1
        for j = s+1:batches
            OBJn = OBJn + J(s)/J(j) + J(j)/J(s);
%             OBJn = OBJn + (J(s) - J(j))^2;
        end
    end    

end

end

