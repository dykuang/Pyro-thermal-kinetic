function [E, obj] = E_VO(csv_file,index)
% This function calculates the activation energy during the
% reaction at each conversion rate by using VO's method.
%
% Inputs:
% csv_file: csv file to be loaded. row alpha, column beta.
% index: which columns to use, row vector 
%
% Outputs:
% E: a col vector containg calculated activation energy.
% obj: value of obj function used in VA, can be used to access
% roughly how "well" the result is. A col vector

%% parameter setting
par = set_parameter();
R = par.R;
beta = par.beta; % heating rates

batches = size(index,2);

%% importdata
data = importdata(csv_file); 
% data = data.data;
% data = data(:,2:end);
sample_density = size(data,1);
E = zeros(sample_density, 1);
obj = E;


[E(1), obj(1)] = fminbnd(@(x) ...
    ComputeOBJ0(x, data(1,index), beta(index)), 10, 1e6);

for i = 2:sample_density    
    [E(i), obj(i)] = fminbnd(@(x) ComputeOBJn(x, ...
        data(i,index), beta(index)),10, 1e6);
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
    J = I./beta;
    
    for s = 1: batches-1
        for j = s+1:batches
            OBJ0 = OBJ0 + J(s)/J(j) + J(j)/J(s);
        end
    end    

end


function OBJn = ComputeOBJn(E, T, beta)
    OBJn=0;

    I = zeros(1,batches);
    for s = 1 : batches
        I(s) = integral(@(T)exp(-E./T/R),0,T(s));
    end
    J = I./beta;

    for s = 1: batches-1
        for j = s+1:batches
            OBJn = OBJn + J(s)/J(j) + J(j)/J(s);
        end
    end    

end

end

