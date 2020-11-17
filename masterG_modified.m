%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script does masterplot method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pdata] = masterG_modified(E, file_T)

 par = set_parameter();
 R = par.R;
%  beta = par.beta;
 alpha = par.alpha;
 j = 2;
 
%  T = importdata(file_T);
 T = file_T;
 T = T(:,j);
% need to check if size of E matches with size of alpha

% Power law
p2 = @(x) x.^(1/2);
p3 = @(x) x.^(1/3);
p4 = @(x) x.^(1/4);
% p5 = @(x) x.^(1/5);
% p6 = @(x) x.^(1/6);
% p7 = @(x) x.^(1/7);

% Limited surface reation between both phases
% r1 = @(x) x;
% r2 = @(x) 1-(1-x).^(1/2);
% r3 = @(x) 1-(1-x).^(1/3);

% geometrical contraction model
r1 = @(x) x.^2;
r2 = @(x) (1-x).*log(1-x);
r3 = @(x) ( (1-x).^(-2)-1 )/2;

% Random
a2 = @(x) (-log(1-x)).^(1/2);
a3 = @(x) (-log(1-x)).^(1/3);
a4 = @(x) (-log(1-x)).^(1/4);
% a5 = @(x) (-log(1-x)).^(1/5);

% % Random -- fraction
% a17 = @(x) (-log(1-x)).^(1);
% a18 = @(x) (-log(1-x)).^(2);
% a19 = @(x) (-log(1-x)).^(3);
% a110 = @(x) (-log(1-x)).^(4);
% a111 = @(x) (-log(1-x)).^(5);

% Diffusion
d1 = @(x) x.^2;
d2 = @(x) x + (1-x).*log(1-x);
d3 = @(x) (1-(1-x).^(1/3)).^2;
d4 = @(x) (1-2*x/3)-(1-x).^(2/3);

% % Order of reaction
% f1 = @(x) -(log(1-x));
% f2 = @(x) 1./(1-x)-1;
% f3 = @(x) ((1-x).^(-2)-1)/2;
f2 = @(x) 1-(1-x).^(1/2);
f3 = @(x) 1-(1-x).^(1/3);
f4 = @(x) 1./(1-x);

P = zeros(size(T,1), 1);

% P(1)=integral(f, E(1)/R./T(1,5), inf, 'AbsTol',1e-6,'RelTol',1e-3);
% for i = 2 : size(T,1)
%     P(i) = P(i-1)+...
%         integral(f, E(i)/R./T(i,5), E(i-1)/R./T(i-1,5), 'AbsTol',1e-6,'RelTol',1e-3); 
% end
% flag = (eye(size(E,1),size(E,1))*E)./R./T;


% This representation gives smoother result
P(1)=integral(@(x) exp(-E(1)/R./x), 0, T(1), 'AbsTol',1e-6,'RelTol',1e-3);
for i = 2 : size(T,1)
    P(i) = P(i-1)+ ...
           integral(@(x) exp((-E(i-1)-E(i))/R./x/2), T(i-1), T(i),'AbsTol',1e-6,'RelTol',1e-3);
%            integral(@(x) exp((-E(i-1))/R./x), T(i-1), T(i), 'AbsTol',1e-6,'RelTol',1e-3);
end
% 
% for i = 1: size(T,1)
%     P(i) = 0.00484*exp(-1.0516*E(i)/R/T(i));
% end

P = P/P(alpha==0.5);

% figure(1)
% plot(alpha, p2(alpha)/p2(0.5),alpha, p3(alpha)/p3(0.5),alpha, p4(alpha)/p4(0.5),...
% alpha, p5(alpha)/p5(0.5),alpha, p6(alpha)/p6(0.5), alpha, p7(alpha)/p7(0.5), alpha, P, 'o-')
% legend('P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'experiment')
% 
% % hold on
% figure(2)
% plot(alpha, r1(alpha)/r1(0.5),alpha, r2(alpha)/r2(0.5),alpha, r3(alpha)/r3(0.5),...
%     alpha, P, 'o-');
% legend('r1', 'r2', 'r3', 'experiment')
% 
% figure(3)
% plot(alpha, a2(alpha)/a2(0.5),alpha, a3(alpha)/a3(0.5), ...
%     alpha, a4(alpha)/a4(0.5),alpha, a5(alpha)/a5(0.5), alpha, P,'o-');
% legend('a2', 'a3', 'a4', 'a5', 'experiment')
% 
% figure(4)
% plot(alpha, a17(alpha)/a17(0.5),alpha, a18(alpha)/a18(0.5), ...
%     alpha, a19(alpha)/a19(0.5),alpha, a110(alpha)/a110(0.5),...
%     alpha, a111(alpha)/a111(0.5),alpha, P,'o-');
% legend('a12', 'a13', 'a14', 'a15', 'a16','experiment')
% 
% figure(5)
% plot(alpha, d1(alpha)/d1(0.5),alpha, d2(alpha)/d2(0.5), ...
%     alpha, d3(alpha)/d3(0.5),alpha, d4(alpha)/d4(0.5), alpha, P,'o-');
% legend('d1', 'd2', 'd3', 'd4', 'experiment')
% 
% figure(6)
% plot(alpha, f1(alpha)/f1(0.5),alpha, f2(alpha)/f2(0.5), ...
%     alpha, f3(alpha)/f3(0.5), alpha, P,'o-');
% legend('f1', 'f2', 'f3', 'experiment')
% hold off
styles = {'y', 'm', 'c', 'r', 'g', 'b', 'k', ...
          'ks--', 'mo--', 'c*-', 'r^-', 'g<:', 'b>--', 'kh:',...
          'm+:', 'g*-'};
funcs = {p2, p3, p4, r1, r2, r3, a2, a3, a4, d1, d2, d3, d4, f2, f3, f4};
labels = {'P2'; 'P3'; 'P4'; 'R1'; 'R2'; 'R3';...
         'A2'; 'A3'; 'A4';...
         'D1'; 'D2'; 'D3'; 'D4';...
         'F2'; 'F3'; 'F4'; 'experiment'};
num_funcs = size(funcs);
% plot(alpha, p2(alpha)/p2(0.5),alpha, p3(alpha)/p3(0.5),alpha, p4(alpha)/p4(0.5),...
%      alpha, r1(alpha)/r1(0.5),alpha, r2(alpha)/r2(0.5),alpha, r3(alpha)/r3(0.5),...
%      alpha, a2(alpha)/a2(0.5),alpha, a3(alpha)/a3(0.5), alpha, a4(alpha)/a4(0.5), ...
%      alpha, d1(alpha)/d1(0.5),alpha, d2(alpha)/d2(0.5), ...
%      alpha, d3(alpha)/d3(0.5),alpha, d4(alpha)/d4(0.5),...
%      alpha, f2(alpha)/f2(0.5), ...
%      alpha, f3(alpha)/f3(0.5),alpha, f4(alpha)/f4(0.5),...
%      alpha, P, '^-')

Pdata = zeros(size(alpha,1), num_funcs(2)+1);
figure()
hold on
for i = 1: num_funcs(2)
    plot(alpha, funcs{i}(alpha)/funcs{i}(0.5), styles{i})
    Pdata(:,i) = funcs{i}(alpha)/funcs{i}(0.5);
end
plot(alpha, P, 'kp-', 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',10)
Pdata(:,end) = P;

legend(labels)
xlim([0.1 0.8])
ylim([0 2])
hold off


% plot_data = [p2(alpha)/p2(0.5) p3(alpha)/p3(0.5) p4(alpha)/p4(0.5),...
%              a2(alpha)/a2(0.5) a3(alpha)/a3(0.5), ...
%              a4(alpha)/a4(0.5) ...
%              d1(alpha)/d1(0.5) d2(alpha)/d2(0.5), ...
%              d3(alpha)/d3(0.5) d4(alpha)/d4(0.5)...
%              r1(alpha)/r1(0.5) r2(alpha)/r2(0.5) r3(alpha)/r3(0.5)...
%              f2(alpha)/f2(0.5)  f3(alpha)/f3(0.5) f4(alpha)/f4(0.5)];
% 
% sigma2 = zeros(size(plot_data,2),1);
% for i = 1: size(plot_data,2)
%     sigma2(i) = sum((P(1:end)-plot_data(1:end,i)).^2);
% end

% Model = ['P2 '; 'P3 '; 'P4 '; 'A2 '; 'A3 '; 'A4 '; ...
%          'D1 ';'D2 ';'D3 ';'D4 '; 'R1 '; 'R2 '; 'R3 ';'F2 '; 'F3 '; 'F4 '];
% disp('the closest model is:')
% disp(Model(sigma2==min(sigma2),:));
% disp('sigma2 == ')
% disp(min(sigma2))

% R2 = 1 - sigma2./sum((P-mean(P)).^2);

end