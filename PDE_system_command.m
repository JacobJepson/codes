clear;
T =200;
mu = 1;
eta = 1;
kappa = 0.7;
s = 0.2;
omega=0.03;
L_0 = 1;

T_span =[0;T];  
N = 200; %number of space steps 
h = 1/N; %space step length 
x = linspace(0,1,N+1).'; 
init_cond = [omega*(1-x.^2);L_0];


[T_out,Y_out] = ode23s(@(t,y) PDE_system_run(t,y,x,h,kappa,s,mu,eta),T_span,init_cond);

all_out = Y_out(:,1:end-1).'; %Matrix of cell volume fraction
L=Y_out(:,end); %Vector of moving boundary
figure(1)
plot(T_out, L, 'k', 'LineWidth', 1);
xlabel('$t$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$L$', 'Interpreter', 'Latex', 'FontSize', 15);
title('Moving boundary', 'Interpreter', 'Latex', 'FontSize', 15);

%Rescales from moving boundary parameterisation
for j=1:size(L)
    x_new(:,j) = x(:)*L(j);
end

figure(2) %% CELL PLOT %%
for i = 1:size(L)
plot(x_new(:,i),all_out(:,i),'Color', 'k', 'LineWidth', 1);
hold on
end
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$n$', 'Interpreter', 'Latex', 'FontSize', 15);
title('Cell volume fraction', 'Interpreter', 'Latex', 'FontSize', 15);
ylim([0 1])
