
k = 0.5; %kappa
init = 1; %choose sensible initial value
E = @(x) phiPP(x,k);
wavespeed = fzero(E,init);


function E=phiPP(U,k)
s=0.2; mu = 1; eta=1;

S=@(x,y) (y*(U+y)+x*(k-x)*(   (  (1-s)*(1-(k-x)-s)/(  (eta*(1-s)+s)*(1+(mu-1)*(k-x)-s))   )*( (1-s)/((1-(k-x)-s)^2)  +eta-1)      ))/((k-x)*y);
n_axis = linspace(1e-15,k);
tol=-1e-15; 
sol = ode23s(@(t,x) S(t,x), n_axis,tol);
X=k-sol.x; %re-scale 

figure(3)
hold on
plot(X,sol.y,'Color',[0.7 0.7 0.7])
axis([0 k -1 0])
box on
xlabel('$n$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$q$', 'Interpreter', 'Latex', 'FontSize', 15);
title('heteroclinic connection $q(n)$', 'Interpreter', 'Latex', 'FontSize', 15);
E=U+sol.y(end);
end
