function output = PDE_system_run(~,input,x,h,kappa,s,mu,eta)
S=input(end); 
n=[input(1:end-2);0]; %n=0 at  the boundary

np = [n(2:end);NaN]; %n_{i+1}
nm = [NaN;n(1:end-1)]; %n_{i-1}

dsdt = -(1/(2*S*h))*(-4*n(end-1)+n(end-2)); %moving boundary


advec = [0;x(2:end-1).*dsdt.*(np(2:end-1)-n(2:end-1))/(S*h);NaN]; %advective term 

nphi = n.*(( (1-s).*(1-n-s))./(  (eta.*(1-s)+s).*(1+(mu-1).*n-s)  )).*( (1-s)./((1-n-s).^2)  +eta-1) ; 

%second derivative: nphi(n)*n_xx
d2ndx2 = [ 2*nphi(1)*(n(2)-n(1))/(h^2*S^2)  ;   nphi(2:end-1).*(np(2:end-1)-2*n(2:end-1)+nm(2:end-1))/(S^2*h^2) ;NaN];
%nonlinear derivative term: (n_x)^2 * (n*phi(n))_n
nonlindx = [0;(1/(4*h^2*S^2)).*(np(2:end-1)-nm(2:end-1)).^2.*(  (2.*n(2:end-1).*(s - 1).^2)./((s - eta.*(s - 1)).*(s + n(2:end-1) - 1).^2.*(n(2:end-1).*(mu - 1) - s + 1)) - ((s - 1).*((s - 1)./(s + n(2:end-1) - 1).^2 - eta + 1).*(s + n(2:end-1) - 1))./((s - eta*(s - 1)).*(n(2:end-1).*(mu - 1) - s + 1)) - (n(2:end-1).*(s - 1).*((s - 1)./(s + n(2:end-1) - 1).^2 - eta + 1))./((s - eta*(s - 1)).*(n(2:end-1).*(mu - 1) - s + 1)) + (n(2:end-1).*(mu - 1).*(s - 1).*((s - 1)./(s + n(2:end-1) - 1).^2 - eta + 1).*(s + n(2:end-1) - 1))./((s - eta.*(s - 1)).*(n(2:end-1).*(mu - 1) - s + 1).^2));NaN];
  

dndt = [advec(1:end-1)+d2ndx2(1:end-1)+nonlindx(1:end-1)+n(1:end-1).*(kappa-n(1:end-1));0];
output=[dndt;dsdt];


end 
