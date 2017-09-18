function out_E = backward_master_BAYES(c,cl,v1,v2,v3,params,Ev,TN_f,dNs_f,dt,Tend,dY)
%Let us first pre-allocate the solution:
U = zeros(9,9);
[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [1;0;0;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_1 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;1;0;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_2 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;1;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_3 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;1;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_4 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;1;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_5 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;1;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_6 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;1;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_7 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;0;1;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_8 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) ad_mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;0;0;1]);
Rho_pre = transpose(Rho_pre);
rho_pre_9 = Rho_pre(:,end);

U = [rho_pre_1 rho_pre_2 rho_pre_3 rho_pre_4 rho_pre_5 rho_pre_6 rho_pre_7 rho_pre_8 rho_pre_9];



%We calculate the number of timesteps
L_end = round(Tend/dt);
%How often we should save?
L_save = 1;

TE_me = zeros(1,(L_end)/L_save+1);
E_me = zeros(9,(L_end)/L_save+1);
%Now we start
TE_me(1) = 0;
E_me(:,1) = Ev;

t0p=0;

dYb=flip(dY);
%dYb=[dYb 0];

for L=1:L_end
    Et = U*Ev;
    E = reshape(Et,3,3);
    Enorm = trace(E);
    E = E/Enorm;
    Et = Et/Enorm;
    
    t0p = t0p +dt;
    if mod(L,1)==0
        L_var = (L)/L_save+1;
        TE_me(L_var) = t0p;
        E_me(:,L_var) = Et;
    end
    

    %r = sqrt(params.eta)*(params.tau*eye(3)+params.xi*(v3*v3' + v2*v2'));
    %G = r'*E*r-E;
    n = v3*v3' + v2*v2';
    HW = n*E + E*n; %-trace(n*E + E*n)*E; %Wiseman H operator
    HW = -sqrt(params.eta)*params.xi*HW;
    
    %dYYb=dYb(L+1);
    dE = dt*dYb(L)*HW;
    E = dE + E;
    E = E/trace(E);
   
    Ev = reshape(E,9,1);
    %rhov_norm = sqrt(rhov'*rhov);
end

%E_me

TE_me=Tend-TE_me;
TE_me = fliplr(TE_me);
E_me = fliplr(E_me);

out_E{1} = TE_me;
out_E{2} = E_me;



function dEv = ad_mastereq(t,Ev,c,cl,params)
E = reshape(Ev,3,3);
id3 = speye(3);%identity
v3 = id3(:,1);%spin-up
v2 = id3(:,2);%spin-down
v1 = id3(:,3);%empty dot
H = 0.5*params.GuessOmega*(v2*v3'+ v3*v2'); %Hamiltonian
%r = (2*params.tau*params.xi*(v3*v3' + v2*v2')+(params.xi*(v3*v3' + v2*v2'))^2);
HH = 1i*(H*E - E*H); %Hamiltonian commutator %- params.eta*0.5*(r*E + E*r'); %trace(r*E + E*r')*E; %Wiseman's H operator
cup = (v1*v3');
cdown = (v2*v1');

Dup = cup'*E*cup-0.5*(cup'*cup)*E-0.5*E*(cup'*cup); %dissipator - annihilation operator
Ddown = cdown'*E*cdown-0.5*(cdown'*cdown)*E-0.5*E*(cdown'*cdown); % dissipator - creation operator
%q=params.tau*eye(3)+params.xi*(v3*v3' + v2*v2');
n=v3*v3' + v2*v2';
Dnon = n'*E*n-0.5*(n'*n)*E-0.5*E*(n'*n);
Dnon = (params.xi*params.xi)*Dnon;
%Dnon = (1-params.eta)*Dnon;

dE = HH + params.Guessgammadown*Ddown + params.Guessgammaup*Dup + Dnon;

dEv = reshape(dE,9,1);
