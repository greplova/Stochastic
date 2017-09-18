function out_me = forward_master_BAYES(c,cl,v1,v2,v3,params,rho,TN,dNs,dt,Tend,Iexp)
%Let us first pre-allocate the solution:
U = zeros(9,9);

[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [1;0;0;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_1 = Rho_pre(:,end);

[T_pre,Rho_pre]=ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [1;0;0;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_1H = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;1;0;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_2 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [0;1;0;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_2H = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;1;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_3 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [0;0;1;0;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_3H = Rho_pre(:,end);


[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;1;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_4 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;1;0;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_4H = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;1;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_5 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;1;0;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_5H = Rho_pre(:,end);


[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;1;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_6 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;1;0;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_6H = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;1;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_7 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;1;0;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_7H = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;0;1;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_8 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;0;1;0]);
Rho_pre = transpose(Rho_pre);
rho_pre_8H = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereq(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;0;0;1]);
Rho_pre = transpose(Rho_pre);
rho_pre_9 = Rho_pre(:,end);

[T_pre,Rho_pre] = ode45(@(t,x) mastereqH(t,x,c,cl,params), 0:dt/10:dt, [0;0;0;0;0;0;0;0;1]);
Rho_pre = transpose(Rho_pre);
rho_pre_9H = Rho_pre(:,end);

U = [rho_pre_1 rho_pre_2 rho_pre_3 rho_pre_4 rho_pre_5 rho_pre_6 rho_pre_7 rho_pre_8 rho_pre_9];
UH = [rho_pre_1H rho_pre_2H rho_pre_3H rho_pre_4H rho_pre_5H rho_pre_6H rho_pre_7H rho_pre_8H rho_pre_9H];

%Now we start
rhov=reshape(rho,9,1); %saving rho as a vector
t0=0;

%We alculate the number of timesteps
L_end = round(Tend/dt);
%How often we should save?
L_save = 1;
T_me = zeros(1,(L_end)/L_save+1);
Rho_me = zeros(9,(L_end)/L_save+1);
EdNUpDown = zeros(1,(L_end)/L_save+1);
EdN0 = zeros(1,(L_end)/L_save+1);
size(dNs)
size(T_me)

T_me(1) = 0;
Rho_me(:,1) = rhov;


for L=1:L_end
    Rhot = U*rhov;
    rhome = reshape(Rhot,3,3);
    rhonorm = trace(rhome);
    
    rhome = rhome/rhonorm;
    Rhot = Rhot/rhonorm;
    
    t0 = t0+dt;
    if mod(L,1)==0
        L_var = (L)/L_save+1;
        T_me(L_var) = t0;
        Rho_me(:,L_var) = Rhot;
        if dNs(L) == 1
            EdN0(L_var) = params.eta*params.tau*params.tau*dt;
            EdNUpDown(L_var) = params.eta*(params.xi*params.xi+2*params.xi*params.tau+params.tau*params.tau)*dt;
        else 
            EdN0(L_var) = 1-params.eta*params.tau*params.tau*dt;   
            EdNUpDown(L_var) = 1-params.eta*(params.xi*params.xi+2*params.xi*params.tau+params.tau*params.tau)*dt;
        end
        %EdN(L_var)=params.tau*params.tau+((params.xi*params.xi-2*params.xi*params.tau)*dNs((L-99):L))/100;
    end
    
    
    r = sqrt(params.eta)*(params.tau*eye(3)+params.xi*(v3*v3' + v2*v2'));
    n = v3*v3' + v2*v2';
    HW = n*rhome + rhome*n-trace(n*rhome + rhome*n)*rhome; %Wiseman H operator
    HW = -sqrt(params.eta)*params.xi*HW;
    dW = -(dt*(Iexp(L)/(sqrt(params.eta)*params.tau) - (sqrt(params.eta)*params.tau + 2*sqrt(params.eta)*params.xi*trace(rhome*n))));
    dW = min(max(dW,-10*sqrt(dt)),10*sqrt(dt));
    if min(diag(rhome))<-0.1
        disp('Your time-step is most likely too large')
        pause()
    end
    cW=-sqrt(params.eta)*params.xi*n; %WiMi2010 4.322
    dY(L)= dW/dt+trace((cW+cW')*rhome); %PQSsuppA.38
    drho = dW*HW;
    
    
    %drho = dNs(L)*G;
    rhome = drho + rhome;
    rhome = rhome/trace(rhome);

    %trace(rho)
    rhov = reshape(rhome,9,1);
    rhov_norm = sqrt(rhov'*rhov);
end

out_me{1} = T_me;
out_me{2} = Rho_me;
out_me{3} = EdNUpDown;
out_me{4} = EdN0;
out_me{11} = dY;

%rho|x=1,dN=0
out_me{5} = rho_pre_1/(rho_pre_1(1)+rho_pre_1(5)+rho_pre_1(9));
%rho|x=1,dN=1
rho_tmp = rho_pre_1/(rho_pre_1(1)+rho_pre_1(5)+rho_pre_1(9));
rho_tmp = reshape(rho_tmp,3,3);
rho_tmp = r*rho_tmp*r'/(trace(r*rho_tmp*r'));
rho_tmp = reshape(rho_tmp,9,1);
out_me{6} = rho_tmp;

% rho|x=5,dN=0
out_me{7} = rho_pre_5/(rho_pre_5(1)+rho_pre_5(5)+rho_pre_5(9));
% rho|x=5,dN=1
rho_tmp = rho_pre_5/(rho_pre_5(1)+rho_pre_5(5)+rho_pre_5(9));
rho_tmp = reshape(rho_tmp,3,3);
rho_tmp = r*rho_tmp*r'/(trace(r*rho_tmp*r'));
rho_tmp = reshape(rho_tmp,9,1);
out_me{8} = rho_tmp;

% rho|x=9,dN=0
out_me{9} = rho_pre_9/(rho_pre_9(1)+rho_pre_9(5)+rho_pre_9(9));
% rho|x9,dN=1
rho_tmp = rho_pre_9/(rho_pre_9(1)+rho_pre_9(5)+rho_pre_9(9));
rho_tmp = reshape(rho_tmp,3,3);
rho_tmp = r*rho_tmp*r'/(trace(r*rho_tmp*r'));
rho_tmp = reshape(rho_tmp,9,1);
out_me{10} = rho_tmp;

out_me{12} = rho_pre_1H/(rho_pre_1H(1)+rho_pre_1H(5)+rho_pre_1H(9));
out_me{13} = rho_pre_5H/(rho_pre_5H(1)+rho_pre_5H(5)+rho_pre_5H(9));
out_me{14} = rho_pre_9H/(rho_pre_9H(1)+rho_pre_9H(5)+rho_pre_9H(9));

function drhov = mastereq(t,rhov, c, cl, params)
rho = reshape(rhov,3,3);
id3 = speye(3);%identity
v3 = id3(:,1);%spin-up
v2 = id3(:,2);%spin-down
v1 = id3(:,3);%empty dot
H = 0.5*params.GuessOmega*(v2*v3'+ v3*v2'); %Hamiltonian
%r = (2*params.tau*params.xi*(v3*v3' + v2*v2')+(params.xi*(v3*v3' + v2*v2'))^2);
HH = -1i*(H*rho - rho*H); %commutator with Hamiltonian %- params.eta*0.5*(r*rho + rho*r');% - trace(r*rho + rho*r')*rho; %Wiseman's H operator
cup = (v1*v3');
cdown = (v2*v1');
Dup = cup*rho*cup'-0.5*(cup'*cup)*rho-0.5*rho*(cup'*cup); %dissipator - annihilation operator
Ddown = cdown*rho*cdown'-0.5*(cdown'*cdown)*rho-0.5*rho*(cdown'*cdown); % dissipator - creation operator
%q=params.tau*eye(3)+params.xi*(v3*v3' + v2*v2');
n=v3*v3' + v2*v2';
Dnon = n*rho*n'-0.5*(n'*n)*rho-0.5*rho*(n'*n);
Dnon = (params.xi*params.xi)*Dnon;
drho = HH + params.Guessgammadown*Ddown + params.Guessgammaup*Dup + Dnon;

drhov = reshape(drho,9,1);

function drhov = mastereqH(t,rhov, c, cl, params)
rho = reshape(rhov,3,3);
id3 = speye(3);%identity
v3 = id3(:,1);%spin-up
v2 = id3(:,2);%spin-down
v1 = id3(:,3);%empty dot
H = 0.5*params.GuessOmega*(v2*v3'+ v3*v2'); %Hamiltonian
%r = (2*params.tau*params.xi*(v3*v3' + v2*v2')+(params.xi*(v3*v3' + v2*v2'))^2);
HH = -1i*(H*rho - rho*H); %commutator with Hamiltonian 
drho = HH;

drhov = reshape(drho,9,1);
