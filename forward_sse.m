function out_sse = forward_sse(c,cl,c1,c2,c3,params,dt,Tend,psi)
%variables that we need
t0 = 0;
R = rand(1);


%Option for the ode-solver
options = odeset('Events',@rand_events);

%Let us first pre-allocate the solution:
U = zeros(3,3);
[T_pre,Psi_pre] = ode45(@(t,x) schrodinger(t,x,c,cl,params), 0:dt/10:dt, [1;0;0]);
Psi_pre = transpose(Psi_pre);
psi_pre_1 = Psi_pre(:,end);
[T_pre,Psi_pre] = ode45(@(t,x) schrodinger(t,x,c,cl,params), 0:dt/10:dt, [0;1;0]);
Psi_pre = transpose(Psi_pre);
psi_pre_2 = Psi_pre(:,end);
[T_pre,Psi_pre] = ode45(@(t,x) schrodinger(t,x,c,cl,params), 0:dt/10:dt, [0;0;1]);
Psi_pre = transpose(Psi_pre);
psi_pre_3 = Psi_pre(:,end);
U = [psi_pre_1 psi_pre_2 psi_pre_3];


%We alculate the number of timesteps
LS_end = round(Tend/dt);
%How often we should save?
L_save = 100;

%now we can pre-allocate our variables:
TN = zeros(1,LS_end);
dNs = zeros(1,LS_end);
Ts = zeros(1,(LS_end)/L_save+1);
Psis = zeros(3,(LS_end)/L_save+1);

%We save the inital point
Ts(1) = 0;
Psis(:,1) = psi;

%We pre-allocate the random numbers used to generate dNs
RNs = rand(1,LS_end);
Iexp = zeros(1,LS_end);

%Keep track of jumps
M_jump = 1;


for LS = 1:LS_end
    psi = U*psi;
    t0 = t0+dt;
    
    %generation of experimental current to creat dW's in SME
    id3 = speye(3);%identity
    v3 = id3(:,1);%spin-up
    v2 = id3(:,2);%spin-down
    %v1 = id3(:,3);%empty dot
    %Iexp(LS)= sqrt(params.eta)*params.tau*(sqrt(params.eta)*params.tau - sqrt(params.eta)*params.xi*psi'*(v3*v3' + v2*v2')*psi/(psi'*psi) - randn(1,1)*sqrt(dt)/dt);
    
    %Here we generate the dNs
    %Notice that we dont draw any random number 
    %or change sign of any variable
    %We also save at 
    RN = RNs(LS);
    if RN<dt*(psi'*(c3'*c3)*psi)/(psi'*psi)   
        TN(LS) = t0;
        dNs(LS) = 1;
    else
        TN(LS) = t0;
        dNs(LS) = 0;
    end
    
    %We check if we jump
    psi_norm = psi'*psi;
    if psi_norm < R;
        T_jumps(M_jump) = t0;
        M_jump = M_jump + 1;
        R_jump = rand(1);
        P_jump = zeros(1,cl);
        for N = 1:cl
            P_jump(N) = psi'*((c{N}'*c{N})*psi); %Probability for Nth jump
        end
        P_jump = P_jump/sum(P_jump); %Normalize probabilities
        P_jump = cumsum(P_jump); %Make the cumulated probability

        M = find(P_jump>R_jump,1);
        % Finds the 1st index, M, of P_jump that is higher than R_jump

        %Do the quantum jump
        psi = c{M}*psi;
        psi = psi/norm(psi);
        R = rand(1);
    end
     
     Iexp(LS)= sqrt(params.eta)*params.tau*(sqrt(params.eta)*params.tau + 2*sqrt(params.eta)*params.xi*psi'*(v3*v3' + v2*v2')*psi/(psi'*psi) - min(max(randn(1,1),-10),10)*sqrt(dt)/dt);

     %Save the times and Psis. We don't save the last as it will also appear
    %in the next loop.
    if mod(LS,L_save)==0
        L_var = (LS)/L_save+1;
        Ts(L_var) = t0;
        Psis(:,L_var) = psi;
    end
end

%create density matrix representation for resulting wavefunction


Psi_norm = zeros(1,length(Psis));
for N=1:length(Psi_norm)
    Psi_norm(N) = Psis(:,N)'*Psis(:,N);
end

Rho_psis = zeros(9,length(Ts));
for n=1:length(Ts)
   psin = Psis(:,n);
   %Rho_psis{n} = psin*psin'; 
   Rho_psis(:,n) = reshape(psin*psin'/Psi_norm(n)^2,9,1); %check that is 9,1 and not 1,9
end

out_sse{1} = Ts;
out_sse{2} = Psis;
out_sse{3} = TN;
out_sse{4} = dNs;
out_sse{5} = Psi_norm;
out_sse{6} = Rho_psis;
out_sse{7} = T_jumps;
out_sse{8} = Iexp;

function dpsi = schrodinger(t,x,c,cl,params)
id3 = speye(3);%identity
v3 = id3(:,1);%spin-up
v2 = id3(:,2);%spin-down
H = 0.5*params.Omega*(v2*v3'+ v3*v2'); %Hamiltonian

dpsi = 0;
for N = 1:cl
    dpsi = dpsi - ((c{N})'*c{N})/2*x; %Nth jump opertor
end
dpsi =  dpsi - 1i*H*x;
