function out_sse = backward_sse(c,cl,c1,c2,c3,params,dt,Tend,psib,T_jumps)
%variables that we need
t0b = 0;

%Option for the ode-solver
options = odeset('Events',@rand_events);

%Let us first pre-allocate the solution:
U = zeros(3,3);
[Tb_pre,Psib_pre] = ode45(@(t,x) backward_schrodinger(t,x,c,cl,params), 0:dt/10:dt, [1;0;0]);
Psib_pre = transpose(Psib_pre);
psib_pre_1 = Psib_pre(:,end);
[Tb_pre,Psib_pre] = ode45(@(t,x) backward_schrodinger(t,x,c,cl,params), 0:dt/10:dt, [0;1;0]);
Psib_pre = transpose(Psib_pre);
psib_pre_2 = Psib_pre(:,end);
[Tb_pre,Psib_pre] = ode45(@(t,x) backward_schrodinger(t,x,c,cl,params), 0:dt/10:dt, [0;0;1]);
Psib_pre = transpose(Psib_pre);
psib_pre_3 = Psib_pre(:,end);
Ub = [psib_pre_1 psib_pre_2 psib_pre_3];


%We alculate the number of timesteps
M_end = round(Tend/dt);
%How often we should save?
M_save = 100;

%now we can pre-allocate our variables:
Tbs = zeros(1,(M_end)/M_save+1);
Psibs = zeros(3,(M_end)/M_save+1);

%We save the inital point
Tbs(1) = 0;
Psibs(:,1) = psib;

%M_jump keeps track of which jump we are at
M_jump = 1;
%and T_jump is the time of that jump
T_jump = T_jumps(M_jump);
for M = 1:M_end
    t0b=t0b+dt;
    psib=Ub*psib;
    psib = psib/norm(psib);
    
    %we check is we have a jump
    if t0b>=T_jump
        M_jump = M_jump+1;
        %check if we had the last jump
        %if so, we make sure that no more jump happens
        if M_jump>length(T_jumps)
            T_jump = Tend+dt;
        else
        %otherwise we save the next tump time in T_jump
            T_jump = T_jumps(M_jump);
        end
        %find out if the spin is in the dot or not and apply respective
        %operator
        norm_dot_t=round(abs(psib(1))^2+abs(psib(2))^2); %is spin in up-state?
        if norm_dot_t == 1
            psib = c2'*psib; %appylying annihilation operator (spin- up leaves the dot)
        else
            psib = c1'*psib; %otherwise there was no spin to start with so we apply creation operator to get the spin into the dot
        end
    end
    
    %norm_dot_t_2=round(abs(psib(1))^2); %check spin-up again
    
    %if norm_dot_t_2 == 1
    %    psi0 = c1*psib; %get flipped spin out of the dot
    %else
    %    psi0 = c2*psib; %otherwise there was no spin to start with so we apply creation operator to get the spin into the dot
    %end
    psib=psib/norm(psib);%normalize the wavefuction
    
    if mod(M,100)==0
        M_var = (M)/M_save+1;
        Tbs(M_var) = t0b;
        Psibs(:,M_var) = psib;
    end

end

Tbs=Tend - Tbs;
Psib_norm = zeros(1,length(Psibs));
for N=1:length(Psib_norm)
    Psib_norm(N) = Psibs(:,N)'*Psibs(:,N);
end

%create density matrix representation for resulting wavefunction
Rho_psibs = zeros(9,length(Tbs));
for n=1:length(Tbs)
   psinb = Psibs(:,n);
   %Rho_psis{n} = psin*psin'; 
   Rho_psibs(:,n) = reshape(psinb*psinb'/Psib_norm(n)^2,9,1); %check that is 9,1 and not 1,9
   %Rho_psibs = [Rho_psibs reshape(psinb*psinb'/Psib_norm(n)^2,9,1)]; %check that is 9,1 and not 1,9
end

out_sse{1} = Tbs;
out_sse{2} = Psibs;
out_sse{3} = Psib_norm;
out_sse{4} = Rho_psibs;

function dpsib = backward_schrodinger(t,x,c,cl,params)
id3 = speye(3);%identity
v3 = id3(:,1);%spin-up
v2 = id3(:,2);%spin-down
H = 0.5*params.Omega*(v2*v3'+ v3*v2'); %Hamiltonian
dpsib = 0;
for N = 1:cl
    dpsib = dpsib - ((c{N})'*c{N})*x; %Nth jump opertor
end
dpsib =  dpsib + 1i*H*x;