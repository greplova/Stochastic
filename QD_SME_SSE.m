function out = QD_SME_SSE(Tin,Rin,params)
%parpool: (used for parallel cluster simulation)
%parpool(20)

global R %Random number used for the algorithm


%Here we define parameters, if they are not specificied by the function
%call.
if exist('params','var') == 0
    params.Omega = 5; 
    params.gammaup = 3;
    params.gammadown = 3;
    params.tau=sqrt(62400);
    params.xi=-sqrt(62400)/100;
    params.eta=0.5;
    params.GuessOmega=5;
    params.Guessgammaup=3;
    params.Guessgammadown=3;
end

%Setting Random-seed
%If Rseed is not set, Grendel will draw the same random number in every run
if exist('Rin','var') == 0
    Rseed = 12;
else
    Rseed = str2num(Rin);
end
s=RandStream('mt19937ar','Seed',Rseed);
RandStream.setGlobalStream(s);

%Setting time parameters
if exist('Tin','var') == 0
    Tend = 10;
else
    Tend = str2num(Tin);
end

%time step length
%dt = 0.0002;
dt = 0.00002;

%Set initial state
psi0 = [0;1;0]; %start in 'spin-down' state

%The c's are our jump operators. Just add more by c{n}
id3 = speye(3);%identity
v3 = id3(:,1);%spin-up
v2 = id3(:,2);%spin-down
v1 = id3(:,3);%empty dot

c{1} = sqrt(params.gammaup)*(v1*v3');%annihilation operator for spin up
c{2} = sqrt(params.gammadown)*(v2*v1');%creation operator for spin down
c{3} = (params.tau+params.xi)*eye(3);
cl = length(c);

c1 = sqrt(params.gammaup)*(v1*v3'); % annihilation operator for spin-up - matrix form
c2 = sqrt(params.gammadown)*(v2*v1'); %creation operator for spin-down - matrix form
c3 = sqrt(params.eta)*(params.tau*eye(3)+params.xi*(v3*v3' + v2*v2'));

R = rand(1);

%%%%%%%%%%%%%%%%%
%% forward SSE %%
%%%%%%%%%%%%%%%%%
tic %<--  start timer

%Now we call the function to run forward_sse
%This will make sure that we only save,
%the variable that we need, so we save memory
out_sse = forward_sse(c,cl,c1,c2,c3,params,dt,Tend,psi0);

%Here we save the output
Ts = out_sse{1};
Psis = out_sse{2};
TN = out_sse{3};
dNs = out_sse{4};
Psi_norm = out_sse{5};
Rho_psis = out_sse{6};
T_jumps = out_sse{7};
Iexp = out_sse{8};


toc %<--   stop timer
%pause(0.1)

figure(1),plot(Ts,abs(Psis(1,:)).^2./Psi_norm) %plot of quantum trajectory

figure(23), plot(min(reshape(TN,[],1000)), mean(reshape(1.6*0.0001*Iexp,[],1000)))

%hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% backward evolution for SSE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare T_jumps for the backward evolution
T_jumps=flip(T_jumps);
T_jumps=Tend-T_jumps;

disp('backward')
tic
psib0=Psis(:,end); %last Psi


%Now we call the function to run forward_sse
%This will make sure that we only save,
%the variable that we need, so we save memory

out_sse = backward_sse(c,cl,c1,c2,c3,params,dt,Tend,psib0,T_jumps);

%Here we save the output

Tbs = out_sse{1};
Psibs = out_sse{2};
Psib_norm = out_sse{3};
Rho_psibs = out_sse{4};
toc

%figure(7), plot(Tbs, abs(Psibs(1,:)).^2./Psib_norm) %plot of the backward quantum trajecotry
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparation for the Master Equations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Set initial state for ME
psi00 = [0;1;0]; %spin-down state inside the dot 
rho=psi00*psi00';

%Prepare elements for backward evolution
%flip dNs
dNs_f=flip(dNs);
dNs_f=[dNs_f 0];
TN_f=flip(TN);
TN_f=TN(end)-TN_f;
TN_f=[TN_f TN(end)];

E=eye(3)/3;
Ev=reshape(E,9,1); %Rho_me(:,end); % the last rho

%This is for the parrallel run for many omegas
%Omegas = linspace((params.Omega*0.5),(params.Omega*1.5),17);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forward Master Equation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Forward master')
tic
out_me = forward_master_BAYES(c,cl,v1,v2,v3,params,rho,TN,dNs,dt,Tend,Iexp);
T_me = out_me{1};
Rho_me = out_me{2};
dY = out_me{11};

toc



 figure(3),plot(T_me,abs(Rho_me(1,:)+abs(Rho_me(5,:))))
 
 Rho_me(:,end)
 pause(0.1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adjoint evolution for SME %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Backward master')

%set initial state
E=eye(3)/3;
Ev=reshape(E,9,1); %Rho_me(:,end); % the last rho

tic
out_E = backward_master_BAYES(c,cl,v1,v2,v3,params,Ev,TN_f,dNs_f,dt,Tend,dY);
TE_me = out_E{1};
E_me = out_E{2};
toc

max(Iexp)
mean(Iexp)
var(Iexp)

tic
out_E_tmp = backward_master_BAYES_TMP(c,cl,v1,v2,v3,params,Ev,TN_f,dNs_f,dt,Tend,dY_tmp);
TE_me_tmp = out_E_tmp{1};
E_me_tmp = out_E_tmp{2};
toc

% figure(8),plot(T_me,abs(E_me(1,:)+abs(E_me(5,:))))
% pause()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluating the probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RHO ZERO TEMPERATURE
T_final=linspace(0,Tend,30000); %preparing the grid
pE1=E_me(1,:); %picking elements of E we need for evaluation of probability <up|E|up>
pE5=E_me(5,:); %<down|E|down>
pE2=E_me(2,:); %<up|E|down>
pE4=E_me(4,:); %<down|E|up>
pE9=E_me(9,:); %<0|E|0>

pE1_final=interp1_failsafe(TE_me,pE1,T_final);%E_me on the grid
pE5_final=interp1_failsafe(TE_me,pE5,T_final);
pE2_final=interp1_failsafe(TE_me,pE2,T_final);
pE4_final=interp1_failsafe(TE_me,pE4,T_final);
pE9_final=interp1_failsafe(TE_me,pE9,T_final);

pRho1=Rho_me(1,:); %<up|Rho|up>
pRho5=Rho_me(5,:); %<down|Rho|down>
pRho2=Rho_me(2,:); %<up|Rho|down>
pRho4=Rho_me(4,:); %<down|Rho|up>
pRho9=Rho_me(9,:); %<0|Rho|0>

pRho1_final=interp1_failsafe(T_me,pRho1,T_final);
pRho2_final=interp1_failsafe(T_me,pRho2,T_final);
pRho5_final=interp1_failsafe(T_me,pRho5,T_final);
pRho4_final=interp1_failsafe(T_me,pRho4,T_final);
pRho9_final=interp1_failsafe(T_me,pRho9,T_final);

Pseparate_final=(pE1_final.*pRho1_final+pE5_final.*pRho5_final)./(pE1_final.*pRho1_final+pE5_final.*pRho5_final+pE9_final.*pRho9_final);

figure(9), plot(T_final,Pseparate_final);


Pjoint2_final=(pE1_final.*pRho1_final+pE5_final.*pRho5_final+pE2_final.*pRho4_final+pE4_final.*pRho2_final)./(pE1_final.*pRho1_final+pE5_final.*pRho5_final+pE2_final.*pRho4_final+pE4_final.*pRho2_final+pE9_final.*pRho9_final);

figure (11), plot(T_final,Pjoint2_final)


%PSI
pRf1=Rho_psis(1,:);%picking elements of density matrix for forward SSE we need for evaluation of probability <up|E|up>
pRf5=Rho_psis(5,:); %<down|E|down>
pRf2=Rho_psis(2,:); %<up|E|down>
pRf4=Rho_psis(4,:); %<down|E|up>
pRf9=Rho_psis(9,:); %<0|E|0>
pRf1_final=interp1_failsafe(Ts,pRf1,T_final);%Rho_psis on the grid
pRf5_final=interp1_failsafe(Ts,pRf5,T_final);
pRf2_final=interp1_failsafe(Ts,pRf2,T_final);
pRf4_final=interp1_failsafe(Ts,pRf4,T_final);
pRf9_final=interp1_failsafe(Ts,pRf9,T_final);
pRb1=Rho_psibs(1,:);%picking elements of density matrix for forward SSE we need for evaluation of probability <up|E|up>
pRb5=Rho_psibs(5,:); %<down|E|down>
pRb2=Rho_psibs(2,:); %<up|E|down>
pRb4=Rho_psibs(4,:); %<down|E|up>
pRb9=Rho_psibs(9,:); %<0|E|0>

pRb1_final=interp1_failsafe(Tbs,pRb1,T_final);%Rho_psibs on the grid
pRb5_final=interp1_failsafe(Tbs,pRb5,T_final);
pRb2_final=interp1_failsafe(Tbs,pRb2,T_final);
pRb4_final=interp1_failsafe(Tbs,pRb4,T_final);
pRb9_final=interp1_failsafe(Tbs,pRb9,T_final);

PPSIseparate_final=(pRb1_final.*pRf1_final+pRb5_final.*pRf5_final)./(pRb1_final.*pRf1_final+pRb5_final.*pRf5_final+pRb9_final.*pRf9_final);

figure(12), plot(T_final,PPSIseparate_final);

save('simulation_oxf_100.mat');



function out = trace_v_mat(Rhov,Ev)
rho1 = Rhov{1};
rho2 = Rhov{2};
rho3 = Rhov{3};
rho4 = Rhov{4};
rho5 = Rhov{5};
rho6 = Rhov{6};
rho7 = Rhov{7};
rho8 = Rhov{8};
rho9 = Rhov{9};

E1 = Ev{1};
E2 = Ev{2};
E3 = Ev{3};
E4 = Ev{4};
E5 = Ev{5};
E6 = Ev{6};
E7 = Ev{7};
E8 = Ev{8};
E9 = Ev{9};

%out = E1.*rho1 + E2.*rho4 + E4.*rho2 + E3.*rho7 + E5.*rho5 + E7.*rho3 + E6.*rho8 + E8.*rho6 + E9.*rho9;
out = E1.*rho1 + E5.*rho5 + E7.*rho3 + E9.*rho9;


function out = trace_v_mat_2(Rhov,Ev)
rho1 = Rhov(1);
rho2 = Rhov(2);
rho3 = Rhov(3);
rho4 = Rhov(4);
rho5 = Rhov(5);
rho6 = Rhov(6);
rho7 = Rhov(7);
rho8 = Rhov(8);
rho9 = Rhov(9);

E1 = Ev{1};
E2 = Ev{2};
E3 = Ev{3};
E4 = Ev{4};
E5 = Ev{5};
E6 = Ev{6};
E7 = Ev{7};
E8 = Ev{8};
E9 = Ev{9};

%out = E1.*rho1 + E2.*rho4 + E4.*rho2 + E3.*rho7 + E5.*rho5 + E7.*rho3 + E6.*rho8 + E8.*rho6 + E9.*rho9;
out = E1.*rho1 + E5.*rho5 + E9.*rho9;




function yy = interp1_failsafe(x,y,xx)
x_tmp = [x(diff(x)~=0) x(end)];
y_tmp = [y(diff(x)~=0) y(end)];
yy = interp1(x_tmp,y_tmp,xx);


function out=save_output(N_omega, Rseed, T_me, Rho_me, TN, dNs, dNs_f, TN_f, TE_me, E_me, Tend, T_final, Pseparate_final, Pjoint2_final, PPSIseparate_final, PPSIjoint2_final)

save(['/home/eliska/QuantumDot/output_PQS' num2str(N_omega) '_T' num2str(Tend) '_R' num2str(Rseed) '.mat'],'T_me', 'Rho_me', 'TN', 'dNs', 'dNs_f', 'TN_f', 'TE_me', 'E_me', 'Tend', 'T_final', 'Pseparate_final', 'Pjoint2_final', 'PPSIseparate_final', 'PPSIjoint2_final') 

out = 1;





function drhopv = mastereq_prob(t,rhopv, c, cl, params) %function for evaluation of survival probability
rhop = reshape(rhopv,3,3);
id3 = speye(3);%identity
v3 = id3(:,1);%spin-up
v2 = id3(:,2);%spin-down
v1 = id3(:,3);%empty dot
H = 0.5*params.Omega*(v2*v3'+ v3*v2'); %Hamiltonian
r = (1i*H+2*params.tau*params.xi*(v3*v3' + v2*v2')+(params.xi*(v3*v3' + v2*v2'))^2);
HH = r*rhop + rhop*r' - trace(r*rhop + rhop*r')*rhop; %Wiseman's H operator
cup = (v1*v3');
cdown = (v2*v1');
Dup = -cup'*cup*rhop-rhop*cup'*cup; %dissipator - annihilation operator
Ddown = -cdown'*cdown*rhop-rhop*cdown'*cdown; % dissipator - creation operator
q=params.tau*eye(3)+params.xi*(v3*v3' + v2*v2');
Dnon = -q'*q*rhop-rhop*q'*q;
drhop = -0.5*HH + params.gammadown*Ddown + params.gammaup*Dup + Dnon;

drhopv = reshape(drhop,9,1);


function [value,isterminal,direction] = rand_events(t,y)
global R
psi = y;
value = psi'*psi - R;     % Detect norm = R
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only

