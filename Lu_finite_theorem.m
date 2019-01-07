function [m,AF,BF,CF,delta,A,B,C,L,Lam1,R]=Lu_finite_theorem()
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper

clear all
clc
% Input:

% Output:
%known parameters
k1=0.1;
k2=0.1;
delta=0.2;
A=cell(1,2);
B=cell(1,2);
C=cell(1,2);
L=cell(1,2);
% 
% A{1}=[-1 0 0;0 -2 1;0 0 -1];
% B{1}=[1 0 0.3]';
% C{1}=[0 0.1 0.3];
% L{1}=[0.5 -0.1 0.1];
% 
% A{2}=[-2 0.3 0;0 -2 0;0 0 -1.8];
% B{2}=[-0.6 0.5 0]';
% C{2}=[0.2 0.2 0];
% L{2}=[0 -0.3 0.2];

if 1
    A{1}=[-9.1 50;-1 -10];
    B{1}=[0 1]';
    C{1}=[1 0];
    L{1}=[1 0];
    
    A{2}=[-0.1 50;-1 -16];
    B{2}=[0 1]';
    C{2}=[1 0];
    L{2}=[1 0];
else
    
    A{1}=[-1 2;-3 0.1];
    B{1}=[-0.2 1]';
    C{1}=[1 0.5];
    L{1}=[0.6 1];
    
    A{2}=[-1 2;-3 0.1];
    B{2}=[-0.2 1]';
    C{2}=[1 0.5];
    L{2}=[0.6 1];
end
pro1=0.2;
pro2=0.8;
beta=0.5;
M1=[0.1 -0.3]';
M2=-0.02;
M3=0.01;
N1=[0.1 0.01];
N2=0.02;
N3=0.01;
%% e_{i}>0,i=1,2,3...13
e1=0.1;
e2=0.05;
e3=0.02;
e4=0.03;
e5=0.12;
e6=0.02;
e7=0.01;
e8=0.5;
e9=0.1;
e10=0.02;
e11=0.02;
e12=0.03;
e13=0.12;

e2=0.5;
e3=0.2;
e4=0.3;
e5=0.2;
e6=0.2;
e7=0.1;
e8=0.5;
e9=0.1;
e10=0.2;
e11=0.2;
e12=0.3;
e13=0.2;


%% Variables assignment
[m,n]=size(B{1})
r=n+1;
m1=m*2;
R=0.1*eye(m1);

%% constant value
T1=0.01;
I1=0.02;
T2=1/T1;
I2=1/I1;
%% SDP variables

S1=sdpvar(m1);
S2=sdpvar(m1);

W1=sdpvar(m1);
W2=sdpvar(m1);

P1=sdpvar(m);
P2=sdpvar(m);
% P1=cell(1,2);
% for i=1:2
%     P1{i}=sdpvar(m);
% end

% for i=1:4
%     Q{i}=sdpvar(m);
% end


%% performance index
Gamma=sdpvar(n);

%% event-triggered scheme parameters
Lam1=sdpvar(1);
%Lam1=sdpvar(1);


%% filter variables and P
Af=cell(2,2);
Bf=cell(2,2);
Cf=cell(2,2);
for i=1:2
    for j=1:2
            Af{i}{j}=sdpvar(m);
            Bf{i}{j}=sdpvar(m,n);
            Cf{i}{j}=sdpvar(n,m);
    end
end

s1=R^(1/2)*S1*R^(1/2);
s2=R^(1/2)*S2*R^(1/2);
w1=R^(1/2)*W1*R^(1/2);
w2=R^(1/2)*W2*R^(1/2);
% s1=S1;
% s2=S2;
% w1=W1;
% w2=W2;
Q=[P1 -P2;-P2 P2];
%[-pro1* e1;pro1*Bf{1}*C{1} e1];
e=zeros(m);
e1_2=zeros(1,m);
%% Theorem 2 LMIs
flag=0;
Phi1=cell(2,2);
for i=1:2
    for j=1:2
            %%%%%i表示对象的下标取值，j表示滤波器下标取值；
            %定义
            Phi1{i}{j}=blkvar;
            
            PA=[P1*A{i} -Af{i}{j};-P2*A{i} Af{i}{j}];
            
            PB1=[-pro1*Bf{i}{j}*C{i} e;pro1*Bf{i}{j}*C{i} e];
            
            PB2=[-pro2*Bf{i}{j}*C{i} e;pro2*Bf{i}{j}*C{i} e];
            
            PB3=[-pro2*Bf{i}{j};pro2*Bf{i}{j}];
            
            PB4=[P1*B{i};-P2*B{i}];
            
            QM1=[-Af{i}{j}*M1;Af{i}{j}*M1];
            QM2=[-Bf{i}{j}*M2;Bf{i}{j}*M2];
            QM3=QM2;
            QM4=QM3;
            
            Naa=[e e;e N1'*N1];
            %  N_1=[e N1]'*[e N1];
            N_1=[e e; N1'*N1 e];
            N_3=[N2*C{i} e1_2]'*[N2*C{i} e1_2];
            N_5=[e1_2 -N3*Cf{i}{j}];
            Eij=[L{i}';-Cf{i}{j}'];
            
            Cij=[delta*C{i}'*Lam1*C{i} e;e e];
            
            
            %列举
            Phi1{i}{j}(1,1)=PA+PA'+s1+s2-beta*Q-T2*w1-I2*w2+(1/e1+5*T1/e5+5*T2/e9)*N_1;
            
%             Phi1{i}{j}(1,2)=PB1+T2*w1; 
%             Phi1{i}{j}(1,4)=I2*w2+PB2;
%             Phi1{i}{j}(1,6)=-PB3;
%             Phi1{i}{j}(1,7)=PB4;
%             
%             Phi1{i}{j}(1,8)=T1*PA';
%             Phi1{i}{j}(1,9)=I1*PA';
%             Phi1{i}{j}(1,10)=Eij;
  

            Phi1{i}{j}(1,2)=PB1+T2*W1;
            Phi1{i}{j}(1,9)=PB2+I2*W2;
            Phi1{i}{j}(1,14)=-PB3;
            Phi1{i}{j}(1,19)=PB4;
            
            Phi1{i}{j}(2,2)=-(1-k1)*s1-2*T2*w1;
            
            Phi1{i}{j}(2,3)=sqrt(5*T1)*PB1';
            
            Phi1{i}{j}(3,3)=-Q-Q'+W1;
            Phi1{i}{j}(3,4)=QM3;
            Phi1{i}{j}(4,4)=-1/e7;
            Phi1{i}{j}(2,5)=sqrt(5*I1)*PB1';
            Phi1{i}{j}(5,5)=-Q-Q'+W2;
            Phi1{i}{j}(5,6)=QM3;
            Phi1{i}{j}(6,6)=-1/e10;
            Phi1{i}{j}(2,7)=[N2*C{i} e1_2]';
            Phi1{i}{j}(7,7)=-sqrt(1/e2+5*T1/e6+5*I1/e10);
            
            Phi1{i}{j}(8,8)=-T2*W1;
            
            Phi1{i}{j}(9,9)=-(1-k2)*s2-2*I2*w2+sqrt(1/e3+5*I1/e7+5*I1/e11)*N_3'*N_3+Cij;
            
            Phi1{i}{j}(9,10)=sqrt(5*T1)*PB2';
            Phi1{i}{j}(10,10)=-Q-Q'+W1;
            Phi1{i}{j}(10,11)=QM3;
            Phi1{i}{j}(11,11)=-1/e7;
            Phi1{i}{j}(9,12)=sqrt(5*I1)*PB2';
            Phi1{i}{j}(12,12)=-Q-Q'+W2;
            Phi1{i}{j}(12,13)=QM3;
            Phi1{i}{j}(13,13)=-1/e11;
          
            
           
            Phi1{i}{j}(14,14)=-Lam1;
            Phi1{i}{j}(14,15)=sqrt(5*T1)*PB3';
            Phi1{i}{j}(15,15)=-Q-Q'+W1;
            Phi1{i}{j}(15,16)=QM4;
            Phi1{i}{j}(16,16)=-1/e8;
            Phi1{i}{j}(14,17)=sqrt(5*I1)*PB3';
            Phi1{i}{j}(17,17)=-Q-Q'+W2;
            Phi1{i}{j}(17,18)=QM4;
            Phi1{i}{j}(18,18)=-1/e12;
            
            
            Phi1{i}{j}(19,19)=-Gamma;
            Phi1{i}{j}(19,20)=sqrt(5*T1)*PB4';
            Phi1{i}{j}(20,20)=-Q-Q'+W1;
            Phi1{i}{j}(19,21)=sqrt(5*I1)*PB4';
            Phi1{i}{j}(21,21)=-Q-Q'+W2;
            
            
            Phi1{i}{j}(1,22)=e1*QM1;
            Phi1{i}{j}(22,22)=-e1;
            Phi1{i}{j}(1,23)=e2*QM2;
            Phi1{i}{j}(23,23)=-e2;
            Phi1{i}{j}(1,24)=e3*QM3;
            Phi1{i}{j}(24,24)=-e3;
            Phi1{i}{j}(1,25)=e4*QM4;
            Phi1{i}{j}(25,25)=-e4;
            Phi1{i}{j}(1,26)=sqrt(5*T1)*PA';
            Phi1{i}{j}(26,26)=-Q-Q'+W1;
            Phi1{i}{j}(26,27)=QM1;
            Phi1{i}{j}(27,27)=-1/e5;
            Phi1{i}{j}(1,28)=sqrt(5*I1)*PA';
            Phi1{i}{j}(28,28)=-Q-Q'+W2;
            Phi1{i}{j}(28,29)=QM1;
            Phi1{i}{j}(29,29)=-1/e9;
            Phi1{i}{j}(1,30)=Eij;
            Phi1{i}{j}(30,30)=-1+e13*M3*M3';
            Phi1{i}{j}(1,31)=N_5';
            Phi1{i}{j}(31,31)=-1/e13;
            
            
            Phi1{i}{j}(9,32)=I2*W2;
            Phi1{i}{j}(32,32)=-I2*W2;
            
            Phi1{i}{j}=sdpvar(Phi1{i}{j});
    end
end

% LL=sdpvar(1);
% LL=
%
%constrain=max(eig(Q1))
%% solution of LMIs

LMIs=[Q>0,P1-P2>0,s1>0,s2>0,w1>0,w2>0,0.1>Gamma>0,Lam1>0,Phi1{1}{1}<0,...,
    Phi1{1}{2}<0,Phi1{2}{1}<0,Phi1{2}{2}<0];
%sol=solvesdp(LMIs,Gamma)
if 0
options=sdpsettings('solver','sdpt3','verbose',0);
%sol=solvesdp(LMIs)
sol=optimize(LMIs)
sol.solvertime
sol.yalmiptime
end
if 1
options=sdpsettings('solver','sdpt3','verbose',0);
sol=optimize(LMIs,[],options); 

flag=0
if sol.problem == 0 
    [primal,~]=check(LMIs); 
    flag=min(primal)>=0
else
    yalmiperror(sol.problem) 
end
end
%% optimize
if 0
    options=sdpsettings('solver','sedumi','verbose',0);
    sol=optimize(LMIs,[],options); 

    if sol.problem == 0
            valobj=value(Phi1);
            if 0
                Gamma=check(LMIs);
                if min(Gamma)<1
                    Gamma=Gamma;
                end
            end
    else
        disp('error');
    end

end

%% solve unknown decesion variables
if 1
    
for i=1:2
    for j=1:2
        Af{i}{j}=value(Af{i}{j});
        Bf{i}{j}=value(Bf{i}{j});
        Cf{i}{j}=value(Cf{i}{j});
    end
end

P1=value(P1)
P2=value(P2)
for i=1:2
        AF{i}=inv(P2)*Af{i}{i}
        BF{i}=inv(P2)*Bf{i}{i}
        CF{i}=Cf{i}{i}
end

S1=value(S1)
S2=value(S2)
W1=value(W1)
W2=value(W2)

% 
% S1=R^(-1/2)*value(S1)*R^(-1/2)
% S2=R^(-1/2)*value(S2)*R^(-1/2)
% W1=R^(-1/2)*value(W1)*R^(-1/2)
% W2=R^(-1/2)*value(W2)*R^(-1/2)

P=R^(-1/2)*[P1 -P2;-P2 P2]*R^(-1/2);
% V0=max(eig(P))+T1*exp(beta*T1)*max(eig(S1))+I1*exp(beta*I1)*max(eig(S2))+...,
%     T1*T1*exp(beta*T1)*max(eig(W1))+I1*I1*exp(beta*I1)*max(eig(W2))
Lam1=value(Lam1);
%Lam1=value(Lam1);ph(1,1)=Gamma*d/beta-c2*exp(-beta*T)*min(eig(Q));
Gamma=value(Gamma)
%xlswrite('filter_parameter.xls',Af);
end
d=0.4;
c2=0.2;
T=2;
c1=0.000001;
% ddd=V0*c1+Gamma*d/beta-c2*exp(-beta*T)*min(eig(P))
% index=Gamma*exp(beta*T)
%end