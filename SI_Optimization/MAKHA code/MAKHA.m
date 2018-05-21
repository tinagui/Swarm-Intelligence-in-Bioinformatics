%% A Hybrid Algorithm Monkey-Krill Herd Hybrid Algorithm(MKHHA or MAKHA) 
%Conducted by Ahmed M.E. Khalil and Seif Eddeen K. Fateen
%Last edited on 30 October 2014
function [u,fval,NumEval,MinVector,NFEval] = MAKHA(cost,Lb,Ub)
global IPMIT NP NVAR
format long;
seed=sum(100*clock);
randn('state', seed)

%%Nomenclature of basic code inputs and important outputs%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost:     Fitness function
%Lb:       Lower boundary of decision variable
%Ub:       Upper boundary of decision variable
%IPMIT:    Maximum number of iterations
%u:        Optimum global minimum solution
%fval:     Estimated global minimum value
%NumEval:  Number of function evaluations
%NVAR:     Number of variables

nsbest=0; 
uo=Lb+rand*(Ub-Lb);
MaxGeneration=IPMIT;
Y=abs((Lb(1)-Ub(1))/2); % R average boundary
mbest=3E50;

%%Calculation section
[u1,fval1,NumEval1,MinVector1,NFEval1,XGBEST1] = MAKHAm(cost,Lb,Ub,NP,Y,uo,mbest);


%% Find the Global
if fval1 < mbest
    mbest = fval1;
    nsbest = u1;
end
u=nsbest; fval=mbest;

% Total number of function evaluations
NumEval=NumEval1;

% NFEval
NFEval=NFEval1;
MinVector=MinVector1;
XGBEST=XGBEST1;

end

%% A Mixed Algorithm (MA-KHA)
function [u,fval,NumEval,MinVector,NFEval,XGBEST] = MAKHAm(cost,Lb,Ub,NP,Y,nsbest,mbest)
global IPMIT NVAR
format long;

Kbest=1E50;
seed=sum(100*clock);
randn('state', seed)

% Initial random guess
u0=Lb+(Ub-Lb).*rand(size(Lb));

% Calcualte dimension
d=length(u0);
d=NVAR;

% generating the initial locations of n Hybrids
ns=init(NP,d,Lb,Ub,u0); %old

 
%Chaotic Search methods
cm=0.1;

% Initial values of an array
mm=ones(NP,1)*10^100;
u=nsbest;
MI=1;  %The loops of foraging of the hybrid (keep it equals 1)

fval=mbest;
MaxGeneration=IPMIT;
MinVector=zeros(MaxGeneration,1);
NFEval=zeros(MaxGeneration,1);
XGBEST=zeros(MaxGeneration,NVAR);
k1=0;
D1=1; C_flag =1;
Kf=3E50;F=zeros(d,NP);
Xf=zeros(d,IPMIT);


% Total number of function evaluations
NumEval=MI*IPMIT*NP+MI*IPMIT+2*NP*IPMIT;


% Iterations or pseudo time marching
for k=1:MaxGeneration,     %%%%% start iterations
    
   if k>2
       if rand < 0.075/Y
           for i=1:NP
               for l=1:d
                   cm=chaos_Ricker(cm);
                   while cm<0||cm>1
                       cm=rand;
                       cm=chaos_Logistic(cm);
                   end
                   ns(i,l)=Lb(l)+(Ub(l)-Lb(l))*cm;
               end
           end
       end  
    end
    
    
    ns=findlimitsRand(NP,ns,Lb,Ub);
    

% Evaluate new solutions (for all n Hybrids)
for l=1:NP
   mm(l)=Fun(cost,ns(l,:));
end

% Ranking Hybrids' positions by their Mountain top value/objectives
[mass,Index]=sort(mm);
ns_tmp=ns;

for i=1:NP,
 nsg(i,:)=ns_tmp(Index(i),:);
 mm1(i)=mass(i);
end

%% Find the local current best
nbest=ns_tmp(Index(1),:); mgood=mass(1);

%% Find the Global
if mgood< mbest
   mbest=mgood;
   nsbest=nbest;
end
u=nsbest; fval=mbest; 

% The Watch-Jump process
[ns2,mg]= MA_WatchJump2(nsg,NP,d,mm1,cost,Lb,Ub,u,Y,k,IPMIT);   

[fval1, A] = min (mg);
u1= ns2(A,:);
% Find the Global
if fval1 < mbest
    mbest = fval1;
    nsbest = u1;
end
u=nsbest; fval=mbest;


[mass,Index]=sort(mg);
ns_tmp=ns2;

for i=1:NP,
 ns2(i,:)=ns_tmp(Index(i),:);
 mg(i)=mass(i);
end

if k==1
    Xb=ns2;
    Kb=mg;
end

for l=1:NP
    if Kb(l)<mg(l)
        Kb(l)=mg(l);
        Xb(l,:)=ns2(l,:);
    end
end

[Kb,Index]=sort(Kb);
Xb_tmp=Xb;

for i=1:NP
    Xb(i,:)=Xb_tmp(Index(i),:);
end


%Hybrid (Krill) herd process (Food Forage)
[ns,u2,fval2,Xib,Kib,Kf,F,Xf]=KHA_Food(ns2',mg,Lb,Ub,fval,u',NVAR,NP,cost,k,MaxGeneration,Xb',Kb,Kf,F,Xf,D1,C_flag,MI);

%% Find the Global
if fval2 < mbest
    mbest = fval2;
    nsbest = u2;
end
u=nsbest; fval=mbest; 
MinVector(k)= mbest;

% Ranking Hybrids' positions by their Mountain top value/objectives
[MM,Index]=sort(Kib);
ns_tmp=Xib';

for i=1:NP,
 ns(i,:)=ns_tmp(Index(i),:);
end

for l=1:NP
    if Kb(l)<MM(l)
        Kb(l)=MM(l);
        Xb(l,:)=ns(l,:);
    end
end

[Kb,Index]=sort(Kb);
Xb_tmp=Xb;

for i=1:NP
 Xb(i,:)=Xb_tmp(Index(i),:);
end

if k==k1 || k1==0
   somersault=1;
   c=-NVAR*Y/10; f=NVAR*Y/10; 
elseif k <(k1+IPMIT/10)
   somersault=2;
    c=-NVAR*Y/10; f=NVAR*Y/10;
elseif k <(k1+round(2*IPMIT/10))
    somersault=3;
     c=-NVAR*Y/10; f=NVAR*Y/10;
end

%The Somersault
if somersault==1
ns=MA_Somersault(ns,NP,c,f,Lb,Ub,u,Y);
end

if somersault==2
ns=MA_Somersault2(ns,NP,c,f,Lb,Ub,u,Y);
end

XGBEST(k,:)=u;
MinVector(k)= fval;

if k==(k1+round(2*IPMIT/10))
    k1=k1+round(2*IPMIT/10);
end
    
IFE=(MI+2)*k*NP+MI*k; 
NFEval(k)=IFE;   
end
end

% -------------------------------------------------------
% ----- All the subfunctions are listed here ------------
% The initial locations of n masses
function [ns]=init(n,d,Lb,Ub,u0)
  % if there are bounds/limits,
if length(Lb)>0,
   for i=1:n,
   ns(i,:)=Lb+(Ub-Lb).*rand(1,d);
   end
else
   % generate solutions around the random guess
   for i=1:n,
   ns(i,:)=u0+randn(1,d);
   end
end
end


% KHA_Foraging motion
function [X,u4,fval4,Xib,Kib,Kf,F,Xf]= KHA_Food(X,K,LB,UB,Kgb,Xgb,NP,NK,cost,nr,MG,Xib,Kib,Kf,F,Xf,D1,C_flag,MI)
Dt = mean(abs(UB-LB))/2;
Y=abs((LB(1)-UB(1))/2);
MIF=MI; D1=0; C_flag=1;
Vf=0.02;  %Optimized parameters
if nr< round(MG/1.8) 
nr=1;
MG=MI;
Xib=X;
Kib=K;
elseif nr<MG
Xib=X;
Kib=K;
end
NN=NK;
Dmax=0.002+(0.01-0.002)*rand; %The same as KH-based % Optimized parameters 
for j = 1:MI
    if j==MIF
        D1=0; C_flag=1;
    end
     for ll = 1:NP;
            Sf(ll) = (sum(X(ll,:)./K));
     end
     Xf(:,nr) = Sf./(sum(1./K)); %Food Location
     Xf(:,nr) =Kfindlimits(Xf(:,nr)',LB,UB,Xgb'); %Bounds Checking
     Kf(nr) = Fun(cost,Xf(:,nr)');
      
        if 2<=nr
            if Kf(nr-1)<Kf(nr)
                Xf(:,nr) = Xf(:,nr-1);
                Kf(nr) = Kf(nr-1);
            end
        end

        Kft(MI)=Kf(nr);
        Xft(:,MI)=Xf(:,nr);
        
        if MI>1
            if Kft(MI-1)<Kft(MI)
                Xf(:,nr) = Xft(:,MI-1);
                Kf(nr) = Kft(MI-1);
            end
        end
        
        [Kw] = max(K);
       
        if rand < 0.3
        if round(rand)==0
        [Kgb,A]=min(K);
         Xgb =X(:,A);
        end
        end

         Kw_Kgb = Kw-Kgb;
         w = 0.1+0.8*(1-nr/MG);
          for i = 1:NK
          %Calculation of distances
            Rf = Xf(:,nr)-X(:,i);
             % % % % % % % % % % % % % Foraging Motion % % % % % % % % % %
            % Calculation of FOOD attraction
             if Kf(nr) < K(i)
                Beta_f=-2*(1-nr/MG)*(Kf(nr) - K(i)) /Kw_Kgb/ sqrt(sum(Rf.*Rf)) * Rf;
            else
                Beta_f=0;
             end
             
            % Calculation of BEST position attraction
            Rib = Xib(:,i)-X(:,i);
            if Kib(i) < K(i)
                Beta_b=-(Kib(i) - K(i)) /Kw_Kgb/ sqrt(sum(Rib.*Rib)) *Rib;
            else
                Beta_b=0;
            end
            % Foraging Motion
            F(:,i) = w*F(:,i)+Vf*(Beta_b+Beta_f);
             % % % % % % % % % % % % % Physical Diffusion % % % % % % % % %
            if D1==1 
            D = Dmax*(1-nr/MG)*floor(rand+(K(i)-Kgb)/Kw_Kgb)*(2*rand(NP,1)-ones(NP,1));
            else
            D=0;
            end
            D1=Dmax*(1-nr/MG);
            DX = Dt*(F(:,i)+D);
                            
     % % % % % % % % % % % % % Crossover % % % % % % % % % % % % %
      if C_flag ==1
          p=randperm(NN,1);
          q=randperm(NN,1);
          if p==i
              p=randperm(NN,1);
          end
          if q==i
              q=randperm(NN,1);
          end
     
             Mu=  0.9 + 0.05*(K(i)-Kgb)/Kw_Kgb;
             C_rate = 0.9 + 0.2*(K(i)-Kgb)/Kw_Kgb;
            Cr = rand;
            if Cr < C_rate ;
                % Random selection of Krill No. for Crossover
                NK4Cr = round(NK*rand+.5);
                % Crossover scheme
                X(:,i)=X(:,NK4Cr).*(1-Cr)+X(:,i).*Cr;
               
                if rand < Mu
                    mue=rand;
                    X(:,i)= Xgb+ mue.*(X(:,p)-X(:,q));
                   
                end
            end
      end
           % Update the position
            X(:,i)=X(:,i)+DX;
            X(:,i)=Kfindlimits(X(:,i)',LB,UB,Xgb');
            X(:,i)=findlimitsRand(1,X(:,i)',LB,UB);
            K(i)=Fun(cost,X(:,i)');
            if K(i)<Kib(i)
                Kib(i)=K(i);
                Xib(:,i)=X(:,i);
            end
            if MI>=1                         
            [Kgb, B]= min(Kib);
            Xgb=Xib(:,B);
            K(i)=Kib(i);
            X(i)=Xib(i);
            end
          end
end
        [fval4, A]= min(Kib);
        u4=Xib(:,A)';
        X=Xib'; 
end


% The Somersault (without abs)
    function S=MA_Somersault(ns,n,c,f,Lb,Ub,u,Y)
     p=(sum(ns))/n;
    ceta=c+rand()*(f-c);
     S=ns;
     d=size(ns,2);
     for i=1:n
         for j=1:d
             S(i,j)=ns(i,j)+ceta*abs(p(j)-ns(i,j));
         end
     end
     [S]=findlimitsRand(n,S,Lb,Ub);
    end
      

function [ns]=Kfindlimits(ns,Lb,Ub,best)
% Evolutionary Boundary Constraint Handling Scheme
n=size(ns,1);
for i=1:n
    ns_tmp=ns(i,:);
    I=ns_tmp<Lb;
    J=ns_tmp>Ub;
    A=rand;
    ns_tmp(I)=A*Lb(I)+(1-A)*best(I);
    B=rand;
    ns_tmp(J)=B*Ub(J)+(1-B)*best(J);
    m=sum(isnan(ns_tmp));
  if m>1 || m==1
      ns_tmp=Lb+(Ub-Lb).*rand(size(Lb));
  end
  ns(i,:)=ns_tmp;
end
end

% Randomizing
function [ns]=findlimitsRand(n,ns,Lb,Ub)
 for i=1:n,
     % Apply the lower bound
  ns_tmp=ns(i,:);
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I)+rand*(Ub(I)-Lb(I));

  % Apply the upper bounds
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J)-rand*(Ub(J)-Lb(J));
  % Update this new move
  ns(i,:)=ns_tmp;
  m=sum(isnan(ns_tmp));
  if m>1 || m==1
      ns_tmp=Lb+(Ub-Lb).*rand(size(Lb));
  end
  ns(i,:)=ns_tmp;
end
end


function K=Fun(cost,u)
K=feval(cost,u);
end

function [ns]=findlimits(n,ns,Lb,Ub)
for i=1:n,
     % Apply the lower bound
  ns_tmp=ns(i,:);
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);

  % Apply the upper bounds
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move
  ns(i,:)=ns_tmp;
end
end

% The Somersault2(with abs)
function S=MA_Somersault2(ns,n,c,f,Lb,Ub,u,Y)
fact=(1)/(n-1); sum=0; diff=0; ceta=c+rand()*(f-c);
for l=1:n
    for i=1:n
        diff=ns(l,:)-ns(i,:);
        sum=sum+diff;
    end
end
p=fact.*sum;
%p=(max(ns)+min(ns))/2;
S=ns;
d=size(ns,2);
for i=1:n
    for j=1:d
        S(i,j)=ns(i,j)+ceta*abs(p(j)-ns(i,j));
    end
end
[S]=findlimitsRand(n,S,Lb,Ub);
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % The Watch-Jump process
    function [ns,zn,V1]=MA_WatchJump2(ns,n,d,zn,cost,Lb,Ub,u,Y,k,IPMIT)
    nsR=ns;znR=zn;ns1=ns; Mut=1;
    Kw_Kgb=zn(n)-zn(1); u=ns(1,:);
       for i=1:n
        Rj=randperm(d,d);
        Rjj=randperm(d,d);
        Rk=randperm(n,n);
        ind=randperm(n,1);
        Ki=randperm(n,n);
        p=randperm(n,1);
        q=randperm(n,1);
        for j=1:d
        m=-1+rand*2;
        while m==0
        m=-1+rand*2;
        end
            while Rk(ind)==i
                Rk(ind)=randperm(n,1);
            end
            m1=1;
        nsR(i,j)= ns(i,Rj(j))+m*m1*(ns(i,Rj(j))-ns(Rk(ind),Rj(j)));       
        end
       end
       if Mut==1
       Mu1= 0.9 + 0.05*(zn(i)-zn(1))/Kw_Kgb;
       Mu2= 0.6 + 0.05*(zn(i)-zn(1))/Kw_Kgb;
       for i=1:n
           p=randperm(n,1);
           q=randperm(n,1);
           if p==q
             q=randperm(n,1);
           end
       if rand<Mu1
            mue=rand;
           nsR(i,:)= u+ mue.*(ns(p,:)-ns(q,:));       
       end
        if rand<Mu2
            mue=rand;
           nsR(i,:)= u+ mue.*(nsR(p,:)-nsR(q,:));       
        end
       end
       
       end
       nsR=findlimitsRand(n,nsR,Lb,Ub);
       srt=1;
        for i=1:n
            znR(i)=Fun(cost,nsR(i,:));
        end
    if srt==1
    [znR1,Index]=sort(znR);
    ns_tmp=nsR;
    
    for i=1:n
      nsR(i,:)=ns_tmp(Index(i),:);
       znR(i)=znR1(i);
    end
    end
    
    for i=1:n
        if znR(i)<(zn(i))
            ns(i,:)=nsR(i,:);
            zn(i)=(znR(i));
        end      
    end
    end
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   
   
   
   