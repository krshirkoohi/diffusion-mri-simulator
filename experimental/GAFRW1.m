% Two-layer structure

%% params
ra    = [2.5,5.0]*1e-6;
gam   = 2.675e8;
gMax  = 0.5;
delta = 50e-3;
Delta = 100e-3;
tf    = 2*Delta;
diffu = [2.0,2.0,2.0]*1e-9;
perma = [1.0,0]*1e-5;
relax = [Inf Inf];

N_ii  = 2^5;
N_rw  = 1e2; 
N_t   = 1e2; 
t     = linspace(0,tf,N_t);

const = Tables_FRW;

%% substrate
[X,Y] = meshgrid(linspace(-ra(2)*1.3,ra(2)*1.3,N_ii));
d = sqrt(X.^2+Y.^2);
I = uint8(d<=ra(2)) + uint8(d<ra(1));

%% sequence
gSteps = 1e3; 
gRes = gMax/gSteps;
G=0*t; G(t<=delta)=1; G(t>=Delta&t<=Delta+delta)=-1;
dir=[1,0]; G=kron(dir,G');
seq.G = G; seq.t = t; seq.G_s = 0:gRes:gMax;
bVal = (gam^2 * seq.G_s.*seq.G_s * delta^2 * (Delta - delta/3));

%% analytical answer
if ~exist ('Sa','var') % debugging - save time
    ML.d = 2; ML.r = ra; ML.D = [diffu(1), diffu(2)]; ML.W = [perma(1),perma(2)]; ML.T = relax;
    Sa = ML_compute(seq.G_s,delta,Delta,ML);
end

%% simulation with GAFRW

% setup
N_i = size(I);
dim = numel(N_i);
N_c = sum(unique(I(:))>0);
r   = X(1,2)-X(1,1); % for I, not quadtree
dt  = t(2)-t(1);
D   = diffu; 
dx  = sqrt(2*dim*dt*D');
Pt  = 1.0e-5* sqrt((dt*pi)/diffu(2)); % Pt
if N_c>1 
        P = [0 0 0; 0 +1 Pt; 0 Pt +1]; % 11 permitted if ECS used
    else 
        P = NaN; 
end
% throw in particles 
N       = round(N_rw*numel(I)/sum(I(:)>0));
pos     = prodsum(rand(N,dim),r.*N_i);
ind_s   = maskPos(pos,N_i,r);
% delete points in non-diffusion regions
pos(I(ind_s)==0,:)    = [];
ind_s(I(ind_s)==0)  = [];
% total particles remaining
N_p     = size(pos,1);
phase   = zeros(N_p,1);

% fast random walk
Tm = max(t);
%for tt=1:N_t
for n=1:N_p % do for each particle separately 
    Tr = Tm;
    while Tr > 0
        curPos = pos(n,:);
        curComp = I(maskPos(curPos,N_i,r));
        [k,dist] = dsearchn(curPos,[X(I~=curComp),Y(I~=curComp)]);
        dist = min(dist);
        tau = gen_tau(const);
        t_unit = dist^2/D(1);
        dt = t_unit * tau;
        u = randCirc(1,dim);

        if dt < Tr
            Tr = Tr - dt;
        else
            dt = Tr;
            Tr = 0;
            tau = Tr / t_unit;
        end

        % subtend search radius
        posStep = prodsum(randCirc(1,dim),dist);

        ind_s = maskPos(curPos,N_i,r);
        ind_p = maskPos(mod(curPos+posStep,N_i(1)*r),N_i,r); 

        Ie = I(ind_s)==I(ind_p); 
        
        if ~isnan(P(1)) 
            if Ie ~= 1
                a = I(ind_s);
                b = I(ind_p);
                if a>0 && b>0
                    Pt = P(a+1,b+1);
                    prob = rand(1,1);
                    if prob <= Pt
                        Ie = 1;
                    end
                end
            end
        end

        newPos = curPos + posStep;
        %phase(n) = phase(n) + sum(prodsum(newPos,G(tt,:)),2);
        dph_nl = (evalPhi1(const, tau)* dist * t_unit * u);
        phase(n)=phase(n)+sum(prodsum(dt * newPos, dph_nl));

%         M1=rem(newPos,N_i(1)*r);
%         M2=mod(newPos,N_i(1)*r);
%         ind_bc=(M2~=newPos).*sign(M1);
%         % (!) NPA read up
%         dph_nl = (evalPhi1(const, tau)* dist * t_unit);% * u);
%         phase(n)=phase(n)+(dt * newPos + dph_nl);
%         %sum(prodsum(-ind_bc*N_i(1)*r,sum(G(1:tt,:))),2);
%         newPos=M2;      

    end
end
%end

% compute signal
if ~isfield(seq,'G_s'),seq.G_s=1;end
S=zeros(length(seq.G_s),1);
for gg=1:length(S),S(gg)=mean(exp(1j*gam*dt*phase*seq.G_s(gg)));end

%% plot
figure('Position', [200 200 2000 800]),
%subplot(1,3,1),pcolor(X,Y,I),axis image,shading interp,
subplot(1,3,2),semilogy(bVal,real(S)),hold on,semilogy(bVal,Sa,'r-'), %axis([0, 2.5e11, 1e-2,1])
legend('monte carlo','grebenkov'),title('Results')
ylabel('Signal attenuation'),xlabel('B-Value') 
subplot(1,3,3),plot(bVal,abs((Sa-real(S'))./Sa)),title('relative error'), %axis([0, 2.5e11, 1e-2,1])
annotation('textbox',[.0 1.0 .0 .0], ...
    'String',sprintf('N_{RW}=%d, N_{ii}=%d, N_t=%d, N_{G_s}=%d',N_rw,N_ii,length(t),length(seq.G_s)),'EdgeColor','none')


%% functions (makeshift)
function tau = gen_tau(const)
    U = rand;
    while U == 1
        U = rand;
    end
    if U < 1e-5
        A = 2*log(U*sqrt(pi)/2);
        t = 1; % initial guess;
        f = t*log(t)+A*t+1/2;
        err = abs(f);
        df = log(t)+1+A;
        while err > eps
            t = t - f/df;
            f = t*log(t)+A*t+1/2;
            err = abs(f);
            df = log(t)+1+A;
        end
        tau = t;
    elseif U < 0.99
        tau = interp1(const.F_tau(:,2),const.F_tau(:,1),U);
    else
        tau = log(2/(1-U))/(pi^2);
    end
end
function const = Tables_FRW()
    load('PhysicalConstants.mat','D','gamma')
    load('Ftau.mat','newFtau')
    load('Phi1.mat','Phi1');
    load('Phi2.mat','Phi2');
    load('R(t).mat','Rtable');
    const.D = D;
    const.gam = gamma;
    const.F_tau = newFtau(2:135,:);
    const.Phi_1 = Phi1(1:800,:);
    const.Phi_2 = Phi2;
    const.R_rms = Rtable(1:200,:);
end
function output = evalPhi1(const, tau)
    if tau < 0.05
        output = tau*(1-2*tau)/2;
    elseif tau < 0.8
        output = interp1(const.Phi_1(:,1),const.Phi_1(:,2),tau);
    else
        output = 0.75/pi^2;
    end
end
function output = evalPhi2(const, tau)
    if tau < 1e-3
        output = sqrt(2*tau^3/3);
    else
        output = interp1(const.Phi_2(:,1),const.Phi_2(:,2),tau);
    end
end

