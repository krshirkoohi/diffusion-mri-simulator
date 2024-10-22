% Two-layer structure

%% params

% physical
ra    = [2.5,5.0]*1e-6;
gam   = 2.675e8;
gMax  = 0.5;
delta = 50e-3;
Delta = 100e-3;
tf    = 2*Delta;
diffu = [2.0,2.0,2.0]*1e-9;
perma = [1.0,0]*1e-5;
relax = [Inf Inf];

% scan
N_ii  = 200;
N_rw  = 1e5;
N_t   = 1e3;
U     = N_rw * N_t;
t     = linspace(0,tf,N_t);

%% substrate

[X,Y] = meshgrid(linspace(-ra(2)*1.3,ra(2)*1.3,N_ii));
d = sqrt(X.^2+Y.^2);
I = uint8(d<=ra(2)) + uint8(d<ra(1));

%% sequence

gRes = 0.001; gSteps = gMax/gRes;
G=0*t; G(t<=delta)=1; G(t>=Delta&t<=Delta+delta)=-1;
dir=[1,0]; G=kron(dir,G');
seq.G = G; seq.t = t; seq.G_s = 0:gRes:gMax;
bVal = (gam^2 * seq.G_s.*seq.G_s * delta^2 * (Delta - delta/3));

%% analytical answer

if ~exist ('Sa','var') % debugging - save time
    ML.d = 2; ML.r = ra; ML.D = [diffu(1), diffu(2)]; ML.W = [perma(1),perma(2)]; ML.T = relax;
    Sa = ML_compute(seq.G_s,delta,Delta,ML);
end

%% simulation

tic
simu.r = X(1,2)-X(1,1);
simu.N = N_rw;
simu.D = diffu;
simu.P = flip([ML.W; flip(ML.W)]);
S = diffSim(I,seq,simu,gam);
elapsedTime = toc;
save(sprintf('Results/Sim1/Nii=%d Nrw=1e%d Nt=1e%d',N_ii,log10(N_rw),log10(N_t))); % save                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          