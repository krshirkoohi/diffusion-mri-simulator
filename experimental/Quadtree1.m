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
N_rw  = 1e3; 
N_t   = 1e3; 
t     = linspace(0,tf,N_t);

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

%% quadtree
% approximate isotropic regions and update relevant tickers
% provides shortcuts to take during a random walk
close all;

% simulation parameters
dim = 2;
r   = X(1,2)-X(1,1); % parent resolution
dt  = t(2)-t(1);
Pt  = perma(1)* sqrt((dt*pi)/diffu(2)); % Pt
D   = diffu; 
P   = [-1 -1 -1; -1 +1 Pt; -1 Pt +1]; % 11 permitted if ECS used
N_i = size(I);
N_p = N_rw;

% quadtree and properties
Q.quadtree = qtdecomp(I,0,2); % indexed by pos in I
[x,y,d_ii] = find(Q.quadtree); blocks = [x+(d_ii/2),y+(d_ii/2),d_ii];
Q.topLeft  = [x,y];
Q.count    = length(blocks);
Q.centre   = [blocks(:,1) blocks(:,2)]; 
Q.dx       = d_ii;
Q.ind      = maskPos(Q.centre*r,N_i,r);
Q.mask     = I(Q.ind); 

%% simulation
% model as a network
Q.adj = zeros(Q.count); % -1 = not connected (0 reserved for permeability)
for a=1:Q.count
    x0 = Q.topLeft(a,1);x1 = Q.topLeft(a,1)+Q.dx(a);
    y0 = Q.topLeft(a,2);y1 = Q.topLeft(a,2)+Q.dx(a);
    boundingBox = [[x0 y0];[x1 y0];[x0 y1];[x1 y1]];
    for b=1:Q.count % get all nodes b connected to a
        if a~=b
            if      all(Q.topLeft(b,:)==boundingBox(1,:)) || ...
                    all(Q.topLeft(b,:)==boundingBox(2,:)) || ...
                    all(Q.topLeft(b,:)==boundingBox(3,:)) || ...
                    all(Q.topLeft(b,:)==boundingBox(4,:))
                Q.adj(a,b) = P(Q.mask(b)+1,Q.mask(a)+1); % weight edges with permeability
                Q.adj(b,a) = P(Q.mask(a)+1,Q.mask(b)+1); 
            end
        end
    end
end

% radius search algorithm with variable search size
%T = graph(Q.quadtree); plot(T);
%A = adjacency(T,'weighted'); [a1,a2,a3] = find(A);

% throw in particles at permitted positions
availPos = Q.centre(Q.mask~=1,:);

% select from availPos, N_p times
ind = zeros(N_p,2);
pos = zeros(N_p,2);
for ii=1:N_p 
    seed = randi(length(availPos),1);
    ind(ii,:) = [availPos(seed,1),availPos(seed,2)];
    pos(ii,1) = X(1,ind(ii,1)); 
    pos(ii,2) = X(2,ind(ii,2));
end
startPos = pos;
phase = zeros(N_p,1);
moves = zeros(N_p,2);
% the random walk
% skip corresponding dt and dx based on box size
for n=1:N_p
    ii = 1;
    tt = 1;
    while tt < length(t) - max(Q.dx)
        oldPos = pos(n,:); 
        oldInd = ind(n,:);
        oldNode = find(Q.centre(:,1)==oldInd(1) & Q.centre(:,2)==oldInd(2));
        jumpSize = Q.dx(oldNode);

        neighbours = findNeighbours(oldNode,Q.adj);
        newNode = neighbours(randi(length(neighbours),1));
        
        newInd = Q.ind(newNode);
        newPos(1) = X(newInd); newPos(2) = Y(newInd);

        Pt = Q.adj(oldNode,newNode);
        prob = rand(1,1);
        if prob<=Pt 
            pos(n,1) = newPos(1); 
            pos(n,2) = newPos(2);
        end

        %moves(n,1) = moves(n,1) + pos(n,1);
        %moves(n,2) = moves(n,2) + pos(n,2);

        tt = tt + jumpSize;
        ii = ii + 1; % number of total steps taken

        displ = [(startPos(n,1) - pos(n,1)),(startPos(n,2) - pos(n,2))];
        
        phase(n) = phase(n) + jumpSize*sum(prodsum(displ,G(tt,:)),2);
    end    
end

%% compute signal
if ~isfield(seq,'G_s'),seq.G_s=1;end
S=zeros(length(seq.G_s),1);
for gg=1:length(S),S(gg)=mean(exp(1j*gam*dt*phase*seq.G_s(gg)));end

%% plot results
figure('Position', [200 200 2000 800]),
%subplot(1,3,1),pcolor(X,Y,I),axis image,shading interp,
subplot(1,3,2),semilogy(bVal,real(S)),hold on,semilogy(bVal,Sa,'r-'), %axis([0, 2.5e11, 1e-2,1])
legend('monte carlo','grebenkov'),title('Results')
ylabel('Signal attenuation'),xlabel('B-Value') 
subplot(1,3,3),plot(bVal,abs((Sa-real(S'))./Sa)),title('relative error'), %axis([0, 2.5e11, 1e-2,1])
annotation('textbox',[.0 1.0 .0 .0], ...
    'String',sprintf('N_{RW}=%d, N_{ii}=%d, N_t=%d, N_{G_s}=%d',N_rw,N_ii,length(t),length(seq.G_s)),'EdgeColor','none')








