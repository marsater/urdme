% Written as part of a Bsc project in Engineering Physics 
% by David Mårsäter and Nils Axelsson, spring 2022
%
% Generates time-space plots for quantities of Luxi(I) or the other
% proteins (AiiA (A), internal AHL (Hi) or external AHL (He)) in a
% quorum sensing network.
%
% Based on PDE model from the article: 
%["A synchronized quorum of genetic clocks", nature 2010, isbn={0028-0836}]


% Model parameters
CA = 1; CI = 4; delta = 1e-03; alfa = 2500; k = 1; k1 = 0.1; 
b = 0.06; gammaA = 15; gammaI = 24; gammaH = 0.01; f = 0.3; 
g = 0.01; d0 = 0.88; D = 2.5; d  = 0.8; mu = 1; dd04 = (d/d0)^4;
D1 = 800; % Diffusion rate for external AHL

tau = 10; % Delay Time

% Discretize time [0,Tend]
Tend = 500; %  
Ntime = 10*Tend;
tspan = linspace(0,Tend,Ntime+1);

% Space discretisation
h = 5; % um Distance between cells (incorporated in diffusion)
Nvoxels = 100; % Nvoxels coupled (diffusion)
L = 5*Nvoxels;   % um length of cell array
volume = d*h^3;  % Total intracellular volume for cells in voxel

% Rate of increase for A and I, used with cascad delay
rateAI = ['((1-dd04) *(delta + alfa *' ...
         '(Hi*Hi/vol/vol))/(1+k1*Hi*Hi/vol/vol))*vol'];

% Def reactions
RA1 = 'A > (gammaA * A) / (1+f*(A/vol+I/vol)) > @ ';    % A decrease

RI1 = 'I > (gammaI * I) / (1+f*(A/vol+I/vol)) > @ ';    % I decrease

RHi1 = 'Hi > (gammaH*Hi*A/vol)/(1+g*A/vol)  > @';       % Hi decrease
RHi2 = '@ > (b*I)/(1+k*I/vol) > Hi' ;                   % Hi increase

RHe1 = 'He > mu*He > @';                                % He decrease
Hiout = 'Hi > D*Hi > He';             % Diffusion out of cell
Hein = 'He > (d/(1-d))*D*He > Hi';    % Diffusion into cell

n = 5; % Number of intermediate reactions in cascade delay      

Nspecies = n +4; % n for cascade, 4 for real substances

% Create diffusion matrix
e = ones(Nvoxels,1);
DD = spdiags([e -2*e e],-1:1,Nvoxels,Nvoxels);
% Periodic boundary conditions
DD(1,end) = DD(1,end)+e(1);
DD(end,1) = DD(end,1)+e(end);

ixR = 4; % Index of species He
DD = kron(DD,sparse(ixR,ixR,1,Nspecies,Nspecies));
DD = DD .* D1./(h^2);

% initial data
A0 = 0;
I0 = 0;  
Hi0 = 0;
He0 = 5*volume; % Iniital for oscillation startup

% Def species
Xn = sprintf('X%d',n);
species = {'A' 'I' 'Hi' 'He' 'X$i' Xn};

% Def rates
rates = {'r' n/tau ... % For cascade delay
        'delta' delta 'alfa' alfa 'tau' tau 'k' k 'k1' k1 ...
        'b' b 'gammaA' gammaA 'gammaI' gammaI 'gammaH' gammaH ...
        'f' f 'g' g 'D' D 'd' d 'mu' mu 'dd04' dd04'};

vmod = rparse([],{RA1 RI1 RHi1 RHi2 RHe1 Hiout Hein ... % Reactions above
    ['@ > ' rateAI ' > X1'] ... % Do a delayed increase for A & I
    'X$i > r*X$i > X$e' ...     % Cascade delay         
    [Xn ' > r*' Xn ' > A + I+I+I+I'] ... % Last step 
    }, ...
    species, rates,'danino_diffusion.c',{['ie']' [1:n-1; 2:n]});

% Put initial data in the middle voxel
voxeldataMid = [A0; I0; Hi0; He0; zeros(n,1);];
initialData = zeros(Nspecies,Nvoxels);
initialData(:,floor(Nvoxels/2)) = voxeldataMid;

% Solve
vmod = urdme(vmod,'tspan',tspan, ...
    'u0',initialData, ...
    'D', DD, ...
    'vol',ones(1,Nvoxels).*volume, ...
    'sd',ones(1,Nvoxels), ...
    'solver','nsm', ...
    'seed', 0);

uu = reshape(vmod.U,Nspecies,Nvoxels,[]);

i = floor(Nvoxels/2); % Plot species for middle voxel
figure(1)
set(gcf,'color','w');
plot(tspan, squeeze(uu(1,i,:)), tspan, squeeze(uu(2,i,:)), ...
     tspan,squeeze(uu(3,i,:)), tspan, squeeze(uu(4,i,:)) );
ylabel('Number of molecules');
xlabel('Time t');
legend('A (AiiA)', 'I (LuxI)', 'Hi (Internal AHL)', 'He (External AHL)');
title("Plot for voxel "+ i )
 

% Space-time plot 
figure(6)
set(gcf,'color','w');
hold on
for i = 1:Nvoxels
    scatter(i.*ones(1,length(tspan(1:30:end))), tspan(1:30:end), 40, ...
            squeeze(uu(2,i,1:30:end))', 'filled','s')
end
colormap('parula')
ylabel('Time [min]', 'FontSize',16)
xlabel('Voxel', 'FontSize',16)
title("Initial data in voxel " + floor(Nvoxels/2));