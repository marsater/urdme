%   This file simulates a quorum sensing model in 2D and
%   generates a movie of the simulation that can be saved.

% Based on the PDE model from the article:
% ["A synchronized quorum of genetic clocks", nature 2010, isbn={0028-0836}]

% Reference:
%    [1] D. Mårsäter, N. Axelsson:
%    "Computational modelling of quorum sensing using cascade delay",
%    Independent project in Engineering Physics, Uppsala university (2022).

close all; clear all; tic

% Create 100x1000 rectangular domain
gd = [3;4;0;1000;1000;0;0;0;100;100]; %
ns = [82;49];
sf = 'R1';

G = decsg(gd,sf,ns);

% Distance between points in mesh
h = 10;

% Create the mesh
[P,E,T] = initmesh(G,'hmax',h);

figure(1), clf,
pdegplot(G), axis equal

figure(2), clf,
pdemesh(P,E,T,'NodeLabels','on'), axis tight, axis equal

% Parameters
CA = 1; CI = 4; delta = 1e-03; alfa = 2500;
k = 1; k1 = 0.1; b = 0.06; mu = 1;
gammaA = 15; gammaI = 24; gammaH = 0.01;
f = 0.3; g = 0.01; d0 = 0.88;
D = 2.5; d  = 0.8; dd04 = (d/d0)^4;

% Diffusion rate for external AHL
D1 = 800;

% Discretize time (approx. simulation time for Tend = 1400: 35-45 minutes)
Tend = 1400;
Ntime = 10*Tend;
tspan = linspace(0,Tend,Ntime+1);

% Delay time
tau = 10;

% Define reactions for species A,I,Hi,He
RA1 = 'A > (gammaA * A) / (1+f*(A/vol+I/vol)) > @ ';
RI1 = 'I > (gammaI * I) / (1+f*(A/vol+I/vol)) > @ ';
RHi1 = 'Hi > (gammaH*Hi*A/vol)/(1+g*A/vol)  > @';
RHe1 = 'He > mu*He > @';
RHi2 = '@ > (b*I)/(1+k*I/vol) > Hi' ;

Diff_out = 'Hi > D*Hi > He';
Diff_in = 'He > (d/(1-d))*D*He > Hi';

% Define increase rate for A & I (used in cascade)
rateAI = '((1-dd04)*(delta+alfa*(Hi*Hi/vol/vol))/(1+k1*Hi*Hi/vol/vol))*vol';

% Number of reactions in cascade
n = 5;

% Number of species ('n' for species in cascade, +4 for A,I,Hi,He)
Nspecies = n + 4;

% Create cascade species
Xn = sprintf('X%d',n);

species = {'A' 'I' 'Hi' 'He' 'X$i' Xn};

% Diffusion rates
diffusions = zeros(1,Nspecies);
diffusions(1,4) = D1;
diffusions = num2cell(diffusions);
umod = pde2urdme(P,T,diffusions);
umod.sd = ceil(umod.sd);

% Initial data
A0 = 0;
I0 = 0;
Hi0 = 0;
He0 = 4*mean(umod.vol);

rates = {'r' n/tau ... % cascade rate
        'CA' CA 'CI' CI 'delta' delta 'alfa' alfa 'tau' tau 'k' k ...
        'k1' k1 'b' b 'gammaA' gammaA 'gammaI' gammaI 'gammaH' gammaH ...
        'f' f 'g' g 'd0' d0 'D' D 'd' d ...
        'D1' D1 'mu' mu ...
        'dd04' dd04'};

umod = rparse(umod,{RA1 RI1 RHi1 RHi2 RHe1 Diff_in Diff_out ...
    ['@ > ' rateAI ' > X1'] ...         % Reaction before cascade
    'X$i > r*X$i > X$e' ...             % Cascade
    [Xn ' > r*' Xn ' > A +I+I+I+I'] ... % last step of cascade
    }, ...
    species, rates,'danino_2D.c',{['ie']' [1:n-1; 2:n]});

% Set initial data in voxel 1, rest are empty
voxelMid = repmat( [A0; I0; Hi0; He0; zeros(n,1);], 1 , 1 );
initialData = zeros(Nspecies,numel(umod.vol));
initialData(:,1) = voxelMid;

% Solve
umod = urdme(umod,'tspan',tspan, ...
    'u0',initialData, ...
    'solver','nsm', ...
    'seed',0);
toc

% Plot each species over time for voxels 1-5
 uu = reshape(umod.U,Nspecies,numel(umod.vol),[]);
 for i = 1:5
    figure(i)
    set(gcf,'color','w');

    plot(tspan, squeeze(uu(1,i,:)), tspan, squeeze(uu(2,i,:)), ...
        tspan, squeeze(uu(3,i,:)), tspan, squeeze(uu(4,i,:)) );

    xlabel('Time t');
    legend('A (AiiA)', 'I (LuxI)', 'Hi (Internal AHL)', 'He (External AHL)');
    title("Plot for voxel "+ i )
    grid on
    hold off
 end

% visualize using PDE Toolbox
umod = urdme2pde(umod);

%%
% Generate movie
%====================================================================
% Initialization
workspace;

% Save every xth frame, video plays faster this way
jumpframe = 10;
i = 1:jumpframe:length(tspan);
numberOfFrames = length(i);

hFigure = figure;
hFigure.Position;

% Size of movie plot
hFigure.Position(3:4) = [1000 200];

% Where to show it on your screen (changes depending on screen size)
hFigure.Position(1:2) = [360/2 278];

% Preallocate movie, which will be an array of structures.
% First get a cell array with all the frames.
allTheFrames = cell(numberOfFrames,1);

% standard height and width
%   vidHeight = 344;
%   vidWidth = 446;

% Adjusted for 1000x100 mesh size
vidHeight = 100;
vidWidth = 1000;

allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};

% Next get a cell array with all the colormaps.
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};

% Now combine these to make the array of structures.
myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);

% Need to change from the default renderer to zbuffer to get it
% to work right. openGL doesn't work and Painters is too slow.
set(gcf, 'renderer', 'zbuffer');

%========================================================================
% Create the movie.

% After this loop starts, BE SURE NOT TO RESIZE THE WINDOW AS IT'S
% SHOWING THE FRAMES, or else you won't be able to save it.
for i = i
	cla reset;

    % Enlarge figure to full screen.
 	% set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);

    set(gcf,'color','w');
    pdesurf(umod.pde.P,umod.pde.T,umod.pde.U(2,:,i)');
	axis('tight')
    view(0,90)
    colormap('parula')
    caption = sprintf('Frame #%d of %d, t = %.1f', i, length(tspan), tspan(i));
	title(caption, 'FontSize', 10);
    xlabel('X [um]'), ylabel('Y [um]')
	drawnow;

    % Works, but doesn't capture axes
    %thisFrame = getframe(gca);

    % To capture labels and axes
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;

    rect = [-ti(1), -ti(2)+2, pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    thisFrame = getframe(ax,rect);

	% Write this frame out to a new video file.
	myMovie(i) = thisFrame;
end
% close(writerObj);

%=======================================================================
% See if they want to replay the movie.
message = sprintf('Done creating movie\nDo you want to play it?');
button = questdlg(message, 'Continue?', 'Yes', 'No', 'Yes');
drawnow;	% Refresh screen to get rid of dialog box remnants.
close(hFigure);
if strcmpi(button, 'Yes')
	hFigure = figure;
    hFigure.Position;
    hFigure.Position(3:4) = [1000 200];
    hFigure.Position(1:2) = [360/2 278];

    % Get rid of extra set of axes that it makes for some reason.
	axis off;

	% Play the movie.
    set(gcf,'color','w');
    movie(myMovie(1:jumpframe:end));
	close(hFigure);
end

%=======================================================================
% See if they want to save the movie to an avi file on disk.
promptMessage = sprintf('Do you want to save this movie to disk? (use mp4)');
titleBarCaption = 'Continue?';
button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
if strcmpi(button, 'yes')
	% Get the name of the file that the user wants to save.
	% Note, if you're saving an image you can use imsave()
    % instead of uiputfile().
	startingFolder = pwd;
	defaultFileName = {'*.avi';'*.mp4';'*.mj2'};
	[baseFileName, folder] = uiputfile(defaultFileName, 'Specify a file');
	if baseFileName == 0
		% User clicked the Cancel button.
		return;
	end
	fullFileName = fullfile(folder, baseFileName);

	% Create a video writer object with that file name.
	% The VideoWriter object must have a profile input argument,
    % otherwise you get jpg.

	% Determine the format the user specified:
	[folder, baseFileName, ext] = fileparts(fullFileName);
	switch lower(ext)
		case '.jp2'
			profile = 'Archival';
		case '.mp4'
			profile = 'MPEG-4';
		otherwise
			% Either avi or some other invalid extension.
			profile = 'Uncompressed AVI';
	end
	writerObj = VideoWriter(fullFileName, profile);
	open(writerObj);
    numberOfFrames = length(myMovie);
	for frameNumber = 1 : jumpframe : numberOfFrames
	   writeVideo(writerObj, myMovie(frameNumber));
    end

	close(writerObj);
	% Display the current folder panel so they can see their newly created file.
	cd(folder);
	filebrowser;
	message = sprintf('Finished creating movie file\n     %s.\n\nDone with demo!', fullFileName);
	uiwait(helpdlg(message));
else
	uiwait(helpdlg('Done!'));
end
