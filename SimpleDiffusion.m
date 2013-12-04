function [varargout] = SimpleDiffusion(varargin)
format compact;format short;close all;
scsz = get(0,'ScreenSize');



% USER ENTERED VALUES
NSteps = 400;			% number of steps (loop time)
Ndots = 10;			% number of particles
Scale = 1/10;			% scale of model (life:model)
TimeStep = 1000;		% time step (ms)
Delay = .01;			% slow animation rate by X%
XYZSize = [10 20];		% size of enclosure (XY in µm)
DiffRateA = 0.2;		% diffusion rate coefficient A
DiffRateB = 0.1;		% diffusion rate coefficient B


% BASE DIFFUSION RATES EQUATIONS
Sc = Scale;		% scale of model (life:model)
t = TimeStep/1000;	% time step (ms)
dm = 2;			% dimensions
Da = DiffRateA*t/Sc;	% Diffusion Rate A (D = L² / 2d*t)
Db = DiffRateB*t/Sc;	% Diffusion Rate B
Dr = Da/Db;		% Ratio of Da:Ds (1/Ls)^2;
Dn = Da/Dr;		% new D after scaling L
k = sqrt(dm*Da);	% stdev of D's step size distribution
L = sqrt(2*dm*Da);	% average diagonal (2D) step size
Lx = L/sqrt(2);		% average linear (1D) step size
Ls = 1/sqrt(Dr);	% scales Lx values for Dn
MSD = 2*dm*Da;		% mean squared displacement


XYZL = zeros(3,Ndots);		% XYZ particle locations
XYZS = zeros(3,Ndots);		% XYZ step sizes

XYZSize = XYZSize./Scale;	% scale enclosures 
XWIDE = XYZSize(1)/2;		% half X enclosure size (will double below)
YHIGH = XYZSize(2)/2;		% half X enclosure size (will double below)
BOX = [-1 2 2 2] ./Scale;	% special area location [X Y W H]


%===============================%
for Nt = 1:NSteps 

	XYZS = (k * randn(3,Ndots));	% generates step sizes
	XYZL = XYZL+XYZS;		% adds step to location
	XYZLp(:,Nt) = XYZL(:,1);	% save step of first dot (for trace)

	% Keep everything inside enclosure  %
	[XYZL] = ENCLOSE(Nt,XYZL,XWIDE,YHIGH,Ndots);

	% Plot live diffusion
	MAINPLOT(Nt,XYZL,XWIDE,YHIGH,BOX);
	TRACEDOT(Nt,XYZLp,XWIDE,YHIGH,BOX);

pause(Delay);
if mod(Nt,100)==0;Nt 
end;
end % for Nt = 1:Nsteps 
%===============================%



varargout = {XYZLp}; % export traced particle location
end % end main function
%=====================================================%




function [XYZL] = ENCLOSE(Nt,XYZL,XWIDE,YHIGH,Ndots)

	for j = 1:Ndots 
		if XYZL(1,j)>(XWIDE) || XYZL(1,j)<(-XWIDE)
			XYZL(1,j) = sign(XYZL(1,j))*(XWIDE);
		elseif XYZL(2,j)>(YHIGH) || XYZL(2,j)<(-YHIGH)
			XYZL(2,j) = sign(XYZL(2,j))*(YHIGH);
		end
	end

end


%-------------------------------------------%
% PLOT Particle Motion
%-------------------------------------------%
function [] = MAINPLOT(Nt,XYZL,XWIDE,YHIGH,BOX)
%-------------------------------------------%
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH];
%---

figure(1)
hold off;
subplot(1,2,1), 
scatter(XYZL(1,:),XYZL(2,:),5,[0 0 1]);
axis([xlim, ylim]);
rectangle('Position',BOX)
hold off;

end




%-------------------------------------------%
% TRACEDOT
%-------------------------------------------%
function [] = TRACEDOT(Nt,XYZLp,XWIDE,YHIGH,BOX)
%-------------------------------------------%
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH];
%---

figure(1)
subplot(1,2,2), 
plot( XYZLp(1,:), XYZLp(2,:) );
axis([xlim, ylim]);
rectangle('Position',BOX)
hold on;
%---

end


%=================================%
%           3D PLOT
%---------------------------------%
function [] = PLOT3D(Nt,XYZL,XWIDE,YHIGH)
%-------------------------------------------%
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH];
%---

figure(1)
hold off;
subplot(1,2,1), 
scatter3(XYZL(1,:),XYZL(2,:),5,[0 0 1]);
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

end





