% EXAMPLE 2: example application of the forward model to predicting migration mechanisms of the Sigma 5 (100) twist GB 

matpath = '../src/matlab/'; 
addpath(matpath) 
addpath([matpath,'OT_toolbox/'])
addpath([matpath,'plot_utils/'])

% DEFINE MATERIAL DEPENDENT CONSTANTS: 
lattice_parameter = 3.52; %angstroms, for Ni in this case
r0 = lattice_parameter*sqrt(2)/2; %nearest neighbor distance in perfect FCC crystal

%% The starting point for this tutorial is a unit cell of the translated dichromatic pattern (TDP) in the sample frame of the bicrystal
% Assumptions for this example: no shear coupling (disconnections of pure step character only), microscopic shift is known  
% For GB ID 2 from the Olmsted survey, the energy optimal microscopic shift has been measured to be p = [0.32,0.25,-0.125];
% This measurement was performed on the 0 K GB structure relative to the Sigma 5 coherent dichromatic pattern (CDP) with p = [0,0,0].

%%% the following data file contains an orthogonal repeat unit of the TDP for FCC Ni with two atom types, one for each orientation %%%
TDP_path = '../data/dichromatic_patterns/GB2/prism_GB2_zero_shear_u0_0_p_default.data';
TDPdata0 = importdata(TDP_path,' ',18);
TDPdata = TDPdata0.data;
TDPtable = array2table(TDPdata,'VariableNames',{'id','type','x','y','z'}); 

%%% read SIM BOX DATA %%%
simboxcell = TDPdata0.textdata(4:6); 
simbox = zeros(3,2); 
for j = 1:3
    splitur = strsplit(simboxcell{j},' ');
    simbox(j,:) = str2double(splitur(1:2));
end
simboxcell2 = TDPdata0.textdata(7); 
for j = 1
    splitur = strsplit(simboxcell2{j},' ');
    offdiag = str2double(splitur(1:3));
end
xlounit = simbox(1,1); xhiunit = simbox(1,2); 
ylounit = simbox(2,1); yhiunit = simbox(2,2); 
zlounit = simbox(3,1); zhiunit = simbox(3,2);
xy = offdiag(1); xz = offdiag(2); yz = offdiag(3); 
Lxunit = xhiunit-xlounit; Lyunit = yhiunit-ylounit; Lzunit = zhiunit-zlounit; 
%%% X,Y coordinate data, where X is crystal 1 and Y is crystal 2%%%
N = length(TDPdata); 
types = TDPdata(:,2); 
Xbasis = TDPdata(types==1,3:5)'; 
Ybasis = TDPdata(types==2,3:5)';

%% APPLY regularized optimal transport algorithm modified from Peyre et. al (Computational Optimal Transport, DOI: 10.48550)

epsset = 0.05; %regularization parameter: play with this value (in range 0.005 and above) to access different path probability distributions
                % small values (0.005 is a good example) give min-shuffle mapping
p = [0.32,0.25,-0.1250]; %precomputed, set to [0,0,0] if unknown and results will match TDP

%%% visualization settings %%%
frame = 'SDP'; %choose 'TDP', 'SDP' or 'CDP': see information for each reference frame in papers
viewx = 90; %adjust these to rotate final figure
viewy = 0; 
usearrows = true; usebridge = false; %only one option should be selected here, second uses brownian bridges as in Peyre et. al, but is slow
sphereplt=true; showcoincidence = true; pbcon = true; 
addshear = false; beta_xz = 0; u0_set = 0; %these options are only relevant if shear coupling factor is nonzero
Kcutoff = 1e-8; %filter out paths with probabilities smaller than this value

            
nX = length(X); nY = length(Y); 
if nX ~= nY
    disp('WARNING: one to one mapping is not possible in dichromatic pattern')
end

%%%% Construct pairwise distance matrix with PBC's
pbcon = true; 
x = X'; y = Y'; %coordinates along rows
nx = length(x); N = nx;
c = zeros(nx,nx);  
for i = 1:nx
    xi = x(i,:); 
    dveci = y-xi; 
    dvecpbci = pbcdist(dveci,Lx,Ly,Lz); %displacement vectors considering PBC's
    di = vecnorm(dvecpbci'); %NOTE: vecnorm operates on columns
    c(i,:) = di; %can save time by not double counting here
end
if pbcon
    c = c.^2; 
else
    c2 = distmat(X,Y).^2;
    c = c2;
end
%%%
% Run Sinkhorn Algorithm

cutoff = 1/(N^2); 
epsilon = epsset; 
mu = ones(N,1)/N; nu = mu;
options.niter = 5000;
tic
[u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);
toc
XYcell{iep,ish} = [x y]; %cols: x1,x2,x3,y1,y2,y3
costcell{iep,ish} = c; 
gammacell{iep,ish} = gamma; 

%%% OUTPUT DISPLACEMENTS in TDP, CDP, SDP 
toosmall = find(gamma < cutoff); 
gammamod = gamma;
gammamod(toosmall) = 0; 
[I,J,NP] = find(gammamod); %find nonzero indices (row,col,entry)
maxnp = max(max(NP)); %set line thickness

%%% COLLECT DISPLACEMENT VECTORS in TDP FROM GAMMA 
dispvecs = zeros(length(I),10); %displacement vectors, weights, coordA,coordB

for s = 1:length(I)
    K = NP(s); %path probability;
    acoords = [X(1,I(s)) X(2,I(s)) X(3,I(s))];
    bcoords = [Y(1,J(s)) Y(2,J(s)) Y(3,J(s))];

    %%% find vector connecting a to b
    dab = bcoords-acoords; 
    if pbcon
        dpab = pbcdist(dab,Lx,Ly,Lz); %assumes orthogonal box
    end
    newb_coords = acoords + dpab; 

    dispvecs(s,:) = [dpab K acoords newb_coords]; 
end

%%% CLASSIFY DIFFERENT REFERENCE FRAMES %%%
ndisps = length(dispvecs(:,1)); 
Kvec = dispvecs(:,4); %column vector
Kvecrenorm = normr(Kvec'); %row vector, renormalized K values 

Xcoords = dispvecs(:,5:7); Ycoords = dispvecs(:,8:10); 
dTDP = dispvecs(:,1:3); XTDP = Xcoords; YTDP = Ycoords; 
dCDP = dTDP - p; XCDP = Xcoords; YCDP = Ycoords-p; 
Dvec = sum(dTDP.*Kvec); %probabilistic expression for total net displacement/atom
dSDP = dTDP-Dvec; XSDP = Xcoords; YSDP = Ycoords-Dvec; 

Mvecest = Dvec-p; 
microvecs = [p Mvecest Dvec]; 


%%% TABULATE # UNIQUE POLES 

[upolelist,~,repid] = unique(round(dTDP,2),'rows'); 
npoles = length(upolelist); 

cumprobs = zeros(npoles,3); 
for i = 1:npoles
    idorig = find(repid==i); 
    multiplicity = length(idorig); 
    probsum = sum(Kvec(idorig));
    aveprob = probsum/multiplicity; 
    cumprobs(i,1) = aveprob; 
    cumprobs(i,2) = multiplicity;
    cumprobs(i,3) = probsum; 
end

upolelistcum = [upolelist(:,1:3) cumprobs]; 
[~,sortids] = sort(cumprobs(:,3));
upolelistcum = flip(upolelistcum(sortids,:));

dispmags = vecnorm(upolelistcum(:,1:3),2,2);
poleprobs = upolelistcum(:,end); 

aveL2cost = (dot(dispmags,poleprobs))^2; 

lentol = 0.05; %tolerance for coincident site
nshortpoles = sum(dispmags < lentol); 
nuniqpoles = npoles-nshortpoles; 

disp('unique poles,TDP: dx,dy,dz,aveprob,multiplicity,probsum:')
disp(upolelistcum)
disp('L2 cost:')
disp(aveL2cost)
disp('#unique poles:')
disp(nuniqpoles)
disp('net displacement/atom:')
disp(Dvec)



if sphereplt
    SPHERE_RES = 20; % resolution for your spheres
    SPHERE_RAD = 0.15; % radius of spheres
    [xb, yb, zb] = sphere(SPHERE_RES);
end

if strcmp(frame,'TDP')
    Xp = XTDP; Yp = YTDP; 
elseif strcmp(frame,'CDP')
    Xp = XTDP; Yp = YCDP; 
else 
    Xp = XTDP; Yp = YSDP; 
end


ndisps = length(dTDP);
gammaflat = gamma(:);

figure(1)

plotcube([Lx Ly Lz],[xlo ylo zlo],0,rgb('black'))
hold on

 Xs = []; Ys = []; %copies of coordinates to overwrite with sheared coordinates
k = 1; 
for i = 1:ndisps
    %%% DEFINE EACH DISPLACEMENT FROM COORDINATES
    pathprob = Kvec(i); 
    if pathprob > Kcutoff %don't draw very improbable paths
        acoords = Xp(i,:); bcoords = Yp(i,:); 

        if addshear %account for nonzero FGB (subtract out shear)
            FGB1 = -beta_xz/2; FGB2 = beta_xz/2; 
            u0str2 = num2str(u0set); %0.23 for CDP, GB 1 
            dz_C1 = (acoords(1)-u0set(1))*FGB1; %correction in z based on x distance from u0
            dz_C2 = (bcoords(1)-u0set(1))*FGB2; 
            acoords = acoords - [0 0 dz_C1]; 
            bcoords = bcoords - [0 0 dz_C2]; 
            acoords = pbcwrap(acoords,xlo,xhi,ylo,yhi,zlo,zhi); 
            bcoords = pbcwrap(bcoords,xlo,xhi,ylo,yhi,zlo,zhi);
        end

        if pbcon
            dpab = pbcdist(bcoords-acoords,Lx,Ly,Lz); 
        else
            dpab = bcoords-acoords; 
        end

        newb_coords = acoords+dpab; 
        Xs(k,:) = acoords; 
        Ys(k,:) = newb_coords; 
        k = k+1; 

        if usearrows
            arrow = mArrow3(acoords,newb_coords,'color','k','stemWidth',pathprob/(maxnp/5)/100,'tipWidth',0.07);
            hold on
        end

        if usebridge      
            nbridges = ceil(pathprob*1000); %number of bridges
            nsteps = 250; %number of steps for brownian bridge 
            sigma = .2*sqrt(eps);
            if eps<.01
                sigma=0;
            end
            BB = plotbb3(nsteps,nbridges,acoords,newb_coords,sigma);
            hold on
        end
    end
end

%%%%%%%%% POINT SET PLOTTING ROUTINE %%%%%%% 
if addshear 
    xp1 = Xs(:,1); xp2 = Xs(:,2); xp3 = Xs(:,3); 
    yp1 = Ys(:,1); yp2 = Ys(:,2); yp3 = Ys(:,3); 
else 
    xp1 = Xs(:,1); xp2 = Xs(:,2); xp3 = Xs(:,3); 
    yp1 = Ys(:,1); yp2 = Ys(:,2); yp3 = Ys(:,3);   
end

npts = length(xp1);

if ~sphereplt
    p0 = plot3(xp1,xp2,xp3,'ko');
    p0.MarkerSize=15; 
    p0.MarkerFaceColor=rgb('light red'); 
    hold on

    p1 = plot3(yp1,yp2,yp3,'ko');
    p1.MarkerSize=15; 
    p1.MarkerFaceColor=rgb('light blue'); 
    hold on
else
    for i = 1:npts
        surf(SPHERE_RAD*xb+xp1(i), SPHERE_RAD*yb+xp2(i), SPHERE_RAD*zb+xp3(i), 'facecolor', rgb('light red'), 'edgealpha', 0);
        hold on
        surf(SPHERE_RAD*xb+yp1(i), SPHERE_RAD*yb+yp2(i), SPHERE_RAD*zb+yp3(i), 'facecolor', rgb('light blue'), 'edgealpha', 0);
        hold on
    end
end

if showcoincidence
    Xtest = [xp1 xp2 xp3]; Ytest = pbcwrap([yp1 yp2 yp3],xlo,xhi,ylo,yhi,zlo,zhi);
    Idx = rangesearch(Ytest,Xtest,0.1); %find all Y points within distance 0.05 of X, i.e. coincident sites
    nid = length(Idx); 
    coincident_ids_Y = [];
    k = 1;
    for i = 1:nid
        idlist = Idx{i};
        if ~isempty(idlist)
            coincident_ids_Y(k) = idlist(1); 
            k = k+1; 
        end
    end

    if ~isempty(coincident_ids_Y)
        if ~sphereplt
            p2 = plot3(yp1(coincident_ids_Y),yp2(coincident_ids_Y),yp3(coincident_ids_Y),'ko');
            p2.MarkerSize=15; 
            p2.MarkerFaceColor=rgb('white');
            hold on
        else
            for i = coincident_ids_Y
                SPHERE_RADC = SPHERE_RAD+0.01; 
                surf(SPHERE_RADC*xb+yp1(i), SPHERE_RADC*yb+yp2(i), SPHERE_RADC*zb+yp3(i), 'facecolor', rgb('tan'), 'edgealpha', 0);
                hold on
            end
        end
    end
end

axis off; axis equal; 

if sphereplt %set lighting
    light;
    lighting gouraud;
    material([0.7 0.9 0.8]); %Set the ambient (direct), diffuse (skylight), and specular (highlight) surface reflectivity
end

view(viewx,viewy);
viewstr = ['view','_',num2str(viewx),'_',num2str(viewy)];

            
            
            
