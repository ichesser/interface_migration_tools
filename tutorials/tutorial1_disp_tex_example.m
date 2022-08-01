
% Example 1: plot a net displacement field from an MD simulation as an inverse pole figure in different reference frames

startup_mtex %feel free to comment this out after MTEX has been loaded once

%% DEFINE SEVERAL CONVENTIONS 

% CONVENTIONS FOR THIS EXAMPLE: The left grain (grain 1) and the right grain (grain 2) are each described by three
% integer vectors.  The first vector on the left is the direction in
% fcc (cubic) coordinates (for the left grain) of the sample x axis.  The x axis
% is the normal to the grain boundary, so this vector gives the grain
% boundary normal.  The second and third vectors similarly give the
% directions of the sample y and z axes.  The three vectors on the right
% give the same information for the right grain.  

O1 = [1 0 0; 0 3 1; 0 -1 3]; 
O2 = [1 0 0; 0 3 -1; 0 1 3]; 
O1n = normr(O1); O2n = normr(O2); %grain rotation matrices, should satisfy det(O) = 1

% DEFINE MATERIAL DEPENDENT CONSTANTS: 
lattice_parameter = 3.52; %angstroms, for Ni in this case
r0 = lattice_parameter*sqrt(2)/2; %nearest neighbor distance in perfect FCC crystal


%% IMPORT NET DISPLACEMENT FIELD DATA 

% The data must contain atomic coordinates and displacement vectors relative to an initial state
% *** example data: displacements swept by Sigma 5 (100) twist GB, T = 800 K, ramped SDF method *** 

data_path = '../data/net_displacement_fields/GB2/example_displacement_field.xyz';
X0 = importdata(data_path,' ',2); 
X = array2table(X0.data,'VariableNames',{'id','type','x','y','z','dx','dy','dz'}); %modify this to fit your data file

%% DEFINE REFERENCE FRAMES AND PLOT DISPLACEMENT MAGNITUDES

disps = [X.dx X.dy X.dz]; 
disp_magnitudes = vecnorm(disps,2,2); 

%by default, we are working in the translated dichromatic pattern since MD data has a microscopic shift arising from GB energy minimization
dTDP = disps; 
ndisps = length(dTDP);

% assuming zero coupling factor (valid for this data), we can move from TDP  --> SDP (shuffling dichromatic pattern) by subtracting average net displacement
ave_net_disp = mean(dTDP,1); 
dSDP = dTDP-ave_net_disp; 

%now, we analyze the displacement orientations in the bicrystal sample frame (S) and the crystal frame (C1) of grain 1 in the TDP and SDP 
nTDP_S = normr(dTDP); nTDP_C1 = nTDP_S*O1n; 
nSDP_S = normr(dSDP); nSDP_C1 = nSDP_S*O1n;
dmagTDP = vecnorm(dTDP,2,2); %displacement magnitudes, these should be the same as earlier
dmagSDP = vecnorm(dSDP,2,2); 

% below, we plot displacement magnitude histogram: notice that most 
% displacements associated with migration are smaller than the 1NN distance
% in the FCC crystal.Displacements in SDP are smaller and
% bimodality is enhanced. 

figure(1) 
hTDP = histogram((dmagTDP)/r0,'Normalization','Probability','FaceColor','r','FaceAlpha',0.5);
hold on
hSDP = histogram((dmagSDP)/r0,'Normalization','Probability','FaceColor','b','FaceAlpha',0.5);
legend([hTDP,hSDP],{'TDP','SDP'})
ylabel('Probability','FontSize',18)
xlabel('Displacement length (normalized by 1NN distance)','FontSize',18)

%% NOW, WE DEFINE DISPLACEMENT POLE FIGURES REPRESENTING SHUFFLING VECTOR ORIENTATION DISTRIBUTIONS

%%%% NEED MTEX (startup_mtex) for the lines below this point!
vtex_TDP_S = vector3d(nTDP_S'); vtex_TDP_C1 = vector3d(nTDP_C1'); 
vtex_SDP_S = vector3d(nSDP_S'); vtex_SDP_C1 = vector3d(nSDP_C1'); 

%%%% ADDITIONAL WEIGHTS CAN BE ADDED TO POLE FIGURES 
w_uniform = ones(1,ndisps); %these weights do nothing
w_TDP = dmagTDP.^2; w_SDP = w_uniform; %by conventions of paper, in TDP weight displacements by squared magnitude, but not in SDP

%%%% displacement pole figures as densities on S^2 
density_TDP_S = calcDensity(vtex_TDP_S,'weights',w_TDP); %spherical harmonic fits
density_TDP_C1 = calcDensity(vtex_TDP_C1,'weights',w_TDP);

density_SDP_S = calcDensity(vtex_SDP_S,'weights',w_SDP); 
density_SDP_C1 = calcDensity(vtex_SDP_C1,'weights',w_SDP);

%% PLOT DISPLACEMENT POLE FIGURES AND IDENTIFY RELEVANT SHUFFLE DIRECTIONS 

rotalign = rotation.byAxisAngle(vector3d([0 0 1]),90*degree); %rotation convention from paper

% In the sample frame, 4 prominent displacement poles are apparent representing four shuffle orientations
% the distribution in a vertical line represents twist type motion (see OVITO file for visualization) 

figure
density_plot = density_SDP_S;
pC1 = plot(rotate(density_plot,rotalign),'FontSize',20,'TR',[]); 
mtexColorMap white2black %LaboTeX

% In the crystal frame, we can define reference directions to identify specific shuffle orientations

cs432 = crystalSymmetry('432'); 
hklcustom = [1 0 3; -3 0 1; 3 0 -1; -1 0 -3];
hklvec0uniq = unique(hklcustom,'rows'); 
hklvec = vector3d(hklvec0uniq');
mHKL = Miller(hklvec,cs432,'uvw');

figure
density_plot = density_SDP_C1;
pC1 = plot(rotate(density_plot,rotalign),'FontSize',20,'TR',[]); 
hold on
annotate(mHKL,'labeled','Color',rgb('dark grey'),'FontSize',20)
mtexColorMap LaboTeX

%% FOR VISUALIZATION, ASSIGN COLORS TO EACH DISPLACEMENT ORIENTATION

cs1 = crystalSymmetry('1'); %crystal symmetry for first color scheme
sstcmap1 = HSVDirectionKey(cs1); 

vtouse = vtex_TDP_S; %can modify this choice 
vcolors = sstcmap1.direction2color(vtouse); %one can add these columns to the XYZ file and color the displacement vectors in OVITO 
                                            %see ../ovito/displacement_field_with_color.ovito
                                            