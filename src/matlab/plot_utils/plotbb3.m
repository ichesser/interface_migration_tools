function BB = plotbb3(P,K,a3,b3,sigma)

Bx = randn(P-1,K)/sqrt(P); 
By = randn(P-1,K)/sqrt(P); 
Bz = randn(P-1,K)/sqrt(P);

B3 = zeros(P,K,3); 
B3(2:end,:,1) = Bx; B3(2:end,:,2) = By; B3(2:end,:,3) = Bz; 

t = linspace(0,1,P)';

abline = (1-t).*a3 + t.*b3; 
alph = 0.8; 
lw = 1.5; 
c0 = [1 0 0]; c1 = [0 0 1];
%c0 = rgb('red'); c1 = rgb('blue'); 
BB = B3; 
for k = 1:K 
    
    randnoise = squeeze(B3(:,k,:));  
    cumnoise = cumsum(randnoise); 
    bb = cumnoise - t.*cumnoise(end,:); %construct each bridge by adding noise to abline 
    bb = sigma*bb; 
    lab = abline + bb;
    
    BB(:,k,:) = lab; 
end

% PLOT via P to access scalar color for each of k bridges
colvec = (1-t)*c0 + t*c1; %color interpolation

for p = 1:(P-1)
    
    c = colvec(p,:); 
    
    %%%% find k coordinates for a given interpolating value
    bbk = (BB(p:p+1,:,:)); 
    
    xp = squeeze(bbk(:,:,1));
    yp = squeeze(bbk(:,:,2)); 
    zp = squeeze(bbk(:,:,3)); 
    
    patch(xp,yp,zp,'r','EdgeColor',c, 'EdgeAlpha',alph,'FaceColor','none','Linesmoothing','on','LineWidth',lw);
    hold on
end


end