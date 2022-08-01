function d = pbcdist(dp,Lx,Ly,Lz)
   %%% d: displacement vectors, vectors along rows
   %%% Lx,Ly,Lz: sim box lengths
   %%% output: displacement vectors with PBC's
   d = dp; 
   n = length(d(1,:)); 
   L = [Lx Ly Lz]; 
   for i = 1:n
       Li = L(i); 
       %check magnitude of displacement values for reduction 
       booli_lo = (d(:,i) < -0.5*Li);
       booli_hi = (d(:,i) >= 0.5*Li);
       
       d(booli_lo,i) = d(booli_lo,i) + Li; 
       d(booli_hi,i) = d(booli_hi,i) - Li; 
       
   end
   
   

end
