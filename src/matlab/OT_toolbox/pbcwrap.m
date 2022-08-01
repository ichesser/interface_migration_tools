function dp = pbcwrap(dp1,xlo,xhi,ylo,yhi,zlo,zhi)
%%% wraps set of input coordinates %%%%
    dp = dp1; 
    Lx = xhi-xlo; 
    Ly = yhi-ylo; 
    Lz = zhi-zlo; 
    
    L = [Lx Ly Lz];
    hilo = [xlo xhi; ylo yhi; zlo zhi];
    
    n = length(dp(:,1)); 
    
    for i = 1:3
        
        di = dp(:,i); 
        Li = L(i); 
        lo = hilo(i,1); 
        hi = hilo(i,2); 
        
        for j = 1:n
            xj = di(j); 
            
            if (xj < lo)
                while (xj < lo)
                    xj = xj + Li; 
                end
            end
            
            if (xj >= hi)
                while (xj >= hi)
                    xj = xj - Li; 
                end
            end
            
            %make assignment
            dp(j,i) = xj; 
        end
        
    end
end