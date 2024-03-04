function k = freq(n,L);
    
     N1 =  floor((n-1)/2);
     N2 = (n/2)*ones(rem(n+1,2));
     freq = [(0:N1)  N2 (-N1:-1)]';
     k  = 2*pi/L * freq;
