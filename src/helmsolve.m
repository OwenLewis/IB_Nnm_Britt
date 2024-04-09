%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Helmoltz equation of the form (a-b*L)*u = f
%
function u = helmsolve(f,a,b,grid)

  switch grid.bcx
    case 'per'  %%% I can expliot alising if I want to 
      k = 0:grid.Nx-1;
      wx = 2*pi*k/grid.Lx;
      transformx = @(v)(fft(v));
      itransformx = @(v)(ifft(v));
    case 'dir'
      k = 1:grid.Nx;
      wx = pi*k/grid.Lx;
      transformx = @(v)(dst_cc(v));
      itransformx = @(v)(idst_cc(v));
    case 'nmn'
      k = 0:grid.Nx-1;
      wx = pi*k/grid.Lx;
      transformx = @(v)(dct(v));
      itransformx = @(v)(idct(v));
    otherwise
      disp(sprintfd('WTF kind of boundary condtion is %s?',grid.bcx)); 
      u=NaN;
      return;
  end

  switch grid.bcy
    case 'per'  %%% I can expliot alising if I want to 
      k = 0:grid.Ny-1;
      wy = 2*pi*k/grid.Ly;
      transformy = @(v)(fft(v.').');
      itransformy = @(v)(ifft(v.').');
    case 'dir'
      k = 1:grid.Ny;
      wy = pi*k/grid.Ly;
      transformy = @(v)(dst_cc(v.').');
      itransformy = @(v)(idst_cc(v.').');
    case 'nmn'
      k = 0:grid.Ny-1;
      wy = pi*k/grid.Ly;
      transformy = @(v)(dct(v.').');
      itransformy = @(v)(idct(v.').');
    otherwise
      disp(sprintfd('WTF kind of boundary condtion is %s?',grid.bcx)); 
      u=NaN;
      return;
  end

  % make the eigenvalues
  %
  [wx,wy]=ndgrid(wx,wy);
  Lx = 2/grid.dx^2 * (cos(wx*grid.dx)-1);
  Ly = 2/grid.dy^2 * (cos(wy*grid.dy)-1);
  L  = a-b*(Lx+Ly);

  % transform
  %
  fhat = transformx(f);
  fhat = transformy(fhat);
%  fhat = transformy(f);
%  fhat = transformx(fhat);

  
  % for sinular problems; b=0 periodic or Neumann, project onto the range
  %  
  if( L(1) == 0 )
    fhat(1) = 0;
    L(1)   = 1;
  end

  uhat = fhat./L;

  % this transform order seems to matter for the dirichelet-periodic case -- not sure why??
  %   am I somehow phase shifting data??? Something with complex part perhaps???
  %     fails with the order is dirichlet first then periodic
  %     maybe the dirichlet part is messing the the periodic part
  %  my dst looks off on the highest wave number and I'm really not sure about the inverse
  %   they are inverses of one another, so something is goofy -- need to fix these 
  %   or be simple minded and replace with matrix operations 
  %
  %  still having the same problem after adjusting the dst???
  %    I'm sure this is a problem with complex data
  %
  u = itransformx(uhat);
  u = real( itransformy(u) );
%  u = itransformy(uhat);
%  u = real( itransformx(u) );

end