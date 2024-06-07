%
% spreadmatrix_scalar.m
%
% compute the matrix for the spreading operator
% for a scalar function in the staggered grid formulation
% so cell-centered points:
%  assumes a discretization of x = xmin + dx*(i-1/2); i=1..Nx
%                              y = ymin + dy*(j-1/2); j=1..Ny
%  it also assumes a periodic domain so that xmin == xmin + dx*Nx
%                                            ymin == ymin + dx*Ny
%
% input,  X  --   matrix with ib point locations (size Nib x 2 )
%         grid -- a struct with the following fields:
%         grid.dx & grid.dy --   grid spacing
%         grid.Nx --   number of grid points in the x-direction
%         grid.Ny --   number of grid points in the y-direction
%         grid.xmin -- left edge of the domain
%         grid.ymin -- bottom edge of the domain
%
% output, S -- scaled spreading operator of size Nx*Ny x Nib
%
function S = spreadmatrix_scalar6ptBspline(X,grid)
     
  % record the number of unknowns
  %
  Nib = size(X,1);

  % record the number of grid points
  %
  Nsq = grid.Nx*grid.Ny;
 
  % convert X to grid coordinates (xg,yg)
  %
  xg = (X(:,1)-grid.xmin)/grid.dx  + 1/2;
  yg = (X(:,2)-grid.ymin)/grid.dy  + 1/2;
  
  % indices of grid point down and to the left
  %give 0 through N    
  %the 0 ones will correspond to the center points on the other side of hte
  %periodic box
  I0 = floor( xg );
  J0 = floor( yg );

   % compute shifts of the indices
  %
  Im2=I0-2;
  Im1 = I0-1;
  I1 = I0+1;
  I2 = I0+2;
  I3=I0+3;
  
  Jm2=J0-2;
  Jm1 = J0-1;
  J1 = J0+1;
  J2 = J0+2;
  J3=J0+3;
  
  % compute the wights
  %
  Wxm2= delta(Im2-xg);
  Wxm1 = delta(Im1 - xg);
  Wx0 = delta(I0 - xg);
  Wx1 = delta(I1 - xg);
  Wx2 = delta(I2 - xg);
  Wx3=  delta(I3 - xg);

  Wym2=delta(Jm2-yg);
  Wym1 = delta(Jm1 - yg);
  Wy0 = delta(J0 - yg);
  Wy1 = delta(J1 - yg);
  Wy2 = delta(J2 - yg);
  Wy3 = delta(J3 - yg);

  
  % done computing weights, make I and J's periodic
  %
  Im2=mod(Im2-1,grid.Nx) +1;
  Im1 = mod(Im1-1,grid.Nx) + 1;
  I0 = mod(I0-1,grid.Nx) + 1;
  I1 = mod(I1-1,grid.Nx) + 1;
  I2 = mod(I2-1,grid.Nx) + 1;
  I3 = mod(I3-1, grid.Nx) +1;
 
  Jm2=mod(Jm2-1,grid.Ny)+1;
  Jm1 = mod(Jm1-1,grid.Ny) + 1;
  J0 = mod(J0-1,grid.Ny) + 1;
  J1 = mod(J1-1,grid.Ny) + 1;
  J2 = mod(J2-1,grid.Ny) + 1;
  J3 = mod(J3-1,grid.Ny) + 1;

  % make four copies of each I,J corresponding to the 16 point stencil of
  % the delta function
  %
  Iv = [repmat(Im2,1,6),repmat(Im1,1,6), repmat(I0,1,6), repmat(I1,1,6), repmat(I2,1,6), repmat(I3,1,6)];
  Jv = repmat([Jm2,Jm1, J0, J1, J2, J3], 1,6);
  
  
  % compute the elements of the matrix
  %
  Wx = [repmat(Wxm2,1,6), repmat(Wxm1,1,6), repmat(Wx0,1,6), repmat(Wx1,1,6), repmat(Wx2,1,6), repmat(Wx3, 1,6)];
  Wy = repmat([Wym2, Wym1, Wy0, Wy1, Wy2, Wy3], 1,6);
  W  = Wx(:).*Wy(:);
  
  % column numbers
  %
  Kc = repmat( (1:Nib)',36,1);
  Kc = Kc(:);


  % row numbers
  %
  Kr = sub2ind([grid.Nx,grid.Ny],Iv(:),Jv(:));
  
  % make the (scaled) spreading matrix
  %
S = sparse(Kr(:),Kc(:),W(:), Nsq, Nib);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% delta -- 
%
function phi = delta(r)
%  phi = 0.25*( 1.0 + cos(0.5*pi*r));

    
ra=abs(r);
phi1=11/20 - 1/2 *ra.^2+1/4*ra.^4-1/12*ra.^5;
phi2=17/40+5/8*ra-7/4*ra.^2+5/4*ra.^3 -3/8*ra.^4 +1/24*ra.^5;
phi3=81/40-27/8*ra+9/4*ra.^2-3/4*ra.^3 +1/8*ra.^4-1/120*ra.^5; 

    phi=phi1.*double(ra<1)+phi2.*double(ra<2).*double(ra>=1)+...
        phi3.*double(ra>=2);

    phi=phi.*double(ra<3);                                                              



% ra=abs(r);
% phi1=2/3-ra.^2+1/2*ra.^3;
% phi2=4/3-2*ra+ra.^2-1/6*ra.^3; 
% 
%     phi=phi1.*double(ra<1)+phi2.*double(ra<2).*double(ra>=1);
% 
%     phi=phi.*double(ra<2);                                                              
% 
% 
% 







  
  
  
  
