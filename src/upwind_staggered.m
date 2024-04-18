function divflux=upwind_staggered(c,U,V,grid)
if size(c) ~= [grid.Nx grid.Ny]
    error('Conc. being advected does not match grid size')
end

if size(U) ~= [grid.Nx+1 grid.Ny]
    error('Horizontal vel. does not match grid size')
end

if size(V) ~= [grid.Nx grid.Ny+1]
    error('Vertical vel. does not match grid size')
end

xpad = zeros(grid.Nx+2,grid.Ny);
ypad = zeros(grid.Nx,grid.Ny+2);

xpad(1,:) = c(end,:);
xpad(2:end-1,:) = c;
xpad(end,:) = c(1,:);
ypad(:,1) = c(:,end);
ypad(:,2:end-1) = c;
ypad(:,end) = c(:,1);


horizflux = (xpad(1:end-1,:).*U).*(U > 0) + (xpad(2:end,:).*U).*(U < 0);

vertflux = (ypad(:,1:end-1).*V).*(V > 0) + (ypad(:,2:end).*V).*(V < 0);

divflux = diff(horizflux,1,1)./grid.dx  + diff(vertflux,1,2)./grid.dy;
end