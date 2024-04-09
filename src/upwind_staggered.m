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


rightflux = (c.*U(2:end,:)).*(U(2:end,:) < 0);
leftflux = (c.*U(1:end-2,:)).*(U(1:end-1,:) > 0);

topflux = (c.*V(:,2:end)).*(V(:,2:end) < 0);
botflux = (c.*V(:,1:end-1)).*(V(:,1:end-1) > 0);

divflux = (rightflux - leftflux)./grid.dx  + (topflux - botflux)./grid.dy;
end