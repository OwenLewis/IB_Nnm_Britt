function udotgradc=upwind_corner(c,U,V,grid,dt)
if size(c) ~= [grid.Nx grid.Ny]
    error('Conc. being advected does not match grid size')
end

cpad = zeros(grid.Nx+2,grid.Ny+2);

cpad(2:end-1,2:end-1) = c;

cpad(1,2:end-1) = c(end,:);
cpad(end,2:end-1) = c(1,:);
cpad(2:end-1,1) = c(:,end);
cpad(2:end-1,end) = c(:,1);
cpad(1,1) = c(end,end);
cpad(end,1) = c(1,end);
cpad(1,end) = c(end,1);
cpad(end,end) = c(1,1);


if max(U(:)) > 0
    DX = (cpad(2:end-1,:) - cpad(1:end-2,:))/grid.dx;
else
    DX = (cpad(3:end,:) - cpad(2:end-1,:))/grid.dx;
end

if max(V(:)) > 0
    DY = (cpad(:,2:end-1) - cpad(:,1:end-2))/grid.dy;
    DD = (DX(:,2:end-1) - DX(:,1:end-2))./grid.dy;
else
    DY = (cpad(:,3:end) - cpad(:,2:end-1))/grid.dy;
    DD = (DX(:,3:end) - DX(:,2:end-1))./grid.dy;
end



udotgradc = U*DX(:,2:end-1) + V*DY(2:end-1,:) - dt*U*V*DD;
end