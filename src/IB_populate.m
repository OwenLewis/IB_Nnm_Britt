function IB = IB_populate(X0)

    Nib = length(X0);
    taum = zeros(Nib,2);
    unitnormals = zeros(Nib,2);
    dsm = zeros(Nib,1);

    taum(1,:) = X0(1,:) - X0(end,:);
    dsm(1) = norm(taum(1,:));
    for i = 2:Nib
        taum(i,:) = X0(i,:) - X0(i-1,:);
        dsm(i) = norm(taum(i,:));
    end
    taup = circshift(taum,-1);
    dsp = circshift(dsm,-1);

    tau0 = (taup+taum)/2;
    dsvec = (dsm+dsp)/2;

    N = [tau0(:,2),-tau0(:,1)];

    for i = 1:Nib
        unitnormal(i,:) = N(i,:)/norm(N(i,:));
    end


    IB.Nib = Nib;
    IB.normals = unitnormal;
    IB.dsvec = dsvec;

    
end