function [localE] = compute_local_energy1(mys, imice, jijMat, hirMat)

localE = - jijMat(imice,:)*(mys==mys(imice)) - hirMat(imice, mys(imice));

end