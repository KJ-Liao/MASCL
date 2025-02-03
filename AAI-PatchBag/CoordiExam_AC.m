function [R,T,eRMSD] = CoordiExam_AC(fromXYZ, toXYZ)
  
    len_1 = size(fromXYZ,1);
    len_2 = size(toXYZ,1);
    massc_1 = mean(fromXYZ,1); 
    massc_2 = mean(toXYZ,1);

    rfromXYZ = fromXYZ-repmat(massc_1,len_1,1);
    rtoXYZ = toXYZ-repmat(massc_2,len_2,1);
    [U,~,W] = svd(rtoXYZ'*rfromXYZ);
    R = W*[1 0 0;0 1 0;0 0 det(U*W')]*U'; 
    % [1 0 0;0 1 0;0 0 det(u*w')] to prevent more than one solution if the all point are on the same plane
  
    if (abs(norm(massc_2))<1.0e-07)
        T=-massc_1*R;
        %in the case, the ref is already centered in (0,0,0)
    else
        T=massc_2-massc_1*R;
    end

    newXYZ=fromXYZ*R + repmat(T,len_2,1);
    eRMSD=sqrt(sum(sum((toXYZ-newXYZ).^2,2))/len_2);
end