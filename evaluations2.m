function [j, d, fmax, fmin] = evaluations2(T, R, Phi1,Phi2)


%     imRegionGray = rgb2gray(T);
    thresh = graythresh(T);
    BW_T = imbinarize(T,thresh);
    
%     [m,n]=size(T);
% 	locx=(m-1)/2;
%     locy=(n-1)/2;
%     mask = false(m,n);
%     mask(locx,locy) = true;
%     
%     W1 = graydiffweight(T, mask, 'GrayDifferenceCutoff', 25);
%     figure(1)
%     imshow(log(W1),[])
%     [BW_T, D1] = imsegfmm(W1, mask, thresh);
    
%     figure(2)
%     imshow(BW_T)
%     figure(3)
%     imshow(D1)

%     W2 = graydiffweight(R, mask, 'GrayDifferenceCutoff', 25);
%     figure(4)
%     imshow(log(W2),[])
%     [BW_R, D2] = imsegfmm(W2, mask, thresh);
    
%     imRegionGray = rgb2gray(R);
    thresh = graythresh(R);
    BW_R = imbinarize(R,thresh);
%     figure(5)
%     imshow(BW_R)
%     figure(6)
%     imshow(D2)
    
    
    j = jaccard(BW_T,BW_R);
    d = dice(BW_T,BW_R);
    
    
    
    
    [f0, ~] = compute_JD_and_Curl(Phi1,Phi2,1);
    fmax=max(max(f0));
    fmin=min(min(f0));
end

