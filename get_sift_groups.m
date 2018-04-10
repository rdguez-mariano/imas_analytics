function [ groups ] = get_sift_groups(W,vec_representants )
%SIFT_GROUPS Summary of this function goes here
%   Detailed explanation goes here

L = build_laplacian(W, 'rw');

Y = spectral_clustering(L);

classes = unique(Y);

groupmax = 0;
groupmin = Inf;
groups = {};
for i = 1:length(classes)
    gp = vec_representants(1:2,Y == classes(i));
    if (~((size(gp,1)==1)||(size(gp,2)==1)))
        groups = [groups gp];        
        groupmax = max(size(gp,2),groupmax);
        groupmin = min(size(gp,2),groupmin);
    end    
end
suptitle(['Group max = ' num2str(groupmax) ' / Group min = ' num2str(groupmin)])

end


function [Y] = spectral_clustering(L)
[U,E] = eig(L);
[eigenvalues_sorted,reorder] = sort(diag(abs(E)));

for k=1:(size(eigenvalues_sorted,1)-1)
    if(abs(eigenvalues_sorted(k+1)-eigenvalues_sorted(k))>0.1)
        break;
    end
end
chosen_eig_indices=1:k;
figure;
subplot(2,2,1);plot(chosen_eig_indices,abs(eigenvalues_sorted(chosen_eig_indices)),'*b');title('Chosen eigenvalues')
subplot(2,2,2);plot(1:size(eigenvalues_sorted,1),abs(eigenvalues_sorted),'*b');title('All eigenvalues')
subplot(2,2,3);plot(-10:10,abs(eigenvalues_sorted((k-10):(k+10))),'*b');title('Last chosen eigenvalue on the origin')
U = U(:,reorder(chosen_eig_indices));
Ui=imag(U);
subplot(2,2,4);imshow(Ui>0);title({['Imaginary part of the eigenvector matrix']; ['max = ' num2str(max(Ui(:))) '; mean = ' num2str(mean(Ui(:))) '; median = ' num2str(median(Ui(:))) '; mean_{+} = ' num2str(mean(Ui(Ui>0)))]})

% we take the real part as a saveguard... in any case, a complex matrix
% only comes from numerical errors as L is positive semi-definite
U = real(U);

Y = kmeans(U, k,'EmptyAction','singleton','Start','uniform','Replicates',3);
end

function [L] =  build_laplacian(W, laplacian_normalization)

if strcmp(laplacian_normalization, 'unn')
    D = diag(sum(W));
    L = D - W;
    L = (L+L')/2;
elseif strcmp(laplacian_normalization, 'sym')
    D = diag(sum(W));
    D_12 = sqrt(pinv(D));
    L = eye(size(W)) - D_12*W*D_12;
    L = (L+L')/2;
elseif strcmp(laplacian_normalization, 'rw')
    D = diag(sum(W));
    L = eye(size(W)) - pinv(D)*W;
else
    error('unkown normalization mode')
end
end
