function [ groups, Y, newW ] = get_sift_groups(W,vec_representants,laplacian )
%SIFT_GROUPS Summary of this function goes here
%   Detailed explanation goes here
tic;
L = build_laplacian(W, laplacian);
disp(['The Laplacian was built in ' num2str(toc) ' seconds']);

Y = spectral_clustering(L,laplacian);

classes = unique(Y);
newW = zeros(size(W));

groupmax = 0;
groupmin = Inf;
groups = {};
pos = 1;
for i = 1:length(classes)
    gp = vec_representants(1:2,Y == classes(i));
    %if (~((size(gp,1)==1)||(size(gp,2)==1)))
        groups = [groups gp];
        groupmax = max(size(gp,2),groupmax);
        groupmin = min(size(gp,2),groupmin);
        finalpos = sum(Y == classes(i))+pos-1;
        newW(pos:finalpos,pos:finalpos) = W(Y == classes(i),Y == classes(i));
        pos = finalpos;
    %end
end
suptitle(['Group max = ' num2str(groupmax) ' / Group min = ' num2str(groupmin)])
figure;
subplot(1,2,1);imshow(newW>0);title('New arrangement')
subplot(1,2,2);imshow(W>0);title('Old arrangement')
suptitle('Similarity matrix');

end


function [Y] = spectral_clustering(L,laplacian)
tic;
[U,E] = eig(L);
[eigenvalues_sorted,reorder] = sort(diag(abs(E)));


for k1=1:(size(eigenvalues_sorted,1)-1)
    if(abs(eigenvalues_sorted(k1+1)-eigenvalues_sorted(k1))>0.1)
        break;
    end
end

eps = 0.1;
for k2=1:(size(eigenvalues_sorted,1)-1)
    if (((eigenvalues_sorted(k2+1)>0.02))||((eigenvalues_sorted(k2+1)>10^(-3))&&(abs((eigenvalues_sorted(k2+1)+eps)/(eigenvalues_sorted(k2)+eps))>1.1)))
        break;
    end
end
disp(['Choosing zero eigenvalues in ' num2str(toc) ' seconds']);

k = min(k1,k2);

chosen_eig_indices=1:k;
figure;
subplot(2,2,1);plot(chosen_eig_indices,abs(eigenvalues_sorted(chosen_eig_indices)),'*b');title('Chosen eigenvalues')
subplot(2,2,2);plot(1:size(eigenvalues_sorted,1),abs(eigenvalues_sorted),'*b');title('All eigenvalues')
eshowmin = max(k-10,1);
eshowmax = min(k+10,size(eigenvalues_sorted,1));
subplot(2,2,3);plot((eshowmin:eshowmax) - k,abs(eigenvalues_sorted(eshowmin:eshowmax)),'*b');title('Last chosen eigenvalue on the origin')
U = U(:,reorder(chosen_eig_indices));
Ui=imag(U);
subplot(2,2,4);imshow(Ui>0);title({['Imaginary part of the eigenvector matrix']; ['max = ' num2str(max(Ui(:))) '; mean = ' num2str(mean(Ui(:))) '; median = ' num2str(median(Ui(:))) '; mean_{+} = ' num2str(mean(Ui(Ui>0)))]})

% we take the real part as a saveguard... in any case, a complex matrix
% only comes from numerical errors as L is positive semi-definite
tic;
U = real(U);
if strcmp(laplacian, 'sym')
   U = U./repmat(sqrt(sum(U.^2,2)),1,size(U,2));
end
if ( max(abs(sqrt(sum(U.^2,2))-1))<=0.000001 )
    disp('Eigenvector matrix with normalised rows !')
else
    disp('Eigenvector matrix with non-normalised rows !')
end
    
disp(['Computing the matrix of associated eigenvectors in ' num2str(toc) ' seconds']);

tic;
if strcmp(laplacian, 'rw')
    Y = kmeans(U, k,'EmptyAction','singleton','Start','uniform','Replicates',3);
else
    centers = eye(k);
    Y = kmeans(U, k,'EmptyAction','singleton','Start',centers,'Replicates',1);
end
disp(['Kmeans in ' num2str(toc) ' seconds']);
end

function [L] =  build_laplacian(W, laplacian_normalization)

if strcmp(laplacian_normalization, 'unn')
    D = diag(sum(W));
    L = D - W;
    L = (L+L')/2;
elseif strcmp(laplacian_normalization, 'sym')
    D = sum(W);    
    D = sqrt(D);
    D(D~=0) = 1./D(D~=0);
    D = diag(D);
    L = eye(size(W)) - D*W*D;
    L = (L+L')/2;
elseif strcmp(laplacian_normalization, 'rw')
    D = sum(W);        
    D(D~=0) = 1./D(D~=0);
    D = diag(D);

    L = eye(size(W)) - D*W;
else
    error('unkown normalization mode')
end
end
