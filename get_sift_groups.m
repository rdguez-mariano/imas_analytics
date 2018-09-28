function [ Y, newW, newVec, newArr ] = get_sift_groups(W,vec_representants,laplacian )
%SIFT_GROUPS Summary of this function goes here
%   Detailed explanation goes here
visual = true;
verbose = false;

tic;
L = build_laplacian(W, laplacian);
if(verbose)
    disp(['The Laplacian was built in ' num2str(toc) ' seconds']);
end
Y = spectral_clustering(L,laplacian,visual,verbose);

classes = unique(Y);
newW = zeros(size(W));

groupmax = 0;
groupmin = Inf;

% New arrangement positions for the W matrix -> will be output in newW
pos = 1;
newArr = {};
for i = 1:length(classes)
    finalpos = sum(Y == classes(i))+pos-1;
    newArr = [newArr pos:finalpos];
    pos = finalpos+1;
end
newVec = -1*ones(size(vec_representants));
% Create newW
for i = 1:length(classes)
    newVec(:,newArr{i}) = vec_representants(:,Y == classes(i));
    
    groupmax = max(size(newArr{i},2),groupmax);
    groupmin = min(size(newArr{i},2),groupmin);
    newW(newArr{i},newArr{i}) = W(Y == classes(i),Y == classes(i));
    
    for j=(i+1):length(classes)
        newW(newArr{i},newArr{j}) = W(Y == classes(i),Y == classes(j));
        newW(newArr{j},newArr{i}) = W(Y == classes(j),Y == classes(i));
    end
end

if(visual)
    suptitle(['Group max = ' num2str(groupmax) ' / Group min = ' num2str(groupmin)])
    figure;
    subplot(1,2,1);imshow(newW>0);title('New arrangement')
    subplot(1,2,2);imshow(W>0);title('Old arrangement')
    suptitle('Similarity matrix');
end

end


function [Y] = spectral_clustering(L,laplacian,visual,verbose)
tic;
[U,E] = eig(L);
[eigenvalues_sorted,reorder] = sort(diag(abs(E)));

eps1 = 0.1;
theta = pi/35;
b = [];
for k1=1:(size(eigenvalues_sorted,1)-1)
    if(eigenvalues_sorted(k1+1)>(k1+1)*tan(theta))
        break;
    end
    
    % forecast eig(k1+1) with last regressed line
    if (~isempty(b))
        y = [1 (k1+1)]*b;
        if (y+eps1<eigenvalues_sorted(k1+1))
            break;
        end
    end
    
    % line regression for x = 1,...,k1+1 and y = eig(1),...,eig(k1+1)
    xvec = (1:(k1+1))';
    lastxvec = (1:k1)';
    Coeff = [ones(size(xvec)) xvec];
    b = Coeff\eigenvalues_sorted(xvec);
    
    % max vertical distance between the line and real eigenvalues up to k1
    % (not k1+1)
    if (~isempty(lastxvec))
        y = [ones(length(lastxvec),1) lastxvec]*b;
        diffmax = max(abs( y - eigenvalues_sorted(lastxvec) ));
        % do not accept eig(k1+1) if the regressed line lies too far
        % from a point to be noise
        if (diffmax>eps1)
            break;
        end
    end
    %     figure;hold on;
    %     plot(xvec,eigenvalues_sorted(xvec),'*b');hold on;
    %     plot(0:0.1:4, [ones(size(0:0.1:4,2),1) (0:0.1:4)']*b);
    %     plot(xvec,y,'ob');
end
k = k1;

% for k1=1:(size(eigenvalues_sorted,1)-1)
%     if(abs(eigenvalues_sorted(k1+1)-eigenvalues_sorted(k1))>0.15)
%         break;
%     end
% end
%
% eps = 0.1;
% for k2=1:(size(eigenvalues_sorted,1)-1)
%     if (((eigenvalues_sorted(k2)>eps)&&(abs((eigenvalues_sorted(k2+1)+eps)/(eigenvalues_sorted(k2)+eps))>1.1)))
%         break;
%     end
% end
% k = min(k1,k2);

if(verbose)
    disp(['Choosing zero eigenvalues in ' num2str(toc) ' seconds']);
end


chosen_eig_indices=1:k;
if (visual)
    figure;
    subplot(2,2,1);plot(chosen_eig_indices,abs(eigenvalues_sorted(chosen_eig_indices)),'*b');title('Chosen eigenvalues')
    %subplot(2,2,2);plot(1:size(eigenvalues_sorted,1),abs(eigenvalues_sorted),'*b');title('All eigenvalues')
    subplot(2,2,2);plot(diff(diff(abs(eigenvalues_sorted(chosen_eig_indices)))),'*b');title('2nd derivate')
    eshowmin = max(k-10,1);
    eshowmax = min(k+10,size(eigenvalues_sorted,1));
    subplot(2,2,3);plot((eshowmin:eshowmax) - k,abs(eigenvalues_sorted(eshowmin:eshowmax)),'*b');title('Last chosen eigenvalue on the origin')
    subplot(2,2,4);plot(diff(abs(eigenvalues_sorted(chosen_eig_indices))),'*b');title('derivate')
    
end
U = U(:,reorder(chosen_eig_indices));
if (visual)
    Ui=imag(U);
    %subplot(2,2,4);imshow(Ui>0);title({['Imaginary part of the eigenvector matrix']; ['max = ' num2str(max(Ui(:))) '; mean = ' num2str(mean(Ui(:))) '; median = ' num2str(median(Ui(:))) '; mean_{+} = ' num2str(mean(Ui(Ui>0)))]})
end

% we take the real part as a saveguard... in any case, a complex matrix
% only comes from numerical errors as L is positive semi-definite
tic;
U = real(U);
if strcmp(laplacian, 'sym')
    U = U./repmat(sqrt(sum(U.^2,2)),1,size(U,2));
end

if(verbose)
    if ( max(abs(sqrt(sum(U.^2,2))-1))<=0.000001 )
        disp('Eigenvector matrix with normalised rows !')
    else
        disp('Eigenvector matrix with non-normalised rows !')
    end
    disp(['Computing the matrix of associated eigenvectors in ' num2str(toc) ' seconds']);
end

tic;
if strcmp(laplacian, 'rw')
    Y = kmeans(U, k,'EmptyAction','singleton','Start','uniform','Replicates',3);
else
    centers = eye(k);
    Y = kmeans(U, k,'EmptyAction','singleton','Start',centers,'Replicates',1);
end
if (verbose)
    disp(['Kmeans in ' num2str(toc) ' seconds']);
end
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
