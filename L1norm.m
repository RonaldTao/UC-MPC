function [L1_norm] = L1norm(G)
%L1NORM numerically computes the L1 function norm of the stable
%transfer function G(s),

%% Determine the size of the TF matrix

l = size(G,1); m = size(G,2); % l outputs, m inputs

L1_mat  = zeros(l,m);

for i = 1:l
    for j = 1:m
        
        [y,t] = impulse(G(i,j));
        
        norm_y = zeros(1,size(t,1));
        
        for k = 1:size(t,1)
            norm_y(k) =  norm(y(k,:));
        end
        
        L1_mat(i,j) = trapz(t,norm_y);
        
    end
end


L1_vec = zeros(l,1);
for i = 1:l
   L1_vec(i,1) = sum(L1_mat(i,:)); 
end

L1_norm = max(L1_vec);

end

