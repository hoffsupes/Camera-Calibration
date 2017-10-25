function eps = calc_error(X,c,P,lam1,lam2)

% XE = P * [X; ones(1,length(X))]; %% check, X c in c,r form
% XE = XE ./ XE(3,:);
% 
% XE(1,:) = XE(1,:) - c(1,:);
% XE(2,:) = XE(2,:) - c(2,:);
% 
p1 = P(1,1:end-1);    p14 = P(1,4);
p2 = P(2,1:end-1);	p24 = P(2,4);
p3 = P(3,1:end-1);    p34 = P(3,4);

% eps = sum(sum(bsxfun(@power,XE,2))) + lam1*( norm(p3) - 1 ) + lam2*( cross(p1,p3) * cross(p2,p3)' );

eps = 0;
for i = 1:length(X)
eps = eps + ( ((p1*X(:,i) + p14)/(p3*X(:,i)+p34)) - c(1,i))^2 + ( ((p2*X(:,i) + p24)/(p3*X(:,i)+p34)) - c(2,i))^2;

end
eps = eps + lam1*( norm(p3) - 1 ) + lam2*( cross(p1,p3) * cross(p2,p3)' );


end