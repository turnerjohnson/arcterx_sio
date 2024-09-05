function dudx=secondderivative(u,dim,dx)
% function dudx=secondderivative(u,dim,dx) calculates the second derivative
% of the variable u along dimension dim where dx is the spacing between
% values of u.
% dudx is the same size as u but the edges of u are returned with nans
%
% D. Rudnick, 10 April 2020

nn=size(u);

newsiz=[prod(nn(1:dim-1)) nn(dim) prod(nn(dim+1:end))];
uu=reshape(u,newsiz);
duudx=NaN(size(uu));
duudx(:,2:nn(dim)-1,:)=(uu(:,1:nn(dim)-2,:)-2*uu(:,2:nn(dim)-1,:)+uu(:,3:nn(dim),:))/(dx*dx);

dudx=reshape(duudx,nn);
