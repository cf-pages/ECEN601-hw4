H = hilb(15);
[q1,r1]=gramschmidt(H);
[q2,r2]=modified_gramschmidt(H);
[q3,r3]=qr(H,0);
max(max(abs(H-q1*r1)))
max(max(abs(H-q2*r2)))
max(max(abs(H-q3*r3)))

I = eye(15);
max(max(abs(q1*q1'-I)))
max(max(abs(q2*q2'-I)))
max(max(abs(q3*q3'-I)))


function [Q,R]=gramschmidt(A)
[~,n] = size(A);
V = zeros(size(A));
Q = zeros(size(A));
R = zeros(n);
for j = 1:n
    V(:,j) = A(:,j);
    for i = 1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        V(:,j) = V(:,j)-R(i,j)*Q(:,i);
    end
    R(j,j) = norm(V(:,j),2);
    Q(:,j) = V(:,j)/R(j,j);
end
end

function [Q,R]=modified_gramschmidt(A)
[~,n] = size(A);
Q = zeros(size(A));
R = zeros(n);
V=A;
for i = 1:n
    R(i,i) = norm(V(:,i),2);
    Q(:,i) = V(:,i)/R(i,i);
    for j = i+1:n
        R(i,j) = Q(:,i)'*V(:,j);
        V(:,j) = V(:,j)-R(i,j)*Q(:,i);
    end
end
end

