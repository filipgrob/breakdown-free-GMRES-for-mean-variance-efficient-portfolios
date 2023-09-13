function [xk]=bfgmres(A,b,tol,maxit)

%breakdown-free gmres function for singular systems:
%takes in input a matrix A, a vector of known terms b,
%tolerance and maximum number of iteration for stopping criteria.

%initialization
V=zeros(length(b),maxit); %matrix with a basis of Krilov space on columns 
U=zeros(length(b),maxit); %matrix that takes into account of breakdowns
H=zeros(maxit,maxit); %matrix for Arnoldi method
G=zeros(maxit,maxit); %matrix for the least square problem minimization
p=1; %number of breakdowns
k=1; %number of iterations
xold=0; 
xk=100;
hard_break=0; %flag for breakdowns: 0 if p-th breakdown already happened, 1 otherwise
beta=norm(b);
V(:,1)=b/beta;

%iterative method for solving linear system
while (norm(A*xk-b)>tol)&&(k<maxit)
    %k
    %one step of Arnoldi process
    w=A*V(:,k); %w ha numero righe=numero righe A (***)
    H(1:k,k)=V(:,1:k)'*w;
    G(1:p,k)=U(:,1:p)'*w; %G(1:p,k) al posto gk ed è vettore colonna di len p (***)
    w=w-V(:,1:k)*H(1:k,k)-U(:,1:p)*G(1:p,k); %G(1:p,k) al posto gk (***)
    H(k+1,k)=norm(w);
    H_cond=H(1:k+1,1:k); %traducendo dal paper sarebbe [Ht_old,hk';zeros(1,n_zeroes(2)),hkk];
    %check if a breakdown is happening
    if (cond(H_cond)>(10^(2*p))/tol)||(norm(xk-xold)<tol*norm(xk))
        if hard_break==1 %sono stato dentro l'if a riga 29 durante l'iterazione precedente se soddisfo questa hp (***)
            p=p+1;
            hard_break=0;
        end
        %if a breakdown occurs, matrices of previous iteration (k-1) are
        %taken into account
        U(:,p)=V(:,k); %traducendo dal paper sarebbe [U,V(:,k)];
        G(p,1:k-1)=H(k,1:k-1); %traducendo dal paper sarebbe [G;Ht_old(k,:)];
        H(k,1:k-1)=0; %H_k-1 sarebbe H(1:k,1:k-1) ma visto che voglio H_k-1(k,:) allora seleziono solo la riga k (***)
        %p
        %i=0
        while hard_break==0
            %i=i+1
            %assignin('base', 'V',V(:,1:k-1));
            %assignin('base', 'U',U(:,1:p));
            %assignin('base', 'xk',xk);
            %orthogonalize last vector that brings to a breakdown
            vt=generate_random(V(:,1:k-1),U(:,1:p),tol);
            V(:,k)=vt;
            %repeat one step of Arnoldi process until conditioning below
            %the threshold
            w=A*V(:,k);
            H(1:k,k)=V(:,1:k)'*w; %traducendo dal paper sarebbe hk=V'*w;
            G(1:p,k)=U(:,1:p)'*w; %G(1:p+1,k) al posto gk xkè nel paper sarebbe U'*w; (***)
            w=w-V(:,1:k)*H(1:k,k)-U(:,1:p)*G(1:p,k); %traducendo dal paper sarebbe w-V*hk-U*gk;
            H(k+1,k)=norm(w);
            H_cond=H(1:k+1,1:k); %traducendo dal paper sarebbe Ht_new=[Ht_old,hk';zeros(1,n_zeroes(2)),hkk];
            %cond(H_cond)
            if cond(H_cond)<1/tol
                hard_break=1;
            end
        end
        %p=p+1; %qua secondo me errore anche se nel paper lo mette qua (***)
    end
    V(:,k+1)=w/H(k+1,k);
    %construction of matrix H_tilde for the least square minimization
    if p>1
        %traducendo dal paper sarebbe G(1:p,k+1)=gk;
        %assignin('base', 'gk',gk);
        %assignin('base', 'H',H(1:k+1,1:k));
        %assignin('base', 'G',G(1:p,1:k+1));
        Htilde=[H(1:k+1,1:k);G(1:p,1:k)];
    else
        Htilde=H(1:k+1,1:k);
    end
    %solving the least square minimization
    [Q,R] = qr(Htilde);
    gm=Q*eye(length(diag(Q)),1)*beta;
    yk = R\gm;
    xold=xk;
    %update solution
    xk =V(:,1:k)*yk;
    k=k+1;
end 

end