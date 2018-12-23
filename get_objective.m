function [fobj, fresidue1, fresidue2, fresidue3, freg] = get_objective(D, S, Y, W, H, A, Q, lambda, mu)

E1 = double(D)*S - double(Y);
E2 = double(A)*S - double(Q);
E3 = double(W)*S - double(H);

fresidue1 = sum(sum(E1.^2)); % reconstruction error
fresidue2 = sum(sum(E2.^2)); % optimal sparse code error
fresidue3 = sum(sum(E3.^2)); % classification error
freg = lambda*sum(sum(abs(S)));

%fobj = 0.5*(fresidue1 + mu1*fresidue2 + mu2*fresidue3) + freg;
%fobj = 0.5*(fresidue1 + mu1*fresidue2 + mu2*fresidue3);
fobj = 0.5*(mu*fresidue2 + (1-mu)*fresidue3);



