function [ Ar, Br, Cr, Dr, Er, siter ] = irka_dae2( A, B, C, D, E, nv, r, loglevel,tol,distr_type,distr )
%IRKA_MIMO Reduces an nxn MIMO dynamical system to an rxr MIMO reduced system using IRKA.
%
%   Usage:
%     
%   Inputs:
%     A, B, C, D, E - state space matrices where
%
%     |E11  0||x1_dot|   |A11 A12||x1|   |B1|
%     |      ||      | = |       ||  | + |  |u(t)
%     | 0   0||x2_dot|   |A21  0 ||x2|   |B2|
%     and
%     y(t) = C1*x1 + C2*x2 + D*u(t)
%
%     n1 - dimension of A11
%     r - size of the reduced order model.  Typically, this should be
%         r << n, and most certainly r < n.  If r >= n, then no reduction 
%         is performed and the original system is returned with no error.
%     loglevel
%     tol - tolerence for stopping criteria
%     distr_type - (1-5) type of distribution for the initial shifts
%     distr - [
%
%
%   Outputs:
%     Ar, Br, Cr, Dr - reduced system state space matrices 
% 
%     
%
% global S0

if nargin < 10
  distr = [-2,2,-2,1];
  if nargin < 9
    distr_type = 2;
    if nargin < 8
      tol = 1e-3;
      if nargin < 7
        loglevel = 2;
      end
    end
  end
end

max_iter = 50;

calc_Hinf_err = 1;

n = size(A,1);


if size(A,2) ~= n
  error('The ''A'' matrix must be square.');
elseif (size(D,1) ~= size(C,1)) || (size(D,2) ~= size(B,2))
  error('The ''D'' matrix is the wrong size.  It should be (num outouts) x (num inputs).');
end

logName = [datestr(now,'mmddyyyy') '.irkad2'];
lf = Msgcl(loglevel,logName);

lf.pmsg(lf.ALL,'==========================================================');
lf.pmsg(lf.ALL,'* Beginning IRKA (DAE2) model reduction');
lf.pmsg(lf.ALL,'*      n = %d ==> r = %d',n,r);
lf.pmsg(lf.ALL,'*   Function: irka_mimo.m');
lf.pmsg(lf.ALL,'*   Log level is set to %d.',loglevel);
lf.pmsg(lf.ALL,'*   sigma tolerence = %e',tol);


if nv < r
  lf.pmsg(lf.WARN, 'Since r > n1, the original system is being returned.');
  Ar = A;
  Br = B;
  Cr = C;
  Dr = D;
  Er = E;
else
  tic;
  % Getting all the sub-matrices
  lf.pmsg(lf.PED, 'Extracting A11,A12,A21,E11,B1,C1,C2');
  [E11, A11, A12, A21, B1, B2, C1, C2] = split_dae_mats(A, B, C, E, nv);

  
  % Generate matrices to be used for reduction
  [Bk,Ck,Dkt,Dkpt] = getBCD(A11, A12, B1, B2, C1, C2, D, E11);
  
  % Setting the initial shifts
  lf.pmsg(lf.PED, 'Creating %d sigma values using logspace.',r);
  if distr_type
    sigma_old = set_init_shifts(distr,r,distr_type,lf);
  else
    % Set the initial shifts by taking the first 10 conjugate pairs
    sigma_old = -1./eigs(E,A,30);
    [~,ord] = sort(abs(imag(sigma_old)),'descend');
    sigma_old = sigma_old(ord);
    sigma_old = sort(sigma_old(1:20))
    figure('name','Initial Shifts','NumberTitle','off');
    plot(sigma_old,'*b');
    title('Initial Shifts - DAE');
    xlabel('Real');
    ylabel('Imaginary');
  end
  %siter(:,1) = sigma_old(:);
  siter = [];
  
  lf.pmsg(lf.PED, 'Creating %d left and right directional vectors.',r);
  
  [b,c] = set_init_direct(A11,Bk,Ck,D,E11,sigma_old,lf);
  
  lf.pmsg(lf.WARN, 'Build the initial hermite sub space');
  [Vr,Wr] = getHermSSp(sigma_old,A11,A12,A21,Bk,b,Ck,c,E11,lf);
  
  lf.pmsg(lf.PED,'Creating initial reduced matrices.');
  Ar = Wr'*A11*Vr;
  Er = Wr'*E11*Vr;
  Br = Wr'*Bk;
  Cr = Ck*Vr;
  Dr = D;
  
  if calc_Hinf_err
    sysr = dss(Ar,Br,Cr,full(Dr),Er);
  end
  
  
  k = 1;
  error_in_s = 1;
  min_err = 1;
  
  while (k <= max_iter) && (error_in_s > tol)
    lf.pmsg(lf.ERR,'Starting irka iteration %d',k);
   
    lf.pmsg(lf.PED,'  Eigen-system computation.')
    [X,S]=eig(full(Ar),full(Er));

    lf.pmsg(lf.PED,'  Build new shifts and directions.');
    sigma_new = -diag(S);   
    neg_ent = find(sigma_new<0);
    sigma_new(neg_ent) = -sigma_new(neg_ent);
        
    EX = X\(Er*X);
    b = (EX\(X\Br))';
    c = Cr*X;   
    
    % sort the entries 
    [sigma_new,ord] = sort(sigma_new);
    c = c(:,ord);
    b = b(:,ord);
    
    siter(:,k) = sigma_new(:);
    
    lf.pmsg(lf.WARN, '  Building the new Hermite subspace.')
    [Vr,Wr] = getHermSSp(sigma_new,A11,A12,A21,Bk,b,Ck,c,E11,lf);
    
    lf.pmsg(lf.PED,'  Creating reduced matrices for this iteration');
    Ar = Wr'*A11*Vr;
    Er = Wr'*E11*Vr;
    Br = Wr'*Bk;
    Cr = Ck*Vr;
    
    if calc_Hinf_err
      lf.pmsg(lf.PED, 'Calculating the relative error of the transfer function.');
      sysr_new = dss(Ar,Br,Cr,full(D),Er);
      error_in_H = norm(sysr-sysr_new,inf)/norm(sysr,inf);
      sysr = sysr_new;
    end

    lf.pmsg(lf.PED, '  Calculating the relative error of the sigma shift vector.');
    error_in_s = norm(sigma_new-sigma_old)/norm(sigma_old);
    sigma_old = sigma_new;

    if error_in_s < min_err
      min_err = error_in_s;
    end
    
    lf.pmsg(lf.ERR, '  Sigma error = %e',error_in_s);
    if calc_Hinf_err
      lf.pmsg(lf.ERR, '      H error = %e',error_in_H);
    end;
    k = k+1;
  end


  tf = toc;
  if k <= max_iter
    lf.pmsg(lf.ALL,'* Converged.');
    lf.pmsg(lf.ALL,'* Total iterations = %d',k-1);
    lf.pmsg(lf.ALL,'* Final sigma error = %e',error_in_s);
  else
    lf.pmsg(lf.ALL,'* Failed to converge.');
    lf.pmsg(lf.ALL,'* Minimum sigma error = %e',min_err);
  end
  lf.pmsg(lf.ALL,'* Total elapsed time = %f',tf);
  lf.pmsg(lf.ALL,'* End IRKA model reduction');
  lf.pmsg(lf.ALL,'==========================================================');
end

end

