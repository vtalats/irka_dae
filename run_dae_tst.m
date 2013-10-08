clear
close all;

set(0, 'defaultaxesfontsize',14,'defaultaxeslinewidth',1.0,...
       'defaultlinelinewidth',1.0,'defaultpatchlinewidth',1.0,...
       'defaulttextfontsize',18);


fprintf('Running test on the Oseen system.\n');
fprintf('This is a Stokes-type descriptor system of index-2.\n\n');


% testnum = 5;
mksiso = 0;

for testnum=2
  
%---------------------------------------
% This loads the following data:
%  system matrices - A, B, C, D, E
%  nv - number of variables in A11
%
switch testnum
  case 1  
    load oseen_test;
  case 2  
    load testdata50;
  case 3  
    load testdata;
  case 4  
    load oseen_test_73;
  case 5  
    load oseen_test_74;
  case 6  
    load stokes_test_1;
  case 7  
    load stokes_test_2;
  case 8  
    load data50_2;
  otherwise
    load oseen_test_73;
end

if mksiso==1
  B= B(:,1);
  C= C(1,:);
  D = 0;
end

% Covert matrices to a state space sytem
% sys = dss(full(A),full(B),full(C),full(D),full(E));



% Set the test run parameters
r = 30;
loglevel = 2;
tol = 1e-3;
distr_type = 1;
distr = [-1,2,0,0];

[Ar,Br,Cr,Dr,Er,siter] = irka_dae2(A,B,C,D,E,nv,r,loglevel,tol,distr_type,distr);

end
% sysred = dss(full(Ar),full(Br),full(Cr),full(Dr),full(Er));
% 
% for k = 1:size(C,1)
%   fig_title = sprintf('DAE Test Full');
%   figure('name', fig_title);
%   if exist('Hs_dae_td50.mat')
%     load Hs_dae_fullsys
%     semilogx(w,H,'k');
%     title('Bode Plot');
%     xlabel('Frequency (rad/s)');
%     ylabel('Singular Values (dB)');
%   else
%     [H,w] = bode_plt(A,B(:,k),C(k,:),D(k,k),E,-2,3,100,'k');
%     save('Hs_dae_td50','H','w');
%   end
%   
%   figure('name', 'DAE Test Reduced');  
%   bode_plt(Ar,Br(:,k),Cr(k,:),Dr(k,k),Er,-2,3,100,'r');
%   
% end
% err = norm(sysred-sys)/norm(sys);
% 
% fprintf('\nAbsolute error of the reduced system is %e\n\n',err);

