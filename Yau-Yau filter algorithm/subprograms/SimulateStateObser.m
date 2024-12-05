% function [state, obser, s] = SimulateStateObser(T, Tinc, f, h, Dim)
% nPaths = 10;
% nT     = T/Tinc + 1;
% state  = zeros(nT, Dim);
% obser  = zeros(nT, Dim);
% % state_mean = zeros(nT, Dim);
% 
% % state_samples = zeros(nT, Dim, nPaths);
% %state_samples = reshape(state_samples,nT,Dim*nPaths);
% sqrtdT = sqrt(Tinc);
% s = rng('default');
% 
% % x0 ~ N(0,1), E[x0] = 0: the initial state
% % for p = 1:nPaths
% %     for t = 2:nT
% %         state_samples(t,:,p) = state_samples(t-1,:,p) + ( f( state_samples(t-1,:,p) ) * Tinc ) + ( sqrt(1)*sqrtdT * randn(1,Dim) );
% %     end
% % end
% for t = 2:nT
%     %state(t,:) = sum(state_samples, 3) + ( f( state_mean(t-1,:) ) * Tinc );
%     %state(t,:) = state_samples(t,:,1) ;
%     %state(t,:) = sum(state_samples(t,:,:), 3)/10 ;
%     state(t,:) = state(t-1,:) + ( f( state(t-1,:) ) * Tinc ) + ( sqrt(0.05)*sqrtdT * randn(1,Dim) );
% end
% for t = 2:nT
%     obht_1= h( state(t-1,:));
%     %a = py.torch.tensor(state(t-1,:));
%     %obht_1 = h( py.torch.unsqueeze(a,dim=py.int(0))).detach().numpy();
%     %obht_1 = h( py.torch.unsqueeze(a,dim=py.int(0)));
%     %obhmean = h( state_mean(t-1,:));
%     %obser(t,1) = obser(t-1,1) + ( obht_1(1) * Tinc ) + ( sqrt(0.5)*sqrtdT * randn(1,1) );
%     %obser(t,2:Dim) = obser(t-1,2:Dim) + ( obhmean(2:Dim)  * Tinc ) + ( sqrt(0.15)*sqrtdT * randn(1,Dim-1) );
%     %disp(state(t-1,:)*B);
%     %obser(t,:) = reshape(state_samples(t-1,:,:),1,Dim*nPaths);
%     obser(t,:) = obser(t-1,:) + ( obht_1  * Tinc ) + (sqrt(0.15)*sqrtdT * randn(1,3) );
% end
% 
function [state, obser, s] = SimulateStateObser(T, Tinc, f, h, Dim, y1)
nT     = T/Tinc + 1;
state  = zeros(nT, Dim);
obser  = zeros(nT, Dim);
sqrtdT = sqrt(Tinc);
s = rng("default");
% x0 ~ N(0,1), E[x0] = 0: the initial state
%state(1,:) = [26.7369,16.6951,14.6727];
%state(1,:) = [11.3584,9.42894,9.27732];
for t = 2:nT
    state(t,:) = state(t-1,:) + ( f( state(t-1,:) ) * Tinc ) + ( sqrtdT * randn(1,Dim) );
    %state(t,:) = state(t-1,:) + ( f( state(t-1,:) ) * Tinc ) ;
end
% A = zeros(3,3);
%  A(1,:) = [4.588,-13.602,-10.783];
%  A(2,:) = [1.840,-5.89,-5.366];
%  A(3,:) = [1.840,-4.792,-4.474];
%  B = orth(A);
%  TAB = B/A;
% obser(1,:) = (TAB*[11.3584,16.6951,14.6727]')';
 %obser(1,:) = [11.3584,9.42894,9.27732];
for t = 2:nT
    obser(t,:) = obser(t-1,:) + ( h( state(t-1,:) ) * Tinc )+( sqrt(1)*sqrtdT * randn(1,Dim) );
    %obser(t,:) = obser(t-1,:) + ( h( y1(t-1,:) ) * Tinc ) + ( sqrt(0.5)*sqrtdT * randn(1,Dim) );
    %obser(t,:) = obser(t-1,:) + ( h( obser(t-1,:) ) * Tinc )+( sqrt(1)*sqrtdT * randn(1,Dim) );
    %obser(t,:) = y1(t-1,:) +( sqrtdT * randn(1,Dim) );
%state = obser;
end