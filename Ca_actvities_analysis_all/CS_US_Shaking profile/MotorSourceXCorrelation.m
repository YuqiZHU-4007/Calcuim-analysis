function [stimcorr,motorcorr,orthonormal_basis] = MotorSourceXCorrelation(C,reg_sens,reg_motor,lagtype)
tStart = tic;
% number of models (this code works the same way for clusters or individual cells)
nModel = size(C,1);

%% method 1:
regs = vertcat(reg_sens,reg_motor);
orthonormal_basis = regs';%Gram_Schmidt_Process(regs'); % actually is transposed?

% % C_norm = normr(C);
% betas = zeros(nModel,size(orthonormal_basis,2)+1);
% X = [ones(size(orthonormal_basis,1),1),orthonormal_basis];
% for i_model = 1:nModel
%     y = C(i_model,:)';
%     betas(i_model,:) = regress(y,X)';
% end
R=[];
for ii=1:size(C,1)   
    switch lagtype
        case 1
            [r ,~]= xcorr(C(ii,:)',orthonormal_basis,'unbiased'); % row of R: each regressor
            R(ii)=max(r);
        case 2
            [R(ii) ,~]= xcorr(C(ii,:)',orthonormal_basis,0,'unbiased'); % row of R: each regressor
    end
% figure,plot(C(ii,:)');hold on;plot(orthonormal_basis)
% figure,plot(lag,R)
end
%         [~,IX] = max(R,[],1);

numMotorRegs = size(reg_motor,1);
% stimcorr = max(betas(:,1:end-numMotorRegs),[],2);
% motorcorr = max(betas(:,end-numMotorRegs+1:end),[],2);
stimcorr = max(R(1:end-numMotorRegs,:),[],1);
motorcorr = max(R(end-numMotorRegs+1:end,:),[],1);

tElapsed = toc(tStart);
% if tElapsed>5
%     disp(['MotorSourceCorrelation tElapsed = ' num2str(tElapsed)]);
% end