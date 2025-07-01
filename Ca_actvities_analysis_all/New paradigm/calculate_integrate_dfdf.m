function [area]=calculate_integrate_dfdf(x,ind,type,ref_win)
%ind=frame.cs_start:(frame.cs_end-1);%计算面积的时间窗
%ref_win:计算时间窗之前的参考时间窗
%用不同的积分下限计算曲线下面积
dim=size(x); %按列计算
%polyfit_win=ceil(4.8/fs.ca+1):frame.cs_start-1;%取拟合时间窗

%[fs,~,frame,~,~,~,~,~,~]=setpara([]);

if nargin<3
    type='1-polyfit';
    ref_win=ceil(4.8/fs.ca+1):frame.cs_start-1;
elseif nargin<4
    ref_win=ceil(4.8/fs.ca+1):frame.cs_start-1;
end
smoothwin=3;
for ii=1:dim(2)
    xx=x(:,ii);
    if ~isempty(xx)
        %figure,plot(xx,'linewidth',2);hold on;
        lower_integrate=[];
        up_integrate=[];
        if ischar(type)
            switch type
                case 'polyfit' %sencond order polyfit为积分下限，取CS前5.2s
                    %[~,~,~,error]=test_smooth(xx,fs,frame,trial,polyfit_win,ind);
                    [p,S] = polyfit(ref_win,xx(ref_win)',4);
                    [lower_integrate,delta]=polyval(p,[ref_win ind],S);%
                    lower_integrate=lower_integrate(length(ref_win)+1:end);
                case 'lowest' %最低点的值作为积分下限
                    xxs=smoothdata(xx,'movmean',smoothwin);
                    lower_integrate=min(xxs(ind(1):ind(end)-smoothwin))*ones(size(ind));
                case 'first' %第一个点的值作为积分下限
                    xxs=smoothdata(xx,'movmean',smoothwin);
                    lower_integrate=xxs(ind(1))*ones(size(ind));
                case 'mean'
                    lower_integrate=mean(xx(ref_win))*ones(size(ind));
                case 'sum'
                    lower_integrate=0*ones(size(ind));
                case 'lowest_all' %x在ind时间窗内所有列的最低点的值作为积分下限
                    xs=smoothdata(x,'movmean',smoothwin);
                    lower_integrate=mean(min(xs(ind(1):ind(end)-smoothwin,:)'))*ones(size(ind));
                case 'event'
                    indc=ind(1:end-1);
                    xxs=smoothdata(xx([ref_win indc]),'movmean',smoothwin); %xss=sort(xxs([ref_win indc]-ref_win(1)+1));xss=xss(1:ceil(length(xss)*0.75));
                    %figure,plot(xx(indc));
                    xxxs=smoothdata(xx,'movmean',smoothwin)-xxs(length(ref_win)+1);
                    [m,indm]=max((xx(indc)-xx(indc(1))));  
                    %[~,locs]=findpeaks(abs(xx(indc)));m=xxxs(indc(locs));indm=locs;%figure,plot(xx);hold on;scatter(indc(locs),pks);improved in 20190624,because the decay phase of CR
                    dxxs= detrend(xxs(ref_win-ref_win(1)+1));
                    thres=mean(dxxs)+5*std(dxxs);
                    k=((m-xxxs([ind(1):indc(indm)]))./(indc(indm)-[ind(1):indc(indm)])');
                    %figure,plot(xxxs)
%                     %old version_just calculated rasing phase CR
%                     indc=ind(1:end-1);
%                     xxs=smoothdata(xx([ref_win indc]),'movmean',smoothwin);
%                     xxxs=smoothdata(xx,'movmean',smoothwin);
%                     [m,indm]=max(xxs(length(ref_win)+1:end));[m,indm]=max(abs(xxs(length(ref_win)+1:end)));%indc-ref_win(1)+1;indm=indc(indm);
%                     dxxs= detrend(xxs(ref_win-ref_win(1)+1));
%                     thres=mean(dxxs)+3*std(dxxs);
%                     k=(m-xxxs([ind(1):indc(indm)])./(indc(indm)-[ind(1):indc(indm)])');
%                     %k=m-xxs([ind(1):indm+ind(1)-1]-ref_win(1)+1)./(indm+ind(1)-1-[ind(1):indm+ind(1)-1])';

                    k(isinf(k))=[];
                    if ~isempty(k)
                        [sm,inds]=max(k);
                        %                 figure,plot(xx);hold on;plot(thres*ones(size(xs)));hold on;
                        %                 scatter([inds+ind(1)-1 indm+ind(1)-1],[xs(inds+ind(1)-1) xs(indm+ind(1)-1)]);hold on;
                        %                 plot([1:50],thres/ length(ind)*[1:50]+xs(indm+ind(1)-1)-(inds+ind(1)-1)*thres/ length(ind))
                        
                        if sm>= (thres-mean(xxs(length(ref_win)+1:end)))/length(ind) && ~isinf(sm)
                            xs=smoothdata(xx,'movmean',smoothwin);
                            %figure,plot(x)
                            %lower_integrate=mean(min(xs(ind(1):ind(end)-smoothwin,:)'))*ones(size(ind));
                            %lower_integrate=xs(inds+ind(1)-1)*ones(size(ind));
                            lower_integrate=(m-sm*(indm))*ones(size(ind));
                        else
                            lower_integrate=xx(ind)';
                        end
                    else
                        lower_integrate=xx(ind)';
                    end
            end
        elseif isnumeric(type)
            lower_integrate=type*ones(size(ind));
        end
        %     plot(xx,'--r','linewidth',2);hold on
        %     plot(ind,lower_integrate,'linewidth',2)
        up_integrate=xx(ind)';
        area(ii,1)=trapz(ind,[up_integrate-lower_integrate]);
        %area(ii,1)=sum(up_integrate-lower_integrate);
    else
        area(ii,1)=0;
    end
end

