clc;clear all;
[imgname,imgpath]=uigetfile('C:\standard brain\20180901_31dpf\*.dcimg','dcimg');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\rec186952424.dcimg';
savepath=uigetdir(imgpath,'savedir');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\';
T=25;
fmt='tif';
Compression='lzw';%deflate

tic;
zpath={};
for jj=1:T
    zpath{jj,1} = checkpath([savepath '/raw tiff/z' num2str(jj, '%02d')]);
end
toc;
imio = dcimgreader([imgpath imgname], T, 1);  %, opt.t0); % auto estimate t
mon = createMonitor;
t=imio.nt;
tic;
kk=1;cpu_total=nan(imio.nframes,1);cpu_matlab=nan(imio.nframes,1);
cpu=updateMeasure(mon);
for ii=1:t
    frame=imio.readvolume(ii);
    for  jj=1:T
        imwrite(frame(:,:,jj), [zpath{jj,1} '/' num2str(ii, '%04d') '.' fmt],'Compression',Compression);
        %     frame_lzw=imread([zpath{jj,1} '/' num2str(ii, '%04d') '.' fmt]);
        %     com(ii,jj)=isequal(frame(:,:,jj), frame_lzw);   
    cpu=updateMeasure(mon);
    cpu_total(kk)=cpu.total;cpu_matlab(kk)=cpu.matlab;
    kk=kk+1;
    end
end
toc;
imio.close();
toc;
figure,
plot(cpu_total,'color','b','linewidth',2);hold on;plot(cpu_matlab,'r','linewidth',2);
legend({'cpu_total','cpu_matlab'},'interpreter','none','fontsize',15);

function cpu=updateMeasure(mon)

%// Calculate the cpu usage
cpu.total = 100 - mon.ProcPerfCounter.cpuIdle.NextValue / mon.NumOfCPU ;
cpu.matlab = mon.ProcPerfCounter.Matlab.NextValue / mon.NumOfCPU ;
end

% Calculate the cpu usage
function mon = createMonitor
MatlabProcess = System.Diagnostics.Process.GetCurrentProcess(); %// "Matlab" process
cpuIdleProcess = 'Idle' ;
mon.NumOfCPU = double(System.Environment.ProcessorCount);
mon.ProcPerfCounter.Matlab  = System.Diagnostics.PerformanceCounter('Process', '% Processor Time', MatlabProcess.ProcessName );
mon.ProcPerfCounter.cpuIdle = System.Diagnostics.PerformanceCounter('Process', '% Processor Time', cpuIdleProcess );
end



