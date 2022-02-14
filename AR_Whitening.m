
% 参数设定：
% 确定滤波器的阶数N，阶数应大于1
N=10;

% 设置绘制功率谱参数，参数详情可参见pwelch函数说明
Fs=10000; % 原始数据采样率
nwindow=100000; % 窗长
coe=1; % 数据乘系数，默认为1
nfft=100000;
overlap=0.5;



% 数据读取
[fileName,filePath]=uigetfile({'*.tdms; *.mat; *.txt', ...
    'Data Files (*.TS *.mat)'; '*.*','All Files (*.*)'},'Pick a file');
oriname=fileName(1:size(fileName,2)-4);
filelast=fileName((length(oriname)+1):end);

if(filelast== 'tdms')
    % Unpack data.
    [ConvertedData,ConvertVer,ChanNames, GroupNames, ChannelIndex]=convertTDMS(false,fileName);
    oriname=fileName(1:size(fileName,2)-5);
    ChannelCount=size(ChanNames{1,1},1);

    MainStruct=ConvertedData.Data.MeasuredData(1,3);
    DataLength=MainStruct.Total_Samples;
    PureDataMat=zeros(ChannelCount, DataLength);

    for i=3:1:(3+ChannelCount-1)
        MainStruct=ConvertedData.Data.MeasuredData(1,i);
        PureDataStruct{i-2}=MainStruct.Data;
        PureDataMat(i-2,:)=MainStruct.Data;
    end

else
    fpath=[filePath fileName];
    PureDataMat=load(fpath);
    if isa(PureDataMat, 'struct')
        PureDataMat=cell2mat(struct2cell(PureDataMat));
    end
end


% 转置判断
[m,n]=size(PureDataMat);
if m>n
    PureDataMat=PureDataMat';
    [m,n]=size(PureDataMat);
end



% 白化滤波
OriData=PureDataMat(1,:);

[LayerNum,DataLength]=size(OriData);

FullResultMat=zeros(LayerNum,DataLength);

for i=1:LayerNum

    ChosenData=OriData(i,:);
    % 计算互相关函数 r(x)
    [Corr,lg]=xcorr(ChosenData);
    Corr(lg<0) = [];
    % CorrDisplay=Corr(1,1:N+1);

    % 由Levinsion-Durbin算法求得反射系数 K
    [A,E,K]=levinson(Corr,N);

    % 由格型结构计算白化结果
    [WhiteningResult,BackwardError]=latcfilt(K,ChosenData);

    FullResultMat(i,:)=WhiteningResult;

end



% 绘制功率谱

% 设置功率谱绘制数据区间
nbegin=1;
nend=DataLength;

x=(OriData(:,nbegin:nend)*coe)';
y=(FullResultMat(:,nbegin:nend)*coe)';
window=hanning(nwindow);

noverlap=overlap*nfft;
[Pxx,f1]=pwelch(x,window,noverlap,nfft,Fs);
[Pyy,f2]=pwelch(y,window,noverlap,nfft,Fs);

xx=sqrt(Pxx*2);
yy=sqrt(Pyy*2);
figure(1); % 图1：原数据功率谱
loglog(f1,xx);
figure(2); % 图2：白化滤波数据功率谱
loglog(f2,yy);
figure(3); plot(FullResultMat); % 图3：白化滤波时域结果
