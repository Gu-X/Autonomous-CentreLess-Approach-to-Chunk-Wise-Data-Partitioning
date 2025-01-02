function [IDX]=ACLDP(data,G,BS)
[L,~]=size(data);
seq0=[1:BS:L,L+1];
ii=1;
[NormC,SumC,S,NC]=clclustering(data(seq0(ii):1:seq0(ii+1)-1,:),G);
[NormC,SumC,S,NC]=clmerge(NormC,SumC,S,NC);
for ii=2:1:length(seq0)-1
    [NormC1,SumC1,S1,NC1]=clclustering(data(seq0(ii):1:seq0(ii+1)-1,:),G);
    [NormC,SumC,S,NC]=clmerge([NormC;NormC1],[SumC;SumC1],[S;S1],NC+NC1);
end
dist1=zeros(L,NC);
X=sum(data.^2,2);
for jj=1:1:NC
    dist1(:,jj)=(X*S(jj)+NormC(jj)-2*data*SumC(jj,:)')./S(jj);
end
[~,IDX]=min(dist1,[],2);
end
function [NormC,SumC,S,NC]=clmerge(NormC,SumC,S,NC)
dist1=zeros(NC);
for ii=1:1:NC
    for jj=ii:1:NC
        dist1(ii,jj)=clusterdist(SumC(ii,:),SumC(jj,:),NormC(ii),NormC(jj),S(ii),S(jj));
    end
end
dist1=dist1+dist1';
MGIDX=1;
while(MGIDX~=0)
    MGIDX=0;
    Q=[];
    P=[];
    seq1=diag(dist1);
    dist2=repmat(seq1,1,NC)+repmat(seq1',NC,1);
    dist2=(dist2-dist1);
    dist2=triu(dist2,1);
    [v,~]=max(dist2,[],'all');
    if v>0
        MGIDX=1;
        [Q,P]=find(dist2==v);
        Q=Q(1);
        P=P(1);
    end
    if MGIDX==1
        S1=S(Q)+S(P);
        S([Q,P])=[];
        S(end+1,1)=S1;
        NormC1=NormC(Q)+NormC(P);
        NormC([Q,P])=[];
        NormC(end+1,1)=NormC1;
        SumC1=SumC(Q,:)+SumC(P,:);
        SumC([Q,P],:)=[];
        SumC(end+1,:)=SumC1;
        dist2=zeros(NC-1);
        dist1([Q,P],:)=[];
        dist1(:,[Q,P])=[];
        dist2(1:1:NC-2,1:1:NC-2)=dist1;
        dist1=dist2;
        jj=NC-1;
        for ii=1:1:NC-1
            dist1(ii,jj)=(NormC(ii)*S(jj)+NormC(jj)*S(ii)-2*SumC(ii,:)*SumC(jj,:)')./(S(ii)*S(jj));
        end
        dist1(jj,:)=dist1(:,jj)';
        NC=NC-1;
    end
end
end
function [NormC,SumC,S,NC]=clclustering(data,G)
[L,W]=size(data);
dist0=pdist(data,'euclidean').^2;
dist1=squareform(dist0);
averdist=mean(dist0);
for ii=1:1:G
    dist0(dist0>averdist)=[];
    averdist=mean(dist0);
end
delta=averdist;
dist10=exp(-1*dist1./(2*delta));
tempseq=sum(dist10,2);
SPM=zeros(L);
SPM(dist1<=averdist)=1;
NN=ceil(mean(sum(SPM,2)));
[~,SPM2]=sort(dist1+eye(L)*(max(dist1,[],'all')+1),2,'ascend');
seq1=max(tempseq(SPM2(:,1:1:NN)),[],2);
seq=find((tempseq-seq1)>0);
NC=length(seq);
Weights1=zeros(NC,L);
Weights2=zeros(NC,L);
%%
for ii=1:1:NC
    Weights1(ii,:)=sum(dist1([seq(ii),SPM2(seq(ii),1:1:NN)],:),1);
    Weights2(ii,:)=min(dist1([seq(ii),SPM2(seq(ii),1:1:NN)],:),[],1);
end
[M1,IDX]=sort(Weights2,1,'ascend');
IDX=IDX(1,:);
seq0=M1(1,:)-M1(2,:)==0;
[~,IDX(seq0)]=min(Weights1(:,seq0),[],1);
%%
NormC=zeros(NC,1);
SumC=zeros(NC,W);
S=zeros(NC,1);
for ii=1:1:NC
    seq2=find(IDX==ii);
    data1=data(seq2,:);
    S(ii)=length(seq2);
    NormC(ii)=sum(data1.^2,'all');
    SumC(ii,:)=sum(data1,1);
end
end
function [dist1]=clusterdist(data1,data2,norm1,norm2,s1,s2)
dist1=(norm1*s2+norm2*s1-2*data1*data2')./(s1*s2);
end