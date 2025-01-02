clc
clear all
close all

load(['ExampleData.mat'])
Granularity=8;
ChunkSize=1000;
[IDX]=ACLDP(data,Granularity,ChunkSize);
figure
for i=unique(IDX)'
    plot(data(IDX==i,1),data(IDX==i,2),'.','markersize',8)
    hold on
end
hold off
evaluation = evalclusters(data, IDX,'CalinskiHarabasz')
evaluation = evalclusters(data, IDX,'DaviesBouldin')