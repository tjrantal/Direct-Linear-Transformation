close all;
clear all;
clc;

addpath('functions');


imageNames = {'sampleFigs/IMG_20201029_093356.JPG','sampleFigs/IMG_20201029_093415.JPG'};
for i = 1:length(imageNames)
    cam(i).image = imread(imageNames{i});
end

%Distort images
for i = 1:length(cam)

     %Distort image
     distorted = distort(cam(i).image,[0.1 0.01 0.001 0 0]);
end

