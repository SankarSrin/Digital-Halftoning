clear, clc, close all

% signal parameters
fs = 44100;         % Hz
T = 5;              % s
N = round(fs*T);    % samples

% noise generation
xpink = pinknoise(N);       % pink noise
xred = rednoise(N);         % red (brown) noise 
xblue = bluenoise(N);       % blue noise
xviolet = violetnoise(N);   % violet noise

% sound presentation
soundsc(xpink, fs)
pause(T+1)
soundsc(xred, fs)
pause(T+1)
soundsc(xblue, fs)
pause(T+1)
soundsc(xviolet, fs)