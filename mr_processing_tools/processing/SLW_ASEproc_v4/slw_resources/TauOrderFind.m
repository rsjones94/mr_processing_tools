%%%%%%% SLW used this to find the order of tauvec that is most trendless.

%%% Need to find a tau order that does not have a trend, so we can 
%%% regress out signal drift
tautime=0:0.5:20;
tauveclen=length(tautime);
slopestore=1;
interceptstore=1;
tauvecstore=tautime;
iterations=0;
tic
while iterations<5000000
    tauvec=tautime(randperm(tauveclen));
    %%% p(1) = slope; p(2) = intercept
    p=polyfit(tautime,tauvec,1);
    if abs(p(1))<slopestore & abs(p(2)-max(tautime)/2)<interceptstore
        slopestore=p(1);
        interceptstore=p(2);
        tauvecstore=tauvec;
        disp('tau update')
        toc
    end
    iterations=iterations+1;
    if mod(iterations,100000)==0
        disp(['iterations is ',num2str(iterations)])
    end
end
disp('doneproc')

%% These results are better than a diamond pattern. Compare to indGen.m

%%%% after 4759051 random permutations the following tauvec was identified
%%%% as being the most trendless, according to a linear fit.
%%%%%%% For tau = 0:0.5:20
% % % tauvec=[19 20 6.5 8.5 1 16.5 1.5 14.5 3.5 13 9 19.5 17 2.5 ...
% % %   4.5 12.5 2 9.5 6 14 7.5 5.5 12 8 11 0.5 16 10.5 13.5 5 ...
% % %   10 0 18.5 17.5 3 11.5 7 18 15.5 15 4];
%%% 
%%% p = -0.0009   10.0087















% % % For tau=0:0.5:22.5 after 1000000 iterations
tauvec=[20.5 11 1 20 4.5 18.5 17 8 9.5 12 21.5 17.5 8.5 22 14.5 19 7 7.5 9 ...
    6.5 2 15 18 5.5 16 10 22.5 10.5 0.5 6 13 4 14 21 0 5 2.5 1.5 11.5 ...
    12.5 16.5 13.5 19.5 3 15.5 3.5];





                      
