%% DataMatrix Creation

disp('Data Matrix Is Being Created')

foldnumb = 7; %number of the video that we are interested in
m = readtable('batch_8004.txt'); %name of the batch of detections.

A = table2array(m);
deneme = A(:,1);
targetvalue = foldnumb*100000;
nontargetvalue = targetvalue + 100000;
firstloc = find(((deneme>targetvalue)&(deneme<nontargetvalue)),1,'first');
lastloc = find(((deneme>targetvalue)&(deneme<nontargetvalue)),1,'last');
A = A((firstloc:lastloc),:);

reals = load('reals.mat');
reals = reals.reals;
reals(:,3)=reals(:,3)-reals(:,5)/2;
reals(:,4)=reals(:,4)-reals(:,6)/2;
orgdata = reals;

deneme = orgdata(:,1);
firstloc = find(((deneme>targetvalue)&(deneme<nontargetvalue)),1,'first');
lastloc = find(((deneme>targetvalue)&(deneme<nontargetvalue)),1,'last');
orgdata = orgdata((firstloc:lastloc),:);
son = orgdata(end,1)-targetvalue;

for k=1:length(orgdata(:,1))
    if(orgdata(k,2)==(0))
        orgdata(k,2) = 1;
    elseif(orgdata(k,2)==(-1))
        orgdata(k,2) = 0;
    end
end

for i=1:2

    DataFolder = zeros(son,20);
    if(i==1)
        objscore = 0.05;
    elseif(i==2)
        objscore = 0.5;
    end
    
for k=1:length(A(:,1))
frameNumb = mod(A(k,1),100000);

if((DataFolder(frameNumb,1)==0)&&(A(k,2)>objscore))
DataFolder(frameNumb,1)=A(k,2);
DataFolder(frameNumb,2)=A(k,3);
DataFolder(frameNumb,3)=A(k,4);
DataFolder(frameNumb,4)=A(k,5);
DataFolder(frameNumb,5)=A(k,6);

elseif((DataFolder(frameNumb,1)~=0)&&(DataFolder(frameNumb,6)==0)&&(A(k,2)>objscore))
DataFolder(frameNumb,6)=A(k,2);
DataFolder(frameNumb,7)=A(k,3);
DataFolder(frameNumb,8)=A(k,4);
DataFolder(frameNumb,9)=A(k,5);
DataFolder(frameNumb,10)=A(k,6);
      
elseif((DataFolder(frameNumb,1)~=0)&&(DataFolder(frameNumb,6)~=0)&&(DataFolder(frameNumb,11)==0)&&(A(k,2)>objscore))
DataFolder(frameNumb,11)=A(k,2);
DataFolder(frameNumb,12)=A(k,3);
DataFolder(frameNumb,13)=A(k,4);
DataFolder(frameNumb,14)=A(k,5);
DataFolder(frameNumb,15)=A(k,6);

elseif((DataFolder(frameNumb,1)~=0)&&(DataFolder(frameNumb,6)~=0)&&(DataFolder(frameNumb,11)~=0)&&(DataFolder(frameNumb,16)==0)&&(A(k,2)>objscore))
DataFolder(frameNumb,16)=A(k,2);
DataFolder(frameNumb,17)=A(k,3);
DataFolder(frameNumb,18)=A(k,4);
DataFolder(frameNumb,19)=A(k,5);
DataFolder(frameNumb,20)=A(k,6);
end

if(k==length(A(:,1)))
    if(i==1)
        LM = DataFolder;
    elseif(i==2)
        HM = DataFolder;
    end
end

end

end

%% Track Creation

disp('Tracks Are Being Created')
clearvars -except LM HM foldnumb orgdata
warning off;

tlow = 0.05;
thigh = 0.5;
PRO = struct;
n = 2;
y = 0;
trackid = 0;

for p=1:5:16
if (HM(1,p)~=0)
y = y+1;
trackid = trackid+1;
PRO(1).N(y).T.a = HM(1,p);
PRO(1).N(y).T.b = [HM(1,[(p):(p+4)]),0];
PRO(1).N(y).T.c = [1,1];
PRO(1).N(y).T.d = trackid;
PRO(1).N(y).T.e = 1;
PRO(1).N(y).T.speedx = 0;
PRO(1).N(y).T.speedy = 0;
end
end

if(y==0)
for p=1:5:16
if (LM(1,p)~=0)
y = y+1;
trackid = trackid+1;
PRO(1).N(y).T.a = LM(1,p);
PRO(1).N(y).T.b = [LM(1,[(p):(p+4)]),0];
PRO(1).N(y).T.c = [1,1];
PRO(1).N(y).T.d = trackid;
PRO(1).N(y).T.e = 1;
PRO(1).N(y).T.speedx = 0;
PRO(1).N(y).T.speedy = 0;
end
end
end
    
for k=2:length(LM(:,1))
%% Part A %%
PRO(k)=PRO(k-1);
loy = 0;
ppcurrent = 0;
clear SM;
SM = zeros(100,12);
clear tempPRO;
tempPRO = struct;

if (y==0)
    for p=1:5:16
        if (HM(k,p)~=0)
            y = y+1;
            trackid = trackid+1;
            PRO(k).N(y).T.a = HM(k,p);
            PRO(k).N(y).T.b = [HM(k,[(p):(p+4)]),0];
            PRO(k).N(y).T.c = [0,k];
            PRO(k).N(y).T.d = trackid;
            PRO(k).N(y).T.e = 1;
            PRO(k).N(y).T.speedx = 0;
            PRO(k).N(y).T.speedy = 0;
        end
    end
    
    if (y==0)
    for p=1:5:16
        if (LM(k,p)~=0)
            y = y+1;
            trackid = trackid+1;
            PRO(k).N(y).T.a = LM(k,p);
            PRO(k).N(y).T.b = [LM(k,[(p):(p+4)]),0];
            PRO(k).N(y).T.c = [0,k];
            PRO(k).N(y).T.d = trackid;
            PRO(k).N(y).T.e = 1;
            PRO(k).N(y).T.speedx = 0;
            PRO(k).N(y).T.speedy = 0;
        end
    end
    end
    
    
elseif (y~=0)
    
    for m=1:y
        cmax = 0;
        pcount = 1;
        ps = 0;
        temparrcount = 0;
        
        for pp=1:loy
            if(PRO(k).N(m).T.b(1:5)==SM(pp,[6:10]))
                ps = 1;
                ppcurrent = pp;
                break
            end
        end
        
        %
        
        if(ps==0)
            loy = loy+1;
            checking = 0;
            
        for p=1:5:16
            if (LM(k,p)~=0)
                
            checking=1;
            temparr = PRO(k).N(m).T.b(1:5);
            t = PRO(k).N(m).T.c(2);
            ctempcolor = 0;
            foldernumber = foldnumb*100000;
            mapa = foldernumber+k;
            mapt = foldernumber+t;
            if(foldnumb>=10)
                strcatt = strcat(num2str(mapt),'.png');
                strcata = strcat(num2str(mapa),'.png');
            else
                strcatt = strcat('0',num2str(mapt),'.png');
                strcata = strcat('0',num2str(mapa),'.png');
            end
            trackss = imread(strcatt);
            alarmss = imread(strcata);
            
            for color=1:3

            alarm = alarmss(:,:,color);
            track = trackss(:,:,color);
        
            d3 = temparr(2);
            d4 = temparr(4);
            d1 = temparr(3);
            d2 = temparr(5);
            track = track((d1):(d2),(d3):(d4));
            
            y3 = LM(k,(p+1));
            y4 = LM(k,(p+3));
            y1 = LM(k,(p+2));
            y2 = LM(k,(p+4));
            ymid = (y1+y2)/2;
            xmid = (y3+y4)/2;
            
            ccxgap = d4 - d3;
            ccygap = d2 - d1;
            
            ctemp = 0;
            maxc = 0;
            c = 0;
            
            for rx=-2:2
                for ry=-2:2
                    
                    xmidtemp=xmid+rx;
                    ymidtemp=ymid+ry;
                    
                    ccxmin = xmidtemp - 2.5 - (ccxgap/2);
                    ccxmax = xmidtemp + 2.5 + (ccxgap/2);
                    ccymin = ymidtemp - 2.5 - (ccygap/2);
                    ccymax = ymidtemp + 2.5 + (ccygap/2);
                    alarms = alarm((max(ccymin,1)):(min(ccymax,360)),(max(ccxmin,1)):min(640,(ccxmax)));       
            
                if(size(track)<size(alarms))
                    if(range(track) ~= 0)
                        ctemp = normxcorr2(track,alarms);
                    else
                        ctemp = 0;
                    end
                else
                    ctemp = 0;
                end
                
                if(max(ctemp(:))>max(c(:)))
                    c = ctemp;
                end
                
                end
            end
            
            if(max(c(:))>max(ctempcolor(:)))
                ctempcolor = c;
            end
            
            end
            
            c = ctempcolor;
            
            if(max(c(:))>cmax)
                cmax = max(c(:));
                pcount = p;
                temparrcount = temparr; 
            end
            
            
            end
                
        end
        
        if((LM(k,1)==0)&&(LM(k,6)==0)&&(LM(k,11)==0)&&(LM(k,16)==0))
            cmax = 0;
            SM(loy,[1:5]) = LM(k,[(pcount):(pcount+4)]);
            SM(loy,[6:10]) = temparrcount;
            SM(loy,11) = cmax;
        end
        
        if(checking==1)
        SM(loy,[1:5]) = LM(k,[(pcount):(pcount+4)]);
        SM(loy,[6:10]) = temparrcount;
        SM(loy,11) = cmax;
        end
        
        %% Track Speed %%
        
        if(cmax>0)
            
            if(PRO(k).N(m).T.speedx == 0)
                PRO(k).N(m).T.speedx(1) = ((SM(loy,2)+SM(loy,4))/2)-((PRO(k).N(m).T.b(2)+PRO(k).N(m).T.b(4))/2);
                PRO(k).N(m).T.speedy(1) = ((SM(loy,3)+SM(loy,5))/2)-((PRO(k).N(m).T.b(3)+PRO(k).N(m).T.b(5))/2);
                trust = 1;
                
            elseif(size(PRO(k).N(m).T.speedx,2)<4)
                spe = size(PRO(k).N(m).T.speedx,2)+1;
                PRO(k).N(m).T.speedx(spe) = 0;
                PRO(k).N(m).T.speedy(spe) = 0;
                PRO(k).N(m).T.speedx = circshift(PRO(k).N(m).T.speedx,1);
                PRO(k).N(m).T.speedy = circshift(PRO(k).N(m).T.speedy,1);
                PRO(k).N(m).T.speedx(1) = ((SM(loy,2)+SM(loy,4))/2)-((PRO(k).N(m).T.b(2)+PRO(k).N(m).T.b(4))/2);
                PRO(k).N(m).T.speedy(1) = ((SM(loy,3)+SM(loy,5))/2)-((PRO(k).N(m).T.b(3)+PRO(k).N(m).T.b(5))/2);
                trust = 1;
                    
            elseif((size(PRO(k).N(m).T.speedx,2)>3) && (size(PRO(k).N(m).T.speedx,2)~=9))
                spe = size(PRO(k).N(m).T.speedx,2)+1;
                Vmeanx = mean(PRO(k).N(m).T.speedx);
                Vmeany = mean(PRO(k).N(m).T.speedy);
                sigmax = sqrt(var(PRO(k).N(m).T.speedx));
                sigmay = sqrt(var(PRO(k).N(m).T.speedy));
                PRO(k).N(m).T.speedx(spe) = 0;
                PRO(k).N(m).T.speedy(spe) = 0;
                PRO(k).N(m).T.speedx = circshift(PRO(k).N(m).T.speedx,1);
                PRO(k).N(m).T.speedy = circshift(PRO(k).N(m).T.speedy,1);
                PRO(k).N(m).T.speedx(1) = ((SM(loy,2)+SM(loy,4))/2)-((PRO(k).N(m).T.b(2)+PRO(k).N(m).T.b(4))/2);
                PRO(k).N(m).T.speedy(1) = ((SM(loy,3)+SM(loy,5))/2)-((PRO(k).N(m).T.b(3)+PRO(k).N(m).T.b(5))/2);
                
                if((abs(PRO(k).N(m).T.speedx(1)-Vmeanx)<(3*sigmax)) && (abs(PRO(k).N(m).T.speedy(1)-Vmeany)<(3*sigmay)))
                    trust = 1;
                else
                    PRO(k).N(m).T.speedx(1) = PRO(k).N(m).T.speedx(1)/10;
                    PRO(k).N(m).T.speedy(1) = PRO(k).N(m).T.speedy(1)/10;
                    trust = 0;
                end
                
            elseif(size(PRO(k).N(m).T.speedx,2)==9)
                Vmeanx = mean(PRO(k).N(m).T.speedx);
                Vmeany = mean(PRO(k).N(m).T.speedy);
                sigmax = sqrt(var(PRO(k).N(m).T.speedx));
                sigmay = sqrt(var(PRO(k).N(m).T.speedy));
                PRO(k).N(m).T.speedx = circshift(PRO(k).N(m).T.speedx,1);
                PRO(k).N(m).T.speedy = circshift(PRO(k).N(m).T.speedy,1);
                PRO(k).N(m).T.speedx(1) = ((SM(loy,2)+SM(loy,4))/2)-((PRO(k).N(m).T.b(2)+PRO(k).N(m).T.b(4))/2);
                PRO(k).N(m).T.speedy(1) = ((SM(loy,3)+SM(loy,5))/2)-((PRO(k).N(m).T.b(3)+PRO(k).N(m).T.b(5))/2);
                
                if((abs(PRO(k).N(m).T.speedx(1)-Vmeanx)<(3*sigmax)) && (abs(PRO(k).N(m).T.speedy(1)-Vmeany)<(3*sigmay)))
                    trust = 1;
                else
                    PRO(k).N(m).T.speedx(1) = PRO(k).N(m).T.speedx(1)/10;
                    PRO(k).N(m).T.speedy(1) = PRO(k).N(m).T.speedy(1)/10;
                    trust = 0;
                end
                
            end
            
        end
     
        %%

        if(size(PRO(k).N(m).T.a,2)<10)
        
        lokum = size(PRO(k).N(m).T.a,2)+1;    
        
        if(cmax>0.5) 
            
        if((cmax<0.9) && (trust==0))
        y = y+1;
        trackid = trackid+1;
        PRO(k).N(y)=PRO(k).N(m);
        PRO(k).N(y).T.a(lokum) = 0;
        PRO(k).N(y).T.a = circshift(PRO(k).N(y).T.a,1);
        PRO(k).N(y).T.c(1) = PRO(k).N(y).T.c(1)+1;
        PRO(k).N(y).T.d = trackid;
        PRO(k).N(y).T.e = 1;
        end
        PRO(k).N(m).T.a(lokum) = LM(k,pcount);
        PRO(k).N(m).T.a = circshift(PRO(k).N(m).T.a,1);
        tempPRO.N(m).T.b = PRO(k).N(m).T.b;
        tempPRO.N(m).T.c = PRO(k).N(m).T.c;
        PRO(k).N(m).T.b = [LM(k,[(pcount):(pcount+4)]),cmax];
        moly = PRO(k).N(m).T.c(1);
        PRO(k).N(m).T.c = [(moly+1),k];
        PRO(k).N(m).T.e = PRO(k).N(m).T.e+1;
        
        elseif(cmax<=0.5)
        PRO(k).N(m).T.a(lokum) = 0;
        PRO(k).N(m).T.a = circshift(PRO(k).N(m).T.a,1);
        PRO(k).N(m).T.c(1) = PRO(k).N(m).T.c(1)+1;
        end
            
        %
        
        elseif(size(PRO(k).N(m).T.a,2)==10)
             
        if(cmax>0.5)
             
        if((cmax<0.5) && (trust==0))
        y = y+1;
        trackid = trackid+1;
        PRO(k).N(y)=PRO(k).N(m);
        PRO(k).N(y).T.a = circshift(PRO(k).N(y).T.a,1);
        PRO(k).N(y).T.a(1) = 0;
        PRO(k).N(y).T.c(1) = PRO(k).N(y).T.c(1)+1;
        PRO(k).N(y).T.d = trackid;
        PRO(k).N(y).T.e = 1;
        end
        PRO(k).N(m).T.a = circshift(PRO(k).N(m).T.a,1);
        PRO(k).N(m).T.a(1) = LM(k,pcount);
        tempPRO.N(m).T.b = PRO(k).N(m).T.b;
        tempPRO.N(m).T.c = PRO(k).N(m).T.c;
        PRO(k).N(m).T.b = [LM(k,[(pcount):(pcount+4)]),cmax];
        moly = PRO(k).N(m).T.c(1);
        PRO(k).N(m).T.c = [(moly+1),k];
        PRO(k).N(m).T.e = PRO(k).N(m).T.e+1;
        
        elseif(cmax<=0.5)
        PRO(k).N(m).T.a = circshift(PRO(k).N(m).T.a,1);
        PRO(k).N(m).T.a(1) = 0;
        PRO(k).N(m).T.c(1) = PRO(k).N(m).T.c(1)+1;
        end
        
        end
        
        %
        
        elseif(ps==1)

        if(size(PRO(k).N(m).T.a,2)<10)

        if(SM(ppcurrent,11)>0.5)
            
        lokum = size(PRO(k).N(m).T.a,2)+1;
            
        if((SM(ppcurrent,11)<0.9) && (trust==0))
        y = y+1;
        trackid = trackid+1;
        PRO(k).N(y)=PRO(k).N(m);
        PRO(k).N(y).T.a(lokum) = 0;
        PRO(k).N(y).T.a = circshift(PRO(k).N(y).T.a,1);
        PRO(k).N(y).T.c(1) = PRO(k).N(y).T.c(1)+1;
        PRO(k).N(y).T.d = trackid;
        PRO(k).N(y).T.e = 1;
        end
        PRO(k).N(m).T.a(lokum) = LM(k,pcount);
        PRO(k).N(m).T.a = circshift(PRO(k).N(m).T.a,1);
        tempPRO.N(m).T.b = PRO(k).N(m).T.b;
        tempPRO.N(m).T.c = PRO(k).N(m).T.c;
        PRO(k).N(m).T.b = [SM(ppcurrent,[1:5]),SM(ppcurrent,11)];
        moly = PRO(k).N(m).T.c(1);
        PRO(k).N(m).T.c = [(moly+1),k];
        PRO(k).N(m).T.e = PRO(k).N(m).T.e+1;
        
        elseif(SM(ppcurrent,11)<=0.5)
        lokum = size(PRO(k).N(m).T.a,2)+1;
        PRO(k).N(m).T.a(lokum) = 0;
        PRO(k).N(m).T.a = circshift(PRO(k).N(m).T.a,1);
        PRO(k).N(m).T.c(1) = PRO(k).N(m).T.c(1)+1;
        end
        
        elseif(size(PRO(k).N(m).T.a,2)==10)
        
        if(SM(ppcurrent,11)>0.5)
            
        if((SM(ppcurrent,11)<0.9) && (trust==0))
        y = y+1;
        trackid = trackid+1;
        PRO(k).N(y)=PRO(k).N(m);
        PRO(k).N(y).T.a = circshift(PRO(k).N(m).T.a,1);
        PRO(k).N(y).T.a(1) = 0;
        PRO(k).N(y).T.c(1) = PRO(k).N(y).T.c(1)+1;
        PRO(k).N(y).T.d = trackid;
        PRO(k).N(y).T.e = 1;
        end
        PRO(k).N(m).T.a = circshift(PRO(k).N(m).T.a,1);
        PRO(k).N(m).T.a(1) = LM(k,pcount);
        tempPRO.N(m).T.b = PRO(k).N(m).T.b;
        tempPRO.N(m).T.c = PRO(k).N(m).T.c;
        PRO(k).N(m).T.b = [SM(ppcurrent,[1:5]),SM(ppcurrent,11)];
        moly = PRO(k).N(m).T.c(1);
        PRO(k).N(m).T.c = [(moly+1),k];
        PRO(k).N(m).T.e = PRO(k).N(m).T.e+1;
        
        elseif(SM(ppcurrent,11)<=0.5)
        PRO(k).N(m).T.a = circshift(PRO(k).N(m).T.a,1);
        PRO(k).N(m).T.a(1) = 0;
        PRO(k).N(m).T.c(1) = PRO(k).N(m).T.c(1)+1;
        end
        
        end
        
        end
        
    end
    
end

%% Merge Tracks %%
ycount = y;
Speedy = zeros(100,4);
mas = 1;

for m=1:y
    if(PRO(k).N(m).T.a(1)~=0)
    mass = 0;
    for mm=1:mas
        if(PRO(k).N(m).T.a(1)==Speedy(mm,1))
            mass = 1;
            ssam = mm;
            break
        end
    end
    if(mass==0)
        mas=mas+1;
        Speedy(mas,1) = PRO(k).N(m).T.a(1);
        Speedy(mas,2) = PRO(k).N(m).T.d;
        Speedy(mas,3) = mean(PRO(k).N(m).T.a);
        Speedy(mas,4) = m;
    elseif(mass==1)
        if(mean(PRO(k).N(m).T.a)<Speedy(ssam,3))
            PRO(k).N(m).T.d=(-1);
        elseif(mean(PRO(k).N(m).T.a)>Speedy(ssam,3))
            PRO(k).N(Speedy(ssam,4)).T.d=(-1);
            Speedy(ssam,3) = mean(PRO(k).N(m).T.a);
            Speedy(ssam,4) = m;
            if(PRO(k).N(m).T.e<11)
                PRO(k).N(m).T.e=11;
            end
        end
    end
    end
end

for m=1:ycount
    if(PRO(k).N(ycount-m+1).T.d==(-1))
       PRO(k).N(ycount-m+1) = [];
       y = y-1;
    end
end

%% Create New Tracks for Unmatched Tracks with Obj > 0.5 %%

for p=1:5:16 
    if (HM(k,p)~=0)
        check = 0;
        for m=1:y
            if(PRO(k).N(m).T.b(1:5) == HM(k,[(p):(p+4)]))
            check = 1;
            end
        end
        
        if(check==0)
            y = y+1;
            trackid = trackid+1;
            PRO(k).N(y).T.a = HM(k,p);
            PRO(k).N(y).T.b = [HM(k,[(p):(p+4)]),0];
            PRO(k).N(y).T.c = [1,k];
            PRO(k).N(y).T.d = trackid;
            PRO(k).N(y).T.e = 1;
            PRO(k).N(y).T.speedx = 0;
            PRO(k).N(y).T.speedy = 0;
        end
    end
end


%% If track average obj is below 0.3, delete it %%

port = 0;
portm = zeros(1,y);

for m=1:y
    if(size(PRO(k).N(m).T.a,2)==10)
        if(mean(PRO(k).N(m).T.a)<0.3)
            port = port+1;
            portm(port) = m;
        end
    end
end

for i=1:port
	yolx = portm(port-i+1);
    PRO(k).N(yolx) = [];
    y = y-1;
end

%% If track is not continueing for the last 10 frames, delete it %%

if(k>10)
    ycurrent = y;
    for i=1:ycurrent
        if(PRO(k).N(ycurrent-i+1).T.c(2)<(k-10))
            PRO(k).N(ycurrent-i+1) = [];
            y = y-1;
        end
    end
end

%% max 1000 track %%

ycurrent = y;

if (ycurrent>1000)
    
    sortm = zeros(1,ycurrent);
    
    for o=1:ycurrent
        sortm(o)= mean(PRO(k).N(o).T.a);
    end
    
    while(y>1000)
                [mm,nn] = min(sortm);
                PRO(k).N(nn) = [];
                sortm(nn) = [];
                y = y-1;
    end
end

%% Make a new track if no tracks are present but detections are %%

if (y==0)
    for p=1:5:16
        if (HM(k,p)~=0)
            y = y+1;
            trackid = trackid+1;
            PRO(k).N(y).T.a = HM(k,p);
            PRO(k).N(y).T.b = [HM(k,[(p):(p+4)]),0];
            PRO(k).N(y).T.c = [1,k];
            PRO(k).N(y).T.d = trackid;
            PRO(k).N(y).T.e = 1;
            PRO(k).N(y).T.speedx = 0;
            PRO(k).N(y).T.speedy = 0;
        end
    end
end

if (y==0)
    for p=1:5:16
        if (LM(k,p)~=0)
            y = y+1;
            trackid = trackid+1;
            PRO(k).N(y).T.a = LM(k,p);
            PRO(k).N(y).T.b = [LM(k,[(p):(p+4)]),0];
            PRO(k).N(y).T.c = [1,k];
            PRO(k).N(y).T.d = trackid;
            PRO(k).N(y).T.e = 1;
            PRO(k).N(y).T.speedx = 0;
            PRO(k).N(y).T.speedy = 0;
        end
    end
end
disp1 = strcat('Currently at frame number:',32,num2str(k));
disp(disp1)
end

%% Track Evaluation
disp('Tracks Are Being Processed')
clearvars -except PRO LM HM foldnumb orgdata

%% RedTrack

TM = struct;
RedTxt = zeros(length(LM(:,1)),6);
rd = 1;

for k=1:length(LM(:,1))
    for p=1:5:16
        TM(k).N(p).a = 0;
        TM(k).N(p).b = 0;
        TM(k).N(p).c = 0;
        TM(k).N(p).d = 0;
    end
end

foldernumber = foldnumb*100000;
for k=1:length(LM(:,1))
    mapa = foldernumber+k;
    if(foldnumb>=10)
        ImgRead = strcat(num2str(mapa),'.png');
    else
        ImgRead = strcat('0',num2str(mapa),'.png');
    end
    if(isfile(ImgRead)==1)
    if(size(PRO(k).N,2)~=0)
    MS1 = zeros(100,5);
    MS2 = 0;
    I = imread(ImgRead);
    for p=1:5:16
        y = size(PRO(k).N,2);
        m = 1;
            if(PRO(k).N(m).T.b(1:5) == LM(k,[(p):(p+4)]))
                if(PRO(k).N(m).T.e>10) 
                    if(TM(k).N(p).a==0)
                        TM(k).N(p).a = PRO(k).N(m).T.d;
                        TM(k).N(p).b = PRO(k).N(m).T.b(6);
                        TM(k).N(p).d = mean(PRO(k).N(m).T.a);
                    end
                elseif(PRO(k).N(m).T.e==10)
                    for i=1:k
                        mapayeni = foldernumber+i;
                        if(foldnumb>=10)
                            ImgReadyeni = strcat(num2str(mapayeni),'.png');
                        else
                            ImgReadyeni = strcat('0',num2str(mapayeni),'.png');
                        end
                        if(isfile(ImgReadyeni)==1)
                            if(size(PRO(i).N,2)~=0)
                                if((PRO(i).N(1).T.d == PRO(k).N(1).T.d)&&(PRO(i).N(1).T.a(1)~=0))
                                    RedTxt(rd,1) = foldernumber+i;
                                    RedTxt(rd,2) = PRO(i).N(1).T.d;
                                    RedTxt(rd,[3:6]) = PRO(i).N(1).T.b(2:5);
                                    rd = rd+1;
                                end
                            end
                        end
                    end
                end
            end
        
        if(TM(k).N(p).a ~= 0)
            RedTxt(rd,1) = mapa;
            RedTxt(rd,2) = PRO(k).N(1).T.d;
            RedTxt(rd,[3:6]) = PRO(k).N(1).T.b(2:5);
            rd = rd+1;
        end
        
    end
        
    end
    end
end

RedTxt = RedTxt([1:(rd-1)],:);

%% Before Click Track Performance

FinalRedMatrix = RedTxt;
FinalRedMatrix(:,2) = [];

denemetext = strcat('Video',num2str(foldnumb),'FinalRedText.txt');
if(isfile(denemetext)==1)
    DenemeEski = readmatrix(denemetext);
    if(isempty(DenemeEski)==0)
        RedTxtTempe = setdiff(FinalRedMatrix(:,1),DenemeEski(:,1),'rows');
        [~, rowsA, rowsB] = intersect(RedTxtTempe(:, 1), FinalRedMatrix(:, 1));
        FinalRedMatrix = FinalRedMatrix(rowsB,:);
    end
end

denemetext2 = strcat('Video',num2str(foldnumb),'FinalMissText.txt');
if(isfile(denemetext2)==1)
    DenemeEski2 = readmatrix(denemetext2);
    if(isempty(DenemeEski2)==0)
        rowMiss = [];
        for p=1:size(DenemeEski2,1)
            firstloc1 = find((FinalRedMatrix(:,1)==DenemeEski2(p,1)),1,'first');
            if(size(firstloc1,1)~=0)
                yd1 = FinalRedMatrix(firstloc1,2);
                yd3 = FinalRedMatrix(firstloc1,4)-FinalRedMatrix(firstloc1,2);
                yd2 = FinalRedMatrix(firstloc1,3);
                yd4 = FinalRedMatrix(firstloc1,5)-FinalRedMatrix(firstloc1,3);
                
                ed1 = DenemeEski2(p,2);
                ed3 = DenemeEski2(p,4)-DenemeEski2(p,2);
                ed2 = DenemeEski2(p,3);
                ed4 = DenemeEski2(p,5)-DenemeEski2(p,3);
                
                matrixYOLOy = [yd1 yd2 yd3 yd4];
                matrixYOLOe = [ed1 ed2 ed3 ed4];
                
                areaORGe = rectint(matrixYOLOe,matrixYOLOe);
                IOUd = rectint(matrixYOLOe,matrixYOLOy);
                IOUd = IOUd/areaORGe;
                
                if(IOUd<0.5)
                    rowMiss = [rowMiss;p];
                end
                
            end
        end
        DenemeEski2(rowMiss,:) = [];
        RedTxtTempe2 = setdiff(FinalRedMatrix(:,1),DenemeEski2(:,1),'rows');
        [~, rowsA2, rowsB2] = intersect(RedTxtTempe2(:, 1), FinalRedMatrix(:, 1));
        FinalRedMatrix = FinalRedMatrix(rowsB2,:);
    end
end

if(isfile(denemetext)==1)
    if(isempty(DenemeEski)==0)
        FinalRedMatrix  = sortrows([FinalRedMatrix; DenemeEski]);
    end
end

disp('Track Performance Before Clicking')

objectpresentfolder = zeros(length(orgdata(:,2)),1);
for k=1:length(orgdata(:,2))
    if(orgdata(k,2)>0)
        objectpresentfolder(k,1)=1;
    end
end
O = cumsum(objectpresentfolder);
frameswithobject = O(end);

tfa=0;
tious=0;
tiou=0;

for o=1:size(FinalRedMatrix,1)
    k = FinalRedMatrix(o,1)-foldernumber;
    deneme2 = orgdata(:,1);
    targetvalue2 = k+foldernumber;
    firstloc2 = find((deneme2==targetvalue2),1,'first');
    fa=0;
    iou=0;
    ious=0;
    matrixORG = orgdata(firstloc2,[3:6]);
    areaORG = rectint(matrixORG,matrixORG);
    
    if(orgdata(firstloc2,2)~=(-1))
        d1 = FinalRedMatrix(o,2);
        d3 = FinalRedMatrix(o,4)-FinalRedMatrix(o,2);
        d2 = FinalRedMatrix(o,3);
        d4 = FinalRedMatrix(o,5)-FinalRedMatrix(o,3);
        
        matrixYOLO = [d1 d2 d3 d4];
        areaYOLO = rectint(matrixYOLO,matrixYOLO);
        
        if(orgdata(firstloc2,2)~=0)
            IOU = rectint(matrixORG,matrixYOLO);
            IOU = IOU/(areaORG+areaYOLO-IOU);
        elseif(orgdata(firstloc2,2)==0)
            IOU = 0;
        end
        
        if(IOU>=0.2)
            ious=1;
            iou=1;
        elseif(IOU==0)
            fa=1;
        elseif(IOU<0.2)
            iou=1;
        end
        
        tious=tious+ious;
        tiou=tiou+iou;
        tfa=tfa+fa;
    end
end

tious = (tious)*100/(frameswithobject);
tiou = (tiou)*100/(frameswithobject);

textfile{1,1} = strcat('Hit Rate:',32,num2str(tious));
textfile{2,1} = strcat('Weak+Hit Rate:',32,num2str(tiou));
textfile{3,1} = strcat('False Alarm:',32,num2str(tfa));

finalperf = strcat('Video',num2str(foldnumb),'BeforeClickPerformanceText.txt');
writecell(textfile,finalperf);

clearvars -except PRO LM HM foldnumb orgdata RedTxt foldernumber

%% 

correlationScore = zeros(1,length(LM(:,1)));
objectnessScore = zeros(1,length(LM(:,1)));
avgObjectnessScore = zeros(1,length(LM(:,1)));
displacement = zeros(1,length(LM(:,1)));
bbox = zeros(length(LM(:,1)),4);

for k = 1:length(LM(:,1))
    if(size(PRO(k).N,2)~=0)
    correlationScore(k) = PRO(k).N(1).T.b(6);
    objectnessScore(k) = PRO(k).N(1).T.b(1);
    avgObjectnessScore(k) = mean(PRO(k).N(1).T.a);
    displacement(k) = sqrt((PRO(k).N(1).T.speedx(1))^2+(PRO(k).N(1).T.speedy(1))^2);
    bbox(k,:) = PRO(k).N(1).T.b([2:5]);
    else
    correlationScore(k) = 0;
    objectnessScore(k) = 0;
    avgObjectnessScore(k) = 0;
    displacement(k) = Inf;
    bbox(k,:) = [0,0,0,0];
    end
end

for k = 1:length(LM(:,1))
    mapa = foldernumber+k;
    if(foldnumb>=10)
        ImgRead = strcat(num2str(mapa),'.png');
    else
        ImgRead = strcat('0',num2str(mapa),'.png');
    end
    if(isfile(ImgRead)==0)
        correlationScore(k) = 1;
        objectnessScore(k) = 1;
        avgObjectnessScore(k) = 1;
        displacement(k) = 0;
    end
end

%% Distinct RedTxt

denemetext = strcat('Video',num2str(foldnumb),'FinalRedText.txt');
if(isfile(denemetext)==1)
    DenemeEski = readmatrix(denemetext);
    if(isempty(DenemeEski)==0)
        RedTxtTempe = setdiff(RedTxt(:,1),DenemeEski(:,1),'rows');
        [~, rowsA, rowsB] = intersect(RedTxtTempe(:, 1), RedTxt(:, 1));
        RedTxt = RedTxt(rowsB,:);
    end
end

%% Distinct MissTxt

denemetext2 = strcat('Video',num2str(foldnumb),'FinalMissText.txt');
if(isfile(denemetext2)==1)
    DenemeEski2 = readmatrix(denemetext2);
    if(isempty(DenemeEski2)==0)
        rowMiss = [];
        for p=1:size(DenemeEski2,1)
            firstloc1 = find((RedTxt(:,1)==DenemeEski2(p,1)),1,'first');
            if(size(firstloc1,1)~=0)
                yd1 = RedTxt(firstloc1,3);
                yd3 = RedTxt(firstloc1,5)-RedTxt(firstloc1,3);
                yd2 = RedTxt(firstloc1,4);
                yd4 = RedTxt(firstloc1,6)-RedTxt(firstloc1,4);
                
                ed1 = DenemeEski2(p,2);
                ed3 = DenemeEski2(p,4)-DenemeEski2(p,2);
                ed2 = DenemeEski2(p,3);
                ed4 = DenemeEski2(p,5)-DenemeEski2(p,3);
                
                matrixYOLOy = [yd1 yd2 yd3 yd4];
                matrixYOLOe = [ed1 ed2 ed3 ed4];
                
                areaORGe = rectint(matrixYOLOe,matrixYOLOe);
                IOUd = rectint(matrixYOLOe,matrixYOLOy);
                IOUd = IOUd/areaORGe;
                
                if(IOUd<0.5)
                    rowMiss = [rowMiss;p];
                end
                
            end
        end
        DenemeEski2(rowMiss,:) = [];
        RedTxtTempe2 = setdiff(RedTxt(:,1),DenemeEski2(:,1),'rows');
        [~, rowsA2, rowsB2] = intersect(RedTxtTempe2(:, 1), RedTxt(:, 1));
        RedTxt = RedTxt(rowsB2,:);
    end
end

%% Track Re-ID
clickDiff = 0;

for k=2:size(RedTxt,1)
    if((RedTxt(k,1)>(RedTxt((k-1),1)+9))&&(RedTxt(k,2)==RedTxt((k-1),2)))
        RedTxt([k:size(RedTxt,1)],2) = RedTxt([k:size(RedTxt,1)],2)+1000;
        clickDiff = clickDiff + 2;
    end
end

%%

TrackCountMatrix = RedTxt(:,2);
TrackCount = size(unique(TrackCountMatrix),1);

[~,out_first,~] = unique(TrackCountMatrix, 'first');

[~,out_last,~] = unique(TrackCountMatrix, 'last');

FinalEvaluationResult = [];
FinalclickIndex = 0;

for i = 1:TrackCount
    
    outfirst(i) = RedTxt(out_first(i),1)-foldernumber;
    outlast(i) = RedTxt(out_last(i),1)-foldernumber;
    
    [outputIndex] = PrepareTrackEvaluation(outfirst(i), outlast(i), correlationScore, objectnessScore, avgObjectnessScore, displacement);

    [evaluationResult,clickIndex] = SimulateEvaluateTrack(outputIndex,bbox,orgdata,foldnumb,RedTxt);
    
    FinalEvaluationResult = [FinalEvaluationResult,evaluationResult];
    FinalclickIndex = FinalclickIndex + clickIndex;
    
end

FinalclickIndex = FinalclickIndex - clickDiff;
RedTxtNew = [];

for i = 1:TrackCount
    
    cmos = RedTxt((out_first(i):out_last(i)),:);
    [rows, columns] = size(cmos);
    column1 = cmos(:, 1)';
    allNumbers = [cmos(1,1) : 1 : cmos(end, 1)];
    cmos_out = zeros(length(allNumbers), columns);
    cmos_out(((column1-foldernumber)-(column1(1)-foldernumber-1)),:) = cmos;
    RedTxtNew = [RedTxtNew;cmos_out];
    
end

RedTxt = RedTxtNew;
FinalEvaluationResult = FinalEvaluationResult';
FinalRedMatrix = [RedTxt,FinalEvaluationResult];

if(isempty(FinalRedMatrix)==1)
    disp('performance did not change compared to the previous iteration');
    return
end

denemematrix = FinalRedMatrix;

% Column Missmatch
TF3 = FinalRedMatrix(:,1) == 0;
FinalRedMatrix(TF3,:) = [];

TF4 = FinalRedMatrix(:,7) == 0;
FinalRedMatrix(TF4,:) = [];

% Miss
FinalRedMatrix(:,2) = [];
TF1 = FinalRedMatrix(:,6) == 3;
FinalMissMatrix = FinalRedMatrix(TF1,:);
%FinalMissMatrix = FinalMissMatrix(:,[1:5]);
FinalRedMatrix(TF1,:) = [];

% Weak
TF2 = FinalRedMatrix(:,6) == 2;
FinalRedMatrix(TF2,:) = [];

FinalRedMatrix = FinalRedMatrix(:,[1:5]);

%% Merge

% FinalRedMatrix
if(isfile(denemetext)==1)
    if(isempty(DenemeEski)==0)
        FinalRedMatrix  = sortrows([FinalRedMatrix; DenemeEski]);
    end
end
finalred = strcat('Video',num2str(foldnumb),'FinalRedText.txt');
writematrix(FinalRedMatrix,finalred);

% FinalMissMatrix
if(isfile(denemetext2)==1)
    if(isempty(DenemeEski2)==0)
        FinalMissMatrix  = sortrows([FinalMissMatrix; DenemeEski2]);
    end
end
finalmiss = strcat('Video',num2str(foldnumb),'FinalMissText.txt');
writematrix(FinalMissMatrix,finalmiss);

%% Track Performance
disp('Evaluating Track Performance')

objectpresentfolder = zeros(length(orgdata(:,2)),1);
for k=1:length(orgdata(:,2))
    if(orgdata(k,2)>0)
        objectpresentfolder(k,1)=1;
    end
end
O = cumsum(objectpresentfolder);
frameswithobject = O(end);

tfa=0;
tious=0;
tiou=0;

  for o=1:size(FinalRedMatrix,1)
        k = FinalRedMatrix(o,1)-foldernumber;
        deneme2 = orgdata(:,1);
        targetvalue2 = k+foldernumber;
        firstloc2 = find((deneme2==targetvalue2),1,'first');
        fa=0;
        iou=0;
        ious=0;
        matrixORG = orgdata(firstloc2,[3:6]);
        areaORG = rectint(matrixORG,matrixORG);
        
        if(orgdata(firstloc2,2)~=(-1))
                    d1 = FinalRedMatrix(o,2);
                    d3 = FinalRedMatrix(o,4)-FinalRedMatrix(o,2);
                    d2 = FinalRedMatrix(o,3);
                    d4 = FinalRedMatrix(o,5)-FinalRedMatrix(o,3);
                    
                    matrixYOLO = [d1 d2 d3 d4];
                    areaYOLO = rectint(matrixYOLO,matrixYOLO);
                    
                    if(orgdata(firstloc2,2)~=0)
                    IOU = rectint(matrixORG,matrixYOLO);
                    IOU = IOU/(areaORG+areaYOLO-IOU);
                    elseif(orgdata(firstloc2,2)==0)
                    IOU = 0;
                    end
                    
                    if(IOU>=0.2)
                        ious=1;
                        iou=1;
                    elseif(IOU==0)
                        fa=1;
                    elseif(IOU<0.2)
                        iou=1;
                    end
                    
                    tious=tious+ious;
                    tiou=tiou+iou;
                    tfa=tfa+fa;
        end
  end
  
  tious = (tious)*100/(frameswithobject);       
  tiou = (tiou)*100/(frameswithobject);
  
  textfile{1,1} = strcat('Hit Rate:',32,num2str(tious));
  textfile{2,1} = strcat('Weak+Hit Rate:',32,num2str(tiou));
  textfile{3,1} = strcat('False Alarm:',32,num2str(tfa));
  textfile{4,1} = strcat('Total Click Number:',32,num2str(FinalclickIndex));
  
  finalperf = strcat('Video',num2str(foldnumb),'PerformanceText.txt');
  writecell(textfile,finalperf);
    