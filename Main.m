clear;clc;close all
load net
load bp
clear input_test input_train output_train s s22 s222 train111 bn inputn input error outputn input an Ltest Ltrain minput moutput N1 N2 NN noutput
clear output_test label inputn_test ninput F test_simu0 test_simu s2 inputGMIMf  inputscps0 inputscps1 inputscps2 rowrank
clear error_train input0 inputGMIM1 inputGMIM2 inputGMIM3 inputGMIM4 inputGMIM5 inputGMIM6 inputGMIM7 inputscps3 w1 w2 w3 wij
filename='P.csv';
P=csvread(filename);%读取神经网络的输出元素
IT=0.2:0.2:10;
DIZJEMDONG=inputGMIM(:,2:end);
DIZJEMDONG=sortrows(DIZJEMDONG,1);
CHANGDI=sortrows(inputscps,1);
PGA=10.^(DIZJEMDONG(:,1));
vs=10.^(CHANGDI(:,1));
NCD=77900;
II=1:2000:100000;
[mi,ni]=size(II);
for j=1:ni
for i=1:50
    it=IT(i);
    P=[it,DIZJEMDONG(II(j),:),CHANGDI(NCD,:)];
    P=P';
    PP=mapminmax('apply',P,inputps,0,1);
    AA=sim(net,PP);
    test_simuAA(1,i)=mapminmax('reverse',AA,outputps);
    HVSR(i)=10.^(test_simuAA(1,i));
end
    HVSRRZ(j,:)=HVSR;
end
mhvsr=mean(HVSRRZ,1);
rems1=sum(HVSRRZ-mhvsr).^2;
rmse=sqrt(((1/(ni))*sum(HVSRRZ-mhvsr).^2));
qiangzhend=[DIZJEMDONG(109151,:);DIZJEMDONG(109203,:);DIZJEMDONG(110130,:);DIZJEMDONG(110880,:);DIZJEMDONG(112548,:);DIZJEMDONG(114603,:)
    ;DIZJEMDONG(117110,:);DIZJEMDONG(118110,:);DIZJEMDONG(118940,:);DIZJEMDONG(119770,:);DIZJEMDONG(120032,:)];

for j=1:11
for i=1:50
    it=IT(i);
    P=[it,qiangzhend(j,:),CHANGDI(NCD,:)];
    P=P';
    PP=mapminmax('apply',P,inputps,0,1);
    AA=sim(net,PP);
    test_simuAA(1,i)=mapminmax('reverse',AA,outputps);
    HVSR(i)=10.^(test_simuAA(1,i));
end
    HVSRQZ(j,:)=HVSR;
end
SSRNL=mhvsr./HVSRQZ;
DNLMN=sum(abs(log10(HVSRQZ)-log10(mhvsr)),2)*0.2;
PGVVS30=qiangzhend(:,7);
PGAQZ=10.^(qiangzhend(:,1));
NN=49;
moninihePGA=polyfit(log10(PGAQZ),DNLMN,1);%计算ANN给出的DNL拟合曲线,PGA
DNLNIHE1=1*moninihePGA(1,1)*log10(PGAQZ)+moninihePGA(1,2);

PGA0=(max(PGAQZ)-min(PGAQZ))/NN;
PGAz0=min(PGAQZ):PGA0:max(PGAQZ);
PGA2=interp1(PGAQZ,PGAQZ,PGAz0);
DNLNIHE1=1*moninihePGA(1,1)*log10(PGA2)+moninihePGA(1,2);
NN=49;
PGVVS300=(max(PGVVS30)-min(PGVVS30))/NN;
PGVVS30z0=min(PGVVS30):PGVVS300:max(PGVVS30);
PGVVS302=interp1(PGVVS30,PGVVS30,PGVVS30z0);
moninihePGVVS30=polyfit(log10(PGVVS30),DNLMN,1);
DNLNIHE2=1*moninihePGVVS30(1,1)*(log10((PGVVS302)))+moninihePGVVS30(1,2);

for i=1:NN
   FF(i)=log10(IT(i+1))-log10(IT(i));
end
FF(1,end+1)=FF(1,end);
ADNLMN=sum((abs(log10(HVSRQZ)-log10(mhvsr)).*FF),2);
moniADLPGA=polyfit(log10(PGAQZ),ADNLMN,1);
ADNLNIHE1=1*moniADLPGA(1,1)*log10(PGA2)+moniADLPGA(1,2);

moninihePGVVS30ADNL=polyfit(log10(PGVVS30),ADNLMN,1);
ADNLNIHE2=1*moninihePGVVS30ADNL(1,1)*(log10(PGVVS302))+moninihePGVVS30ADNL(1,2);
AAA=mhvsr.*FF;

for i=1:1
A2(i,:)=abs(trapz(log10(IT),log10((AAA(i,:)))));
end
PNL=100*ADNLMN./A2;
PNL(1,1)=PNL(2,1);
ex=['a*(tanh(log10((x))-b)+1)'];
PNLRen=23.77*(tanh(log(PGA2)-6.20)+1);
PNLreginer=13.82*(tanh(log10(PGA2)-3.56)+1);
startpoint=[1 2];
ft=fittype(ex);
PNLMH1=(61.84*(tanh(log10(PGA2)-3.377)+1))*1;

figure(1)
plot(log10(PGA2),DNLNIHE1)
title('PGA-DNL')
xlabel('log10(PGA)(cm/s^2)')
ylabel('DNL')
hold on
figure(2)
plot(log10(PGA2),ADNLNIHE1)
title('PGA-ADNL')
xlabel('log10(PGA)(cm/s^2)')
ylabel('ADNL')
hold on
figure(3)
plot(log10(PGA2),PNLMH1)
title('PGA-PNL')
xlabel('log10(PGA)(cm/s^2)')
ylabel('PNL')
hold on









