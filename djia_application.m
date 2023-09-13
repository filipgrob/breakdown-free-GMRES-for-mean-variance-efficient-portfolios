clear all
close all
clc

%main script replicating results of two articles on
%a variant of gmres applied to a financial problem. 
%Useful links about the method, the application in
%finance and the dataset are reported below:
%https://www.sciencedirect.com/science/article/pii/S0165188912000991
%https://www.math.kent.edu/~reichel/publications/bfgmres.pdf
%https://www.kaggle.com/datasets/szrlee/stock-time-series-20050101-to-20171231?resource=download

%import data (stock data for dija 30 companies during 2017)
data = readtable('all_stocks_2017-01-01_to_2018-01-01.csv');
prices=table2array(data(:,2));
prices=reshape(prices,251,31);
idx=find(any(isnan(prices),2));
prices(idx,:)=[];
returns=zeros(248,31);

%returns computation
for j=1:31
    for i=1:248
        returns(i,j)=((prices(i+1,j)-prices(i,j))/prices(i,j));
    end
end

%covariance computation until november (excluded)
%indexes <=208--->observation before november 2017
%indexes in [209:230]---> november 2017 observations
%indexes >230--->december 2017 observations
mu=zeros(31,1);
cov=zeros(31,31);
for i=1:31
    mu(i)=mean(returns(1:208,i));
end
for i=1:31
    for j=1:31
        cov(i,j)=(sum((returns(1:208,i)-mu(i)).*(returns(1:208,j)-mu(j))))/207;
    end
end

%computation of minimum variance portfolio weights and returns
tol=1e-8;
maxit=100;
[x]=bfgmres(cov,ones(31,1),tol,maxit);
xstar=x/norm(x);

%I decided to observe returns until november excluded and then compute
%minimum variance portfolio weights. This portfolio is then adopted for
%the following two months. Returns on november and december are observed 
%and compared to the market return that is represented by djia index.
open_nov17=23334.00;
open_dic17=24235.00;
divisor_nov17=sum(mean(prices(209:230,:)))/open_nov17; %find november divisor of djia index
divisor_dic17=sum(mean(prices(231:end,:)))/open_dic17; %find december divisor of djia index
djia_nov=mean(returns(209:230,:),2)/divisor_nov17;
djia_dic=mean(returns(231:end,:),2)/divisor_dic17;
djia=[djia_nov;djia_dic];
minvar=returns(209:end,:)*xstar;

%plots of minimum variance portfolio and market portfolio
figure;
plot(djia,'red');
hold on
plot(minvar,'blue');
legend('Market portfolio returns (DJIA)', 'MinVar portfolio returns',...
'Location', 'southwest');
title('DJIA and minimum variance portfolio comparison');
xlabel('November & Dicember 2017');
ylabel('Returns');


