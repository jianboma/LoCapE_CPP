


% this script is used to calibrate the C++ code
% and show that C++ implementation has generate the same results as
% the algorithm implemented in matlab.


%% 
% demostrate using exponentional signal
% % 
% % clear all;
% % clc;
% % 
% % 
% % % s = rng;
% % % s=load('seed.mat');
% % rng(s.s);
% % f1 = 0.1;
% % N = 1024;
% % n = 0:N-1;
% % snr = 25;
% % sigma = sqrt(0.5*10^(-snr/10));
% % x0 = exp(1i*2*pi*f1*n)+sigma*(randn(1,N)+1i*randn(1,N));
% % writedata(x0, './exponential.txt');
% % 
% % 
% % 
% % 
% % pn = 8;
% % fN = pn*N;
% % n = 0:N-1;
% % % figure;plot(0:1/fN:1-1/fN,abs(fft(x0, pn*length(x0))));
% % fv = 0:1/fN:1-1/fN;
% % 
% % %specify damping grid
% % 
% % eta = (0:0.00005:0.0005);
% % eta_N = length(eta);
% % 
% % 
% % 
% % % start estimating spectrum
% % a_cl_c1 = zeros(eta_N,fN); a_cl_c2 = zeros(eta_N,fN);
% % a_cl_c_all1 = zeros(1,fN);
% % %set localised size 2*step + 1
% % p = 20;
% % 
% % fsta = 0.09; fend = 0.11;
% % jL = floor(fsta*fN)+1;
% % jH = ceil((fend*fN))+1 ;
% % [~, ~, a_cl_c1,~, ~] = LoCapE(x0,fsta,fend,p,floor(N/2),fN,eta);
% % a_cl_c_all1(jL:jH) = max(abs(a_cl_c1));   
% % 
% % 
% % a_cl_c_all1 = a_cl_c_all1./max(a_cl_c_all1);
% % 
% % 
% % 
% % % figure;plot(fv(jL:jH),max(abs(a_cl_c1)));
% % [a,ind] =max(abs(a_cl_c1));
% % for n=1:length(ind)
% %     aa(n) = a_cl_c1(ind(n),n);
% % end
% % figure;plot(fv(jL:jH),real(aa),'--o'); hold on;plot(fv(jL:jH),imag(aa),'--o');
% % hold on;
% % 
% % % data calculated by C++ code
% % command = ['your command'];
% % % system(command);
% % [data_real, data_imag] = read_data('./data/spectrum.txt');
% % plot(fv(jL:jH),data_real,'-x'); hold on;plot(fv(jL:jH),data_imag,'-x');
% % axis([fsta fend -0.05 0.05]);




%%
% demostrate using simulated NMR data

% clear all;
% clc;

% load the signal
data = load('../../data/ss1_8.mat'); 
x0 = data.data;
x0 = x0(1:64*1024); % 64Kb is used in Mnova
N = length(x0);
snr = 25;
% s = rng;
s=load('../../seed_analysis_different_slicewindow.mat');
rng(s.s);
sigma = sqrt(0.5*(x0*x0')/N*10^(-snr/10));
x0 = x0+sigma*(randn(1,N)+1i*randn(1,N));
% writedata(x0, './data/Simulated_NMR.txt');


fs = 500.13*10^6;
fl = -3000.78;
tvec = 0:1/fs:length(x0)*1/fs;

pn = 2;
fN = pn*N;
n = 0:N-1;
fv = 0:1/fN:1-1/fN;

%specify damping grid
nu = 1; % first moment

% eta = (0:0.0001:0.002);
eta = (0:0.00005:0.0005);
% eta = (0:0.00002:0.001);
eta_N = length(eta);

%set localised size 2*step + 1
p = 20;

fsta = 0.2465;
fend = 0.2532;
jL = floor(fsta*fN)+1;
jH = ceil((fend*fN))+1 ;
% start estimating spectrum
a_cl_c1 = zeros(eta_N,jH-jL+1);
datam = zeros(1,jH-jL+1);
tic
[~, ~, a_cl_c1,~, ~] = LoCapE(x0,fsta,fend,p,floor(N/2),fN,eta);
toc
[~, ind] = max(a_cl_c1);
for num = 1:length(ind)
    datam(num) = a_cl_c1(ind(num),num);
end
datam = datam./max(datam);
figure;plot(fv(jL:jH),real(datam),'--o');% hold on;plot(fv(jL:jH),imag(datam),'--o');
hold on;

% data calculated by C++ code
command = ['your command'];
% system(command);
[data_real, data_imag] = read_data('./data/spectrum_simulated_NMR.txt');
plot(fv(jL:jH),data_real,'-x'); %hold on;plot(fv(jL:jH),data_imag,'-x');
axis([fsta fend -0.5 1.2]);
grid on; 


%%
function writedata(x, wfilename)
fid = fopen(wfilename, 'wt');

for num = 1:length(x)
    fprintf(fid,'%f %f\n',real(x(num)), imag(x(num))); %, imag(data(num))
end
fclose(fid);
end

%% read data
function [data_real, data_imag] = read_data(xfilename)
fid = fopen(xfilename,'rt');
C = textscan(fid,'%s %s');
fclose(fid);

data_real = [];
data_imag = [];
for num = 1:size(C{1,1},1)
    data_real = [data_real;str2num(C{1,1}{num})];
    data_imag = [data_imag;str2num(C{1,2}{num})];
end
if nargout == 0
    figure;plot(data_real);hold on; plot(data_imag);
end

end
