close all

%Test 0
load('file64_0','-mat'); load('file256_0','-mat'); load('file1024_0','-mat'); load('file4096_0','-mat');


%Test 1
load('file64_1','-mat'); load('file256_1','-mat'); load('file1024_1','-mat'); load('file4096_1','-mat');

%Test 2
load('file64_2','-mat'); load('file256_2','-mat'); load('file1024_2','-mat'); load('file4096_2','-mat');


plot(f64_0(:,1),f64_0(:,2))
hold on
plot(f64_0(:,1),f256_0(:,2),'r')
hold on
plot(f64_0(:,1),f1024_0(:,2),'g')
hold on
plot(f64_0(:,1),f4096_0(:,2),'c')
legend('N=64','N=256','N=1024', 'N=4096');
hold off


figure;

plot(f64_0(:,1),f64_1(:,2))
hold on
plot(f64_0(:,1),f256_1(:,2),'r')
hold on
plot(f64_0(:,1),f1024_1(:,2),'g')
hold on
plot(f64_0(:,1),f4096_1(:,2),'c')
hold off
legend('N=64','N=256','N=1024', 'N=4096');
figure;


plot(f64_0(:,1),f64_2(:,2))
hold on
plot(f64_0(:,1),f256_2(:,2),'r')
hold on
plot(f64_0(:,1),f1024_2(:,2),'g')
hold on
plot(f64_0(:,1),f4096_2(:,2),'c')
hold off
legend('N=64','N=256','N=1024', 'N=4096');