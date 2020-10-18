rand('seed', 0);
randn('seed', 0);

n = 2;
m = 200;
N = m/2;
M = m/2;

% positive examples
Y = [1.5+0.9*randn(1,0.6*N), 1.5+0.7*randn(1,0.4*N);
2*(randn(1,0.6*N)+1), 2*(randn(1,0.4*N)-1)];

% negative examples
X = [-1.5+0.9*randn(1,0.6*M),  -1.5+0.7*randn(1,0.4*M);
2*(randn(1,0.6*M)-1), 2*(randn(1,0.4*M)+1)];

x = [X Y];
y = [ones(1,N) -ones(1,M)];
%y=[1 -1 1 -1 1 -1 1 1 -1 -1]
ANoLable =  -((ones(n,1)*y).*x)';
A = [ -((ones(n,1)*y).*x)' -y']
xdat = x';
lambda = 0.05;

% partition the examples up in the worst possible way
% (subsystems only have positive or negative examples)
p = zeros(1,m);
p(y == 1)  = sort(randi([1 10], sum(y==1),1));
p(y == -1) = sort(randi([11 20], sum(y==-1),1));

%size(A)
fileID1= fopen('A.dat','w');
fprintf(fileID1, '%i 1:%d 2:%d\n', A(:,[3,1,2])');
fclose(fileID1);

fileID2= fopen('p.dat','w');
fprintf(fileID2, '%d\n', p);
fclose(fileID2);

% AA=zeros(m,n+1);
% fileID3= fopen('A.txt','r');
% formatSpec = '%d';
% for i=1:m
%     for j=1:(n+1)
%         AA(i,j)= fscanf(fileID3, formatSpec);
%     end
% end
% fclose(fileID3);
% AA
% fileID4= fopen('p.txt','r');
% formatSpec = '%d';
% pp=fscanf(fileID4, formatSpec)
% fclose(fileID4);

fileID3= fopen('output.txt','w');

[x history] = linear_svm(A, lambda, p, 1.0, 1, fileID3);
testLable = (ANoLable*x([1,2]) + x(3))>0;
hits = sum(testLable == (y>0)')
fclose(fileID3);

