clear;
clc;
%C:\Users\dell\Desktop\Numerical\project\input.xls
path = input('Input path of CMB data: ','s');
concentrationData = readtable(path,'sheet',1);
[row1,column1] = size(concentrationData);
concentration = table2array(concentrationData(1:row1,3:column1));

uncertaintyData = readtable(path,'sheet',2);
[row2,column2] = size(uncertaintyData);
sigmaConcentration = table2array(uncertaintyData(1:row2,3:column2));

profileData = readtable(path,'sheet',3);
[row3,column3] = size(profileData);
profile = table2array(profileData(1:row3,2:column3));

[row4, column4] = size(concentration);
% calculate sigma sigma*F sigma*C
%row4 = 3; % CHANGE THERE NOT USE ALL DATA TO TEST!!!
T = zeros(1, row4);
MyResult = [];
error = [];
for i = 1:row4
    %sigmaConcentration(i,:) C = concentration(i,:)
    sigma = zeros(column4);
    for j = 1:column4
        sigma(j,j) = sigmaConcentration(i,j);
    end
    F = profile';
    C = concentration(i,:)';
    A = sigma*F;
    b = sigma*C;
    % calculate ||Ax - b||min
    tic
    %A
    %b
    %[x] = lsqnonneg(A,b);
    [x]=nnls2(A,b);
    MyResult = [MyResult,x];
    e = norm(A*x - b,2);
    error = [error,e];
    T(i) = toc;
end
tMul = sum(T)
norm(error)