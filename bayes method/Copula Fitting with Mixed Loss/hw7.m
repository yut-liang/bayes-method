clear
clc;
format long e;
addpath('..');
path_to_files = '.\';

%--- load from Excel file ----
%--- marginal distributions
[vals,nms]=xlsread('Distributions_and_correlations.xls','Distributions','C3:E1002');
InvPDF_W=sort(vals(:,1));
InvPDF_X=sort(vals(:,2));
InvPDF_Y=sort(vals(:,3));

%--- sum distribution
[vals,nms]=xlsread('Distributions_and_correlations.xls','Distributions','K3:K18');
InvPDF_Zwxy=sort(vals(:,1));
%------------------------------------------------

NW=length(InvPDF_W);
NX=length(InvPDF_X);
NY=length(InvPDF_Y);

NZwxy=length(InvPDF_Zwxy);

%% Construct Constraints

N=10;

NC=NZwxy; %number of constraints
matrix_constraint=zeros(NC,1+N^3+1+NC);
constraint_no=0;

IndicatorMatrix = zeros(N,N,N,NZwxy);

for k=1:(NZwxy)
    constraint_no=constraint_no+1;
   %
   %------------------ Calculate weights In Indicator matrix ---------- 
    Z = InvPDF_Zwxy(k);
    for h=1:N; h1=floor((h-1)/N*NW+1.5); h2=floor(h/N*NW+0.5);
        for i=1:N; i1=floor((i-1)/N*NX+1.5); i2=floor(i/N*NX+0.5);
            for j=1:N; j1=floor((j-1)/N*NY+1.5); j2=floor(j/N*NY+0.5);
                if (InvPDF_W(h2)+InvPDF_X(i2)+InvPDF_Y(j2)<=Z)
                    IndicatorMatrix(h,i,j,k)=1;
                elseif (InvPDF_W(h1)+InvPDF_X(i1)+InvPDF_Y(j1)>Z)
                    IndicatorMatrix(h,i,j,k)=0;
                else
                    IndicatorMatrix(h,i,j,k)=0.5;    %----- Stan's variant ----
                end
            end
        end
    end
    %}
    
    ConstraintRow=zeros(1, N^3);
    for h=1:N
        for i=1:N
            for j=1:N
                ConstraintRow((h-1)*N^2+(i-1)*N+j)=IndicatorMatrix(h,i,j,k);
            end
        end
    end
    matrix_constraint(constraint_no,1)=k;
    matrix_constraint(constraint_no,2:(N^3+1))=ConstraintRow/N;
    matrix_constraint(constraint_no,1+N^3+1)=-k/NZwxy;           % benchmark
end

matrix_sum=zeros(3*N,N^3+2);
for h=1:N
    for i=1:N
        for j=1:N
            matrix_sum(    h,(h-1)*N^2+(i-1)*N+j+1)=1;
            matrix_sum(  N+i,(h-1)*N^2+(i-1)*N+j+1)=1;
            matrix_sum(2*N+j,(h-1)*N^2+(i-1)*N+j+1)=1;
        end
    end
end

matrix_sum(:,1)=(1:(3*N))';
matrix_sum(:,N^3+2)=-ones(3*N,1);  %(-1/N)*ones(3*N,1);

%% Construct Problem

problem_statement=sprintf('%s\n', ...
'problem_findcopula_Case1_05, minimize', ...
'  0.7*meanabs_err(matrix_scenarios_05)', ...
'  0.3*meansquare_err(matrix_scenarios_05)', ...
'Constraint: = 0', ...
'  linearmulti(matrix_sum)', ...
'Box: >= 0');




header_matrix{1} = 'id';
for h=1:N
    for i=1:N
        for j=1:N; 
            header_matrix{1+(h-1)*N^2+(i-1)*N+j} = sprintf('h_%g_%g_%g',h,i,j); 
        end
    end
end
header_matrix{1+N^3+1} = 'scenario_benchmark';

clear iargstruc_arr;

count=1;
iargstruc_arr(count) = matrix_pack('matrix_scenarios_05', [matrix_constraint(:,1:(N^3)+1) -matrix_constraint(:,(N^3)+2)],header_matrix);
count = count + 1;
iargstruc_arr(count) = matrix_pack('matrix_sum', matrix_sum,header_matrix);

%------------ Export problem in TXT format -------
%mpsg_problem_exporttotext(problem_statement, iargstruc_arr, '.\');

%------------ Export problem into PSG_TOOLBOX -------
iargstruc_arr_new=tbpsg_convert_old_new(iargstruc_arr);
% Open Toolbox
tbpsg_toolbox(problem_statement, iargstruc_arr_new);

% Solve using Toolbox
%[a1,a2]=tbpsg_run(problem_statement, iargstruc_arr_new);
%out_structure = tbpsg_solution_struct(a1,a2);

%------------ Solve problem in Matlab ------
[solution_str, outargstruc_arr] = mpsg_solver(problem_statement, iargstruc_arr);
optimal_point = get_solution(solution_str, outargstruc_arr, 'optimal_point');
%opt_point_values=optimal_point(:,2);
disp('Optimal_point');
disp(optimal_point);
