% load the data
data = load('cd2d_exact_maps.mat');
exact_maps = data.exactMapArray;

sequence_of_A = load('cd2d_matrices.mat');
sequence_of_A = sequence_of_A.matrices;

data = load("cd2d_approx_maps_init.mat");
approxMaps = data.approxMapArray;
init_res_nrm = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        init_res_nrm(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(init_res_nrm(i - 1));
    end
end

plot(init_res_nrm, LineWidth=2, Marker="+", MarkerSize=8);
hold on;

data = load("cd2d_approx_maps_init_squared.mat");
approxMaps = data.approxMapArray;
init_res_nrm_cubed = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        init_res_nrm_cubed(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(init_res_nrm_cubed(i - 1));
    end
end

plot(init_res_nrm_cubed, LineWidth=2, Marker="+", MarkerSize=8);

% data = load('cd2d_approx_maps_global_squared_0.01.mat');
% approxMaps = data.approxMapArray;
% global_res_nrm_1 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%         A0_nrm = norm(A0, "fro");
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         global_res_nrm_1(i - 1) = norm(residual, "fro") / A0_nrm;
%         disp("Residual Norm of the difference: ");
%         disp(global_res_nrm_1(i - 1));
%     end
% end
% 
% plot(global_res_nrm_1, LineWidth=2, LineStyle="-.", Marker="diamond", MarkerSize=8);

% data = load('cd2d_approx_maps_global_squared_0.0001.mat');
% approxMaps = data.approxMapArray;
% global_res_nrm_2 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%         A0_nrm = norm(A0, "fro");
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         global_res_nrm_2(i - 1) = norm(residual, "fro") / A0_nrm;
%         disp("Residual Norm of the difference: ");
%         disp(global_res_nrm_2(i - 1));
%     end
% end
% 
% plot(global_res_nrm_2, LineWidth=2, Marker="o", MarkerSize=8);

% data = load('cd2d_approx_maps_global_squared_0.001.mat');
% approxMaps = data.approxMapArray;
% global_res_nrm_3 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%         A0_nrm = norm(A0, "fro");
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         global_res_nrm_3(i - 1) = norm(residual, "fro") / A0_nrm;
%         disp("Residual Norm of the difference: ");
%         disp(global_res_nrm_3(i - 1));
%     end
% end
% 
% plot(global_res_nrm_3, LineWidth=2, Marker="hexagram", MarkerSize=8, LineStyle="--");

% data = load('cd2d_approx_maps_global_cubed_0.01.mat');
% approxMaps = data.approxMapArray;
% global_res_nrm_4 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%         A0_nrm = norm(A0, "fro");
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         global_res_nrm_4(i - 1) = norm(residual, "fro") / A0_nrm;
%         disp("Residual Norm of the difference: ");
%         disp(global_res_nrm_4(i - 1));
%     end
% end
% 
% plot(global_res_nrm_4, LineWidth=2, LineStyle="-.", Marker="diamond", MarkerSize=8);

% data = load('cd2d_approx_maps_column_squared_0.5.mat');
% approxMaps = data.approxMapArray;
% col_res_nrm_1 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         col_res_nrm_1(i - 1) = norm(residual, "fro");
%         disp("Residual Norm of the difference: ");
%         disp(col_res_nrm_1(i - 1));
%     end
% end
% 
% plot(col_res_nrm_1, LineWidth=2, LineStyle="--", Marker="v", MarkerSize=8);

% data = load('cd2d_approx_maps_column_squared_0.6.mat');
% approxMaps = data.approxMapArray;
% col_res_nrm_2 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%         A0_nrm = norm(A0, "fro");
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         col_res_nrm_2(i - 1) = norm(residual, "fro") / A0_nrm;
%         disp("Residual Norm of the difference: ");
%         disp(col_res_nrm_2(i-1));
%     end
% end
% 
% plot(col_res_nrm_2, Marker="pentagram", LineWidth=2, MarkerSize=8);

% data = load('cd2d_approx_maps_column_squared_0.7.mat');
% approxMaps = data.approxMapArray;
% col_res_nrm_3 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%         A0_nrm = norm(A0, "fro");
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         col_res_nrm_3(i - 1) = norm(residual, "fro") / A0_nrm;
%         disp("Residual Norm of the difference: ");
%         disp(col_res_nrm_3(i-1));
%     end
% end
% 
% plot(col_res_nrm_3, Marker=">", MarkerSize=8, LineWidth=2, LineStyle="-.");

data = load('cd2d_approx_maps_column_squared_0.8.mat');
approxMaps = data.approxMapArray;
col_res_nrm_4 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        col_res_nrm_4(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(col_res_nrm_4(i-1));
    end
end

plot(col_res_nrm_4, Marker=">", MarkerSize=8, LineWidth=2, LineStyle="-.");

data = load('cd2d_approx_maps_column_cubed_0.8.mat');
approxMaps = data.approxMapArray;
col_res_nrm_5 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        col_res_nrm_5(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(col_res_nrm_5(i-1));
    end
end

plot(col_res_nrm_5, Marker=">", MarkerSize=8, LineWidth=2, LineStyle="-.");

% data = load('cd2d_approx_maps_lfil_squared_3.mat');
% approxMaps = data.approxMapArray;
% lfil_res_nrm_1 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%         A0_nrm = norm(A0, "fro");
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         lfil_res_nrm_1(i - 1) = norm(residual, "fro") / A0_nrm;
%         disp("Residual Norm of the difference: ");
%         disp(lfil_res_nrm_1(i-1));
%     end
% end
% 
% plot(lfil_res_nrm_1, Marker="*", LineWidth=2, LineStyle="--", MarkerSize=8);
% 
% data = load('cd2d_approx_maps_lfil_cubed_3.mat');
% approxMaps = data.approxMapArray;
% lfil_res_nrm_3 = zeros(numel(sequence_of_A) - 1, 1);
% for i = 1:numel(sequence_of_A)
%     if i == 1
%         A0 = sequence_of_A{i};
%         A0_nrm = norm(A0, "fro");
%     else
%         residual = sequence_of_A{i} * approxMaps{i-1} - A0;
%         lfil_res_nrm_3(i - 1) = norm(residual, "fro") / A0_nrm;
%         disp("Residual Norm of the difference: ");
%         disp(lfil_res_nrm_3(i-1));
%     end
% end
% 
% plot(lfil_res_nrm_3, Marker="*", LineWidth=2, LineStyle="--", MarkerSize=8);

data = load('cd2d_approx_maps_lfil_squared_5.mat');
approxMaps = data.approxMapArray;
lfil_res_nrm_2 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        lfil_res_nrm_2(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(lfil_res_nrm_2(i - 1));
    end
end

plot(lfil_res_nrm_2, Marker="hexagram", LineWidth=2, MarkerSize=8);

data = load('cd2d_approx_maps_lfil_cubed_5.mat');
approxMaps = data.approxMapArray;
lfil_res_nrm_4 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        lfil_res_nrm_4(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(lfil_res_nrm_4(i - 1));
    end
end

plot(lfil_res_nrm_4, Marker="hexagram", LineWidth=2, MarkerSize=8);
hold off;
title("Residual Plot for 2D Convection Diffusion Equation", "FontSize", 16);
lgd = legend("$S(A_0)$", "$S(A_0)$, level 1", "column threshold $\tau = 0.8$, level 1", "column threshold $\tau = 0.8$, level 2", "fixed non-zero $lfil = 5$, level 1", "fixed non-zero $lfil = 5$, level 2");
lgd.FontSize = 12;
lgd.Interpreter = "latex";
lgd.Location = "best";
lgd.NumColumns = 2;
xlabel("Matrices in the sequence","FontSize",14)
ylabel("Relative Residual Norm $\frac{||R_k||_F}{||A_0||_F}$","Interpreter","latex", FontSize=14);