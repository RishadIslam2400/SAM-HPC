% load the data
data = load('bubble_exact_maps.mat');
exact_maps = data.exactMapArray;

sequence_of_A = load("sequence_of_A.mat");
sequence_of_A = sequence_of_A.matArray;

data = load("bubble_approx_maps_init.mat");
approxMaps = data.approxMapArray;

% for i = 1:numel(exact_maps)
%     diff = exact_maps{i} - approxMaps{i};
%     nrm = norm(diff, "fro");
% 
%     disp("Norm of the difference: ");
%     disp(nrm);
% end
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

plot(init_res_nrm, LineWidth=2, Marker="o", MarkerSize=12, LineStyle="-");
hold on;

data = load('bubble_approx_maps_simple.mat');
approxMaps = data.approxMapArray;
simple_res_nrm = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        simple_res_nrm(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(simple_res_nrm(i - 1));
    end
end

plot(simple_res_nrm, LineWidth=2, LineStyle="--", Marker="x", MarkerSize=10);

data = load('bubble_approx_maps_global_squared_0.1.mat');
approxMaps = data.approxMapArray;
global_res_nrm_1 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        global_res_nrm_1(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(global_res_nrm_1(i - 1));
    end
end

plot(global_res_nrm_1, LineWidth=2, LineStyle="-.", Marker="diamond", MarkerSize=8);

data = load('bubble_approx_maps_global_squared_0.01.mat');
approxMaps = data.approxMapArray;
global_res_nrm_2 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        global_res_nrm_2(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(global_res_nrm_2(i - 1));
    end
end

plot(global_res_nrm_2, LineWidth=2, Marker=".", MarkerSize=8);

data = load('bubble_approx_maps_global_squared_0.0001.mat');
approxMaps = data.approxMapArray;
global_res_nrm_3 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        global_res_nrm_3(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(global_res_nrm_3(i - 1));
    end
end

plot(global_res_nrm_3, LineWidth=2, Marker="o", MarkerSize=8);

data = load('bubble_approx_maps_col_squared_0.8.mat');
approxMaps = data.approxMapArray;
col_res_nrm_1 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        col_res_nrm_1(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(col_res_nrm_1(i - 1));
    end
end

plot(col_res_nrm_1, LineWidth=2, LineStyle="--", Marker="v", MarkerSize=8);

data = load('bubble_approx_maps_col_squared_0.6.mat');
approxMaps = data.approxMapArray;
col_res_nrm_2 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        col_res_nrm_2(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(col_res_nrm_2(i-1));
    end
end

plot(col_res_nrm_2, Marker="pentagram", LineWidth=2, MarkerSize=8);

data = load('bubble_approx_maps_col_squared_0.7.mat');
approxMaps = data.approxMapArray;
col_res_nrm_3 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        col_res_nrm_3(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(col_res_nrm_3(i-1));
    end
end

plot(col_res_nrm_3, Marker=">", MarkerSize=8, LineWidth=2, LineStyle="-.");

% data = load('bubble_approx_maps_lfil_squared_3.mat');
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

data = load('bubble_approx_maps_lfil_squared_5.mat');
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

data = load('bubble_approx_maps_lfil_squared_7.mat');
approxMaps = data.approxMapArray;
lfil_res_nrm_3 = zeros(numel(sequence_of_A) - 1, 1);
for i = 1:numel(sequence_of_A)
    if i == 1
        A0 = sequence_of_A{i};
        A0_nrm = norm(A0, "fro");
    else
        residual = sequence_of_A{i} * approxMaps{i-1} - A0;
        lfil_res_nrm_3(i - 1) = norm(residual, "fro") / A0_nrm;
        disp("Residual Norm of the difference: ");
        disp(lfil_res_nrm_3(i - 1));
    end
end

plot(lfil_res_nrm_3, Marker="<", LineWidth=2, MarkerSize=8, LineStyle="--");
hold off;
lgd = legend("$S(A_0)$", "$S(A_k)$", "$S(A_k)$, global $thresh = 0.1$", "$S(A_k)$, global $thresh = 0.01$", "$S(A_k)$, global $thresh = 0.0001$", "$S(A_k)$, column threshold $\tau = 0.8$", "$S(A_k)$, column threshold $\tau = 0.6$", "$S(A_k)$, column threshold $\tau = 0.7$", "$S(A_k)$, fixed non-zero $lfil = 5$", "$S(A_k)$, fixed non-zero $lfil = 7$");
lgd.FontSize = 12;
lgd.Interpreter = "latex";
lgd.Location = "best";
lgd.NumColumns = 2;
title("Relative Residual Plot for NWP Model", "FontSize",16)
xlabel("Matrices in the sequence","FontSize",14)
ylabel("Relative Residual Norm $\frac{||R_k||_F}{||A_0||_F}$","Interpreter","latex");