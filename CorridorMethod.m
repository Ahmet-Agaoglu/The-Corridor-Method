function [V,cntrd,optSol]=CorridorMethod(x,y,Lb,Ub)
% CorridorMethod: Determines the Feasible Parameter Set (FPS) for a line fitting problem,
% computes the centroid (for unknown error distributions), and finds the optimal 
% solution (for normally distributed errors).

% INPUT:
% x  : Independent variable values
% y  : Dependent variable values
% Lb : Lower bounds of the dependent variable values
% Ub : Upper bounds of the dependent variable values

% OUTPUT:
% V       : Vertices of the FPS (2D matrix, each row is a vertex [intercept, slope])
% cntrd   : Centroid of the FPS (2-element vector [intercept, slope])
% optSol  : Optimal solution (2-element vector [intercept, slope])

% Initialization
% Sort values with respect to x.
Mat = [x(:),y(:),Lb(:),Ub(:)];
Mat = sortrows(Mat,1,"ascend");
x = Mat(:,1); y = Mat(:,2); Lb = Mat(:,3); Ub = Mat(:,4);

%% Main Script
% Determine the vertices of the lower bound (VL)
[VL,a,b] = findvertices_VL(x,Lb,Ub);
% Determine the vertices of the upper bound (VU)
[VU,c,d] = findvertices_VU(x,Lb,Ub);
% Find the extreme lines (lmin and lmax) enclosing the feasible parameter set
[lmin,lmax] = find_lminmax(x,Lb,Ub,VL,VU,a,b,c,d);
% Construct the complete vertex list V for the FPS
V = [VL;lmin;VU;lmax];

% Define the polygon in the parameter space (FPS)
S = polyshape(V(:,1),V(:,2));
% Find the centroid of the polygon (used when distributional information is unavailable)
[B0_cen,B1_cen] = centroid(S);
cntrd = [B0_cen,B1_cen]';
% Find the optimal solution (for normally distributed errors)
[B0_hat,B1_hat]=find_optimal(V,x,y,Lb,Ub);
optSol = [B0_hat,B1_hat]';

% Plot the results
plot(S); hold on;
plot(B0_cen,B1_cen,'ob'); % Centroid
OLS = polyfit(x,y,1); % Ordinary Least Squares estimate
plot(OLS(2),OLS(1),'*k'); % OLS point
plot(B0_hat,B1_hat,'or'); % Optimal point for normal distribution
legend('FPS','Centroid','OLS','Optimal Point')
xlabel('\beta_0'); ylabel('\beta_1');
hold off;
end

%% Determination of VL (Lower Bound Vertices)
function [VL,a,b] = findvertices_VL(x,Lb,Ub)
    VL = []; a = []; b = [];
    errTol = -1e-12; % Tolerance to handle numerical precision issues
    i = 1;
    while i<length(x)
        % Calculate slopes between current and subsequent points
        m = (Lb(i+1:end)-Lb(i))./(x(i+1:end)-x(i));
        [MaxSlope, ind] = max(m); % Find the maximum slope
        J = ind+i;
        n = Lb(i)-MaxSlope*x(i); % Calculate intercept
        % Check if the line satisfies the upper bound condition
        if Ub-(MaxSlope*x+n)>errTol
            VL = [VL;n MaxSlope];
	        if isempty(a)
		        a = i; % Record the starting index
	        end
            b = J; % Record the ending index
        end
        i = J;
    end
end

%% Determination of VU (Upper Bound Vertices)
function [VU,c,d] = findvertices_VU(x,Lb,Ub)
    VU = []; c = []; d = [];
    errTol = -1e-12; % Tolerance to handle numerical precision issues
    i = 1;
    while i<length(x)
        % Calculate slopes between current and subsequent points
        m = (Ub(i+1:end)-Ub(i))./(x(i+1:end)-x(i));
        [MinSlope, ind] = min(m); % Find the minimum slope
        J = ind+i;
        n = Ub(i)-MinSlope*x(i); % Calculate intercept
        % Check if the line satisfies the lower bound condition
        if (MinSlope*x+n)-Lb>errTol
            VU = [VU;n MinSlope];
	        if isempty(c)
		        c = i; % Record the starting index
	        end
            d = J; % Record the ending index
        end
        i = J;
    end
end

%% Determination of l_min and l_max (Extreme Lines)
function [lmin,lmax] = find_lminmax(x,Lb,Ub,VL,VU,a,b,c,d)
    % Determine extreme lines (lmin and lmax) depending on the availability of VL and VU
    if ~isempty(VL) && ~isempty(VU)
        P1 = [x(b) Lb(b)]; P2 = [x(c) Ub(c)];
        P3 = [x(a) Lb(a)]; P4 = [x(d) Ub(d)];
    elseif isempty(VU)
        m = (Ub(a+1:b-1)-Lb(a))./(x(a+1:b-1)-x(a));
        [~, ind] = min(m);
        w = ind+a;
        P1 = [x(b) Lb(b)]; P2 = [x(w) Ub(w)];
        P3 = [x(a) Lb(a)]; P4 = [x(w) Ub(w)];
    elseif isempty(VL)
        m = (Ub(c)-Lb(c+1:d-1))./(x(c)-x(c+1:d-1));
        [~, ind] = max(m);
        z = ind+c;
        P1 = [x(z) Lb(z)]; P2 = [x(c) Ub(c)];
        P3 = [x(z) Lb(z)]; P4 = [x(d) Ub(d)];
    end
    % Calculate slopes and intercepts of lmin and lmax
    lmin_slope = (P2(2)-P1(2))/(P2(1)-P1(1));
    lmin_intercept = P1(2)-lmin_slope*P1(1);

    lmax_slope = (P4(2)-P3(2))/(P4(1)-P3(1));
    lmax_intercept = P3(2)-lmax_slope*P3(1);

    lmin = [lmin_intercept lmin_slope];
    lmax = [lmax_intercept lmax_slope];
end

%% Find Optimal Point (Normal Distribution)
function [B0_hat,B1_hat]=find_optimal(V,x,y,Lb,Ub)
    C = [];
    OLS = polyfit(x,y,1); % Calculate OLS solution
    OLS_Vals = OLS(1)*x+OLS(2); % Evaluate OLS solution
    if OLS_Vals >= Lb & OLS_Vals <= Ub
        % If OLS lies within FPS, it is the optimal solution
        B0_hat=OLS(2); B1_hat=OLS(1);
    else
        % Iterate over each edge of the polygon
        for i = 1:size(V,1)
            if i<size(V,1)
                V1 = V(i,:); V2 = V(i+1,:);
            else
                V1 = V(i,:); V2 = V(1,:);
            end
            % Minimize objective along the edge
            p = (V2(2)-V1(2))/(V2(1)-V1(1));
            q = V1(2)-p*V1(1);
            B0_bar = sum((y-q*x).*(1+p*x))/sum((1+p*x).^2);
            B1_bar = sum(p*B0_bar+q);
    
            % Check if the point lies within the edge
            if B0_bar<=max(V1(1),V2(1)) && B0_bar>=min(V1(1),V2(1)) && B1_bar<=max(V1(2),V2(2)) && B1_bar>=min(V1(2),V2(2))
                C = [C;B0_bar B1_bar sum((y-B0_bar-B1_bar*x).^2)];
            else
                % Evaluate objective at edge endpoints
                obj1 = sum((y-V1(1)-V1(2)*x).^2);
                obj2 = sum((y-V2(1)-V2(2)*x).^2);
                if obj1<obj2
                    C = [C; V1 obj1];
                else
                    C = [C; V2 obj2];
                end
            end
        end
        % Find the optimal point with the minimum objective value
        C = sortrows(C,3,"ascend");
        B0_hat= C(1,1); B1_hat = C(1,2);
    end
end
