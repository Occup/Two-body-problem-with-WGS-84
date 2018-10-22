clear;

drawcurve = 1;
animation = 2;

task = drawcurve;
%task = animation;

orbit_init = [
    6371004 + 400000;
    0.008;
    eps;
    52 / 57.3;
    eps;
    eps;
    6371004 + 400000;
    eps;
];

if task == drawcurve
    options = odeset('AbsTol',1e-24,'RelTol',1e-12);
    tspan = [0 365*24*3600];
elseif task == animation
    options = odeset('AbsTol',1e-12,'RelTol',1e-6);
    tspan = [0 7*24*3600];
end

tic;
[t, orbit6elem] = ode45(@refreshorbit, tspan, orbit_init,options);
toc;

[r, V] = elem2stVector(orbit6elem);
E = sum(V.^2,2)/2 - 398600.5e9 ./ sqrt(sum(r.^2,2));
H = cross(r,V,2);

if task == drawcurve
    string = [
        'p - t';
        'e - t';
        'O - t';
        'i - t';
        'w - t';
        'f - t';
        'a - t';
        'M - t';
    ];

    for k=1:8
        figure(k);
        plot(t,orbit6elem(:,k));
        grid on;
        title(string(k,:));
    end
    ORBIT6ELEM = [t,orbit6elem,E,H];
    save 'Orbit6elem.txt' ORBIT6ELEM -ascii
elseif task == animation
    StVector=[t,V,r,orbit6elem(:,6)];
    save 'stvector.txt' StVector -ascii
    figure(9);
    grid on;
    animate(r(:,1),r(:,2),r(:,3),0.001);
    title('trace of satellite');
end

function d_orbit = refreshorbit(t, orbit)

%%%%%%%%%%%%%%%%定义常量%%%%%%%%%%%%%%%%%%%
    J2 = 1082.63e-6;        %J2项系数
    Re = 6371004;           %地球半径
    Rho = 1.534e-12;        %大气密度
    C_D = 2.4;              %阻力系数
    S = 9;                  %航天器迎风面积
    m = 10000;              %航天器质量
    u = 398600.5e9;         %地球引力常数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%定义变量%%%%%%%%%%%%%%%%%%%%
    p = orbit(1);           %焦点到准线距离
    e = orbit(2);           %轨道偏心率
    O = orbit(3);           %升交点赤经
    ii = orbit(4);          %轨道倾角
    w = orbit(5);           %近地点幅角
    f = orbit(6);           %真近点角
    a = orbit(7);           %轨道半长轴
    M = orbit(8);           %轨道平近点角
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%计算辅助参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sigma = C_D * S / 2 / m;                            %大气阻力系数
    r = p / (1 + e * cos(f));                           %航天器位置矢量的模
    V = sqrt(u / p * (1 + e ^ 2 + 2 * e * cos(f)));     %航天器速度
    tmp = e * u / p * sin(f) / V;
    if tmp > 1 
        tmp = 1;
    elseif tmp < - 1
        tmp = -1;
    end
    Gamma = asin(tmp);                                  %速度滚转角
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%J2摄动加速度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g_r = - u / r ^ 4 * 1.5 * J2 * Re ^ 2 * (1 - 3 * sin(ii) ^ 2 * sin(w + f) ^ 2);
    g_u = - u / r ^ 4 * 1.5 * J2 * Re ^ 2 * sin(ii) ^ 2 * sin(2 * (w + f));
    g_h = - u / r ^ 4 * 1.5 * J2 * Re ^ 2 * sin(2 * ii) * sin(w + f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%大气阻力摄动加速度%%%%%%%%%%%%%%%%%%%%%%
    D_r = - Sigma * Rho * V ^ 2 * sin(Gamma);
    D_u = - Sigma * Rho * V ^ 2 * cos(Gamma);
    D_h = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%摄动加速度%%%%%%%%%%%%%%%%%%%%%%%%
    a_r = g_r + D_r;                    %径向摄动加速度
    a_u = g_u + D_u;                    %横向摄动加速度
    a_h = g_h + D_h;                    %副法向摄动加速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%进化微分方程组%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dp = 2 * sqrt(p / u) * r * a_u;
    de = sqrt(p / u) * (a_r * sin(f) + a_u * (1 + r / p) * cos(f) + a_u * e * r / p);
    dO = a_h * r * sin(w + f) / sqrt(u * p) / sin(ii);
    di = a_h * r * cos(w + f) / sqrt(u * p);
    dw = sqrt(p / u) * (a_u * (1 + r / p) * sin(f) - a_r * cos(f) - a_h * e * r * sin(w + f) * cot(ii) / p) / e;
    df = sqrt(u * p) / r ^ 2 + cos(f) * sqrt(p / u) * a_r / e - sin(f) * (1 + r / p) * sqrt(p / u) * a_u / e;
    da = 2 * a ^ 2 / sqrt(u * p) * (a_r * e * sin(f) + a_u * (1 + e * cos(f)));
    dM = sqrt(u / a ^ 3) + sqrt(1 - e ^ 2) / e * (cos(f) - 2 * e * r / p) * sqrt(p / u) * a_r - sqrt(1 - e ^ 2) / e * (1 + r / p) * sqrt(p / u) * sin(a_u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    d_orbit = [dp; de; dO; di; dw; df; da; dM]; %输出微分矢量


end 

function [r,V] = elem2stVector(orbit6elem)

    u = 398600.5e9;
    
%%%%%%%%%%%%%%%%定义变量%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = orbit6elem(:,1)';           %焦点到准线距离
    e = orbit6elem(:,2)';           %轨道偏心率
    O = orbit6elem(:,3)';           %升交点赤经
   ii = orbit6elem(:,4)';           %轨道倾角
    w = orbit6elem(:,5)';           %近地点幅角
    f = orbit6elem(:,6)';           %真近点角
    a = orbit6elem(:,7)';           %轨道半长轴
    M = orbit6elem(:,8)';           %轨道平近点角
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Q = [
    %     cos(O) .* cos(w) - sin(O) .* sin(w) .* cos(ii),    - cos(O) .* sin(w) - sin(O) .* cos(w) .* cos(ii),      sin(O) .* sin(ii);
    %     sin(O) .* cos(w) + cos(O) .* sin(w) .* cos(ii),    - sin(O) .* sin(w) - cos(O) .* cos(w) .* cos(ii),      - cos(O) .* sin(ii);
    %             sin(ii) .* sin(w),                                 sin(ii) .* cos(w),                                  cos(ii);
    % ];
    Q(1, 1, :) = cos(O) .* cos(w) - sin(O) .* sin(w) .* cos(ii);
    Q(1, 2, :) = - cos(O) .* sin(w) - sin(O) .* cos(w) .* cos(ii);
    Q(1, 3, :) = sin(O) .* sin(ii);
    Q(2, 1, :) = sin(O) .* cos(w) + cos(O) .* sin(w) .* cos(ii);
    Q(2, 2, :) = - sin(O) .* sin(w) - cos(O) .* cos(w) .* cos(ii);
    Q(2, 3, :) = - cos(O) .* sin(ii);
    Q(3, 1, :) = sin(ii) .* sin(w);
    Q(3, 2, :) = sin(ii) .* cos(w);
    Q(3, 3, :) = cos(ii);

    r_f(1, :) = p ./ (1 + e .* cos(f)) .* cos(f);
    r_f(2, :) = p ./ (1 + e .* cos(f)) .* sin(f);
    r_f(3, :) = 0;

    V_f(1, :) = -sqrt(u ./ p) .* sin(f);
    V_f(2, :) = sqrt(u ./ p) .* (e + cos(f));
    V_f(3, :) = 0;

    % r = Q * r_f;
    % V = Q * V_f;
    r = zeros(3, length(orbit6elem));
    V = zeros(3, length(orbit6elem));
    for k=1:length(orbit6elem)
        r(:, k) = Q(:, :, k) * r_f(:,k);
        V(:, k) = Q(:, :, k) * V_f(:,k);
    end
    r=r';V=V';
end

function animate(x,y,z,p)

ax = newplot();
if ~ishold(ax)
    [minx,maxx] = minmax(x);
    [miny,maxy] = minmax(y);
    [minz,maxz] = minmax(z);
    axis(ax,1.1*[minx maxx miny maxy minz maxz])
end

co = get(ax,'colororder');


if size(co,1)>=3
    colors = [ co(1,:);co(2,:);co(3,:)];
    lstyle = '-';
else
    colors = repmat(co(1,:),3,1);
    lstyle ='--';
end


m = length(z);
k = round(p*m);

head = line('parent',ax,'color',colors(1,:),'marker','o', ...
    'xdata',x(1),'ydata',y(1),'zdata',z(1),'tag','head');
   
% Choose first three colors for head, body, and tail
body = animatedline('parent',ax,'color',colors(2,:),'linestyle',lstyle,...
                    'MaximumNumPoints',max(1,k),'Tag','body');
tail = animatedline('parent',ax,'color',colors(3,:),'linestyle','none',...
                    'MaximumNumPoints',1+m, 'Tag','tail');

if ( length(x) < 2000 )
    updateFcn = @()drawnow;
else
    updateFcn = @()drawnow('update');
end

% Grow the body
for i = 1:k
    set(head,'xdata',x(i),'ydata',y(i),'zdata',z(i))
    addpoints(body,x(i),y(i),z(i));
    updateFcn();
end
drawnow;

% Primary loop
m = length(x);
for i = k+1:m
    set(head,'xdata',x(i),'ydata',y(i),'zdata',z(i))
    addpoints(body,x(i),y(i),z(i));
    addpoints(tail,x(i-k),y(i-k),z(i-k));
    updateFcn();
end
drawnow;
% Clean up the tail
for i = m+1:m+k
    addpoints(tail, x(i-k),y(i-k),z(i-k));
    updateFcn();
end
drawnow;
end
% same subfunction as in comet
function [minx,maxx] = minmax(x)
minx = min(x(isfinite(x)));
maxx = max(x(isfinite(x)));
if minx == maxx
    minx = maxx-1;
    maxx = maxx+1;
end
end
