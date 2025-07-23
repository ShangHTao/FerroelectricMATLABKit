function [x_c,y_c] = find_axis_crossings(x, y)
x_c = [];
y_c = [];
BOOL = x.*y>0;
Trip_Point = diff(BOOL);
change_idx = find(Trip_Point ~= 0);

for i = 1:length(change_idx)
    idx = change_idx(i);
    x1 = x(idx);
    x2 = x(idx+1);
    y1 = y(idx);
    y2 = y(idx+1);

    k = (y2-y1)/(x2-x1);


    if Trip_Point(idx)>0
        x_ = x1-y1/k;
        x_c = [x_c;x_];
        y_c = [y_c;0];

    elseif Trip_Point(idx)<0
        y_ = y1-k*x1;
        x_c = [x_c;0];
        y_c = [y_c;y_];
    else
        disp("error");
        return;
    end
end
