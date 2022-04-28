%数据接收IP&定位算法IP接口；
%version-20220428

%初始化CAN通道，设置相关参数；
chans = canChannelList;
list_row = find(chans.DeviceModel ~= 'Virtual', 1);
canch = canChannel(chans.Vendor(list_row),chans.Device(list_row),1);
NULLmessage = receive(canch,1);
start(canch);

%各基站坐标；
BS1 = [-825,1450];
BS2 = [825,1450];
BS3 = [825,-2500];
BS4 = [-825,-2500];
BS5 = [0,0,600];
BS6 = [0,-1250];
BS = [BS1;BS2;BS3;BS4;BS5;BS6]./10;

sys_start = 1;
%初始化距离缓存；
range_buf = [];
theta_buf = [];
range2ro_buf = [];
coe_buf = [];
vol_r = 7;
vol_t = 7;
vol_ro = 7;
theta_f = [];
sys_cnt = 0;

while (sys_start)
    %调用CAN接收函数；    
    range = canrecive(canch,NULLmessage);
    
    %调用buffer；
    [flag,range_buf] = buf(range,range_buf,vol_r);

    %修正距离偏移；
    if flag == 1
        for i=1:1:vol_r
            if isnumber(0,range_buf(i,:))
                range_buf(i,:) = range_buf(i-1,:);
            end
            range_buf(i,:) = range_buf(i,:)-2;  %基站测距固定偏移；
        end
    end

    %SG filter；
    if flag == 1
        range_buf_f = sgolayfilt(range_buf,2,vol_r);
        range_f = range_buf_f((vol_r+1)/2,:);
    end
    
    %基站选择算法；
    if flag == 1
        [~,lab] = sort(range_f);
        range_f = range_f(lab(1:4));
        BS_s = BS(lab(1:4),1:2);
    end

    %调用chan算法；
    noise = 1;
    if flag == 1
        theta = chan2d(BS_s,range_f,noise);
        [~,theta_buf] = buf(theta,theta_buf,vol_t);
    end
    
    %计算到原点的距离；
    if flag == 1
        range2o = theta(1)^2 + theta(2)^2;
        coe = theta(1)/theta(2);
        [flag_0,range2o_buf] = buf(range2o,range2o_buf,vol_ro);
        [~,coe_buf] = buf(coe,coe_buf,vol_ro);
    end
    
    %SG filter;
    if flag_0 == 1
        range2o_buf_f = sgolayfilt(range2o_buf,2,vol_ro);
        range2o_f = range2o_buf_f((vol_ro+1)/2,:);
        theta_f(1) = (range2o_f*(coe^2)/((coe^2)+1))^(1/2);
        theta_f(2) = theta_f(1)/coe;
        theta_f(1) = sign(theta_buf((vol_t+1)/2),1)*theta_f(1);
        theta_f(2) = sign(theta_buf((vol_t+1)/2),2)*theta_f(2);      
        theta_f = round(theta_f);
        [flag_1,theta_f_buf] = buf(theta_f,theta_f_buf,2);
    end

    %draw;
    if flag_1 == 1
        plot([sys_cnt,sys_cnt+1],[theta_f_buf(1),theta_f_buf(2)],"r",".");
        sys_cnt = sys_cnt + 1;
    end

end
