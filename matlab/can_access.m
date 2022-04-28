%接收CAN数据并处理；
function data = can_access(canch,NULLmessage)
frame_cnt = 0;
while frame_cnt <= 6
    message = receive(canch,1);
    if ~isequal(message,NULLmessage)%剔除空帧；
        if ismember(message.ID,[9,10,11,12,13,14])%仅接收ID为9~14的CAN消息，对应基站1~6的CAN数据；
            can_data = double(message.Data);%统一数据类型；
            if can_data(6) == 0%判断有效位是否为0；
                temp = can_data(7)*256 + can_data(8);%数据形式转换；
            else
                temp = 0;%标记无效帧；
            end
            if temp == 10000%剔除基站错误值10000；
                temp = 0;
            end
            switch message.ID
                case 9
                    data(1) = temp;
                case 10
                    data(2) = temp;
                case 11
                    data(3) = temp;
                case 12
                    data(4) = temp;
                case 13
                    data(5) = temp;
                case 14
                    data(6) = temp;
            end
            frame_cnt = frame_cnt + 1;%有效帧计数；
        end
    end
end
end