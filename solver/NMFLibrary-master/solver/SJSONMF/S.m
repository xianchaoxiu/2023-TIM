function A=S(A);
global m k
A = zeros(m,k); %用来生成非负正交阵的矩阵
B = zeros(1,k); %用来保存生成的矩阵的每一列的平方和

zero_Col_sta = zeros(1,k); %记录空列状况，0表示空列，1表示列非空
col_rand = k; %前n次随机，强制保证每一次都在新的列写值，保证不会出现全0的列。
for row = 1:m %对行数进行遍历，每一行只填充一个数字
    %value_temp = rand(1,1) 
    col_temp = unidrnd(k); %随机获取一个列数（用于填入随机数，该列要满足该列上的平方和不等于1）
    if col_rand > 0
        while zero_Col_sta(1,col_temp) == 1 %该列不是空列，继续随机。
            col_temp = unidrnd(k); 
        end
        col_rand = col_rand - 1;
        zero_Col_sta(1,col_temp) = 1;
    %end
    else
        while B(1,col_temp) == 1 %该列的平方和是1时，继续随机。
            col_temp = unidrnd(k); 
        end
    end
 
    A(row,col_temp) = unifrnd(0,sqrt(1 - B(1,col_temp)) ); % 填入一个在（0，1 - 平方和）的随机数
    B(1,col_temp) = B(1,col_temp) + A(row,col_temp)*A(row,col_temp); %记录平方和的改变
   
end

%zero_Col_sta
    
%最后看一下哪一列还没有满足平方和为1
for col = 1:k %遍历每一列
    if B(1,col) ~= 1 %找到平方和不是1点那一列
        value_remainder = 1 - B(1,col); %计算还差多少等于1
        row = m; %从最后一行开始向上遍历，找到最后一个非零元素
        
        while (row > 0) && (A(row,col) == 0)
            row = row - 1;
        end
        %col
        %row
        %value_remainder
        B(1,col) = B(1,col)  + value_remainder; %平方和更新
        A(row,col) =   sqrt(value_remainder + A(row,col)*A(row,col)); %加上差值的开方
        
    end    
end


