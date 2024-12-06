clc,clear
%--------------------选课示例部分-------------------------------------------
 example = [1 1 1 0 0 0;1 1 0 1 0 0;1 0 1 1 0 0;0 1 1 1 0 0;...
 1 1 0 0 1 0;1 0 1 0 1 0;0 1 1 0 1 0;1 0 0 1 1 0;0 1 0 1 1 0;...
 0 0 1 1 1 0;1 1 0 0 0 1;1 0 1 0 0 1;0 1 1 0 0 1;1 0 0 1 0 1;...
 0 1 0 1 0 1;0 0 1 1 0 1;1 0 0 0 1 1;0 1 0 0 1 1;0 0 1 0 1 1;0 0 0 1 1 1];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%--------------读取.xls文件的数据作为研究对象-------------------------------
a = xlsread('database.xls','3门6门百分比','G2:L302');
[sn,va] = size(a);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%-------------随机模拟选课情况（与读取.xls文件数据二选一）-------------------
%a = zeros(200,6);           %此处的200为学生规模，可更改
%sn = 200                
% for j = 1:200
%         r = randperm(6);
%         a(j,r(1:3))=1;
% end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%----------------------将随机生成的结果固定，若有需求，可选用----------------
% fid=fopen(['random.txt'],'w');
% [r,c]=size(a);            
%  for i=1:r
%   for j=1:c
%   fprintf(fid,'%f\t',a(i,j));
%   end
%   fprintf(fid,'\r\n');
%  end
% fclose(fid);
%a=load('random.txt')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%-----为方便后续计算，给予每个选课一个指标并将选课矩阵转换为数列--------------
sign = [1 10 100 1000 10000 100000];
a_s = a * sign';
for times = 1:20
    choice(times) = -1;
end
%------给每一种选课方案标号-------------------------------------------------
mark = 1;
for o = 1:4
    for p = o + 1:5
        for q = p+1:6
            choice(mark) = choice(mark)+1;
            mark = mark + 1;
        end
    end
end

%-------------用枚举法统计出每种选课方案的选课人数---------------------------
for order = 1:sn
    if a_s(order) == 111000
        choice(20) = choice(20) + 1;
    end
    if a_s(order) == 110100
        choice(19) = choice(19) + 1;
    end
    if a_s(order) == 110010
        choice(18) = choice(18) + 1;
    end
    if a_s(order) == 110001
        choice(17) = choice(17) + 1;
    end
    if a_s(order) == 101100
        choice(16) = choice(16) + 1;
    end
    if a_s(order) == 101010
        choice(15) = choice(15) + 1;
    end
    if a_s(order) == 101001
        choice(14) = choice(14) + 1;
    end
    if a_s(order) == 100110
        choice(13) = choice(13) + 1;
    end
    if a_s(order) == 100101
        choice(12) = choice(12) + 1;
    end
    if a_s(order) == 100011
        choice(11) = choice(11) + 1;
    end
    if a_s(order) == 11100
        choice(10) = choice(10) + 1;
    end
    if a_s(order) == 11010
        choice(9) = choice(9) + 1;
    end
    if a_s(order) == 11001
        choice(8) = choice(8) + 1;
    end
    if a_s(order) == 10110
        choice(7) = choice(7) + 1;
    end
    if a_s(order) == 10101
        choice(6) = choice(6) + 1;
    end
    if a_s(order) == 10011
        choice(5) = choice(5) + 1;
    end
    if a_s(order) == 1110
        choice(4) = choice(4) + 1;
    end
    if a_s(order) == 1101
        choice(3) = choice(3) + 1;
    end
    if a_s(order) == 1011
        choice(2) = choice(2) + 1;
    end
    if a_s(order) == 111
        choice(1) = choice(1) + 1;
    end
end
%----根据每种方案的选课人数，选取选课人数最多的八种方案作为开设的班级---------
[maxium,index] = maxk(choice,8);
for qeue = 1:8  %此处的8可随期望班级修改
    class_order = index(:,qeue);
    class(qeue,:) = example(class_order,:);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%-----------运用topsis评估法，对每个人对八个班级的适配度进行评估-------------
f = 1;
for l = 1:sn
    for e = 1:8
        ac = a(l,:) - class(e,:);
        distance(f) = norm(ac);
        f=f+1;
    end
end
evalue = reshape(distance,[8,sn]);%评估矩阵，单个元素值越小则与正理想越接近
%----------------得到每个人所有适配的班级并初步进行分班----------------------
[min_idea,class_number] = min(evalue);
repeat_degree = [1,10,100,1000,10000,100000,1000000,100000000];
for h = 1:sn
    if min_idea(:,h) == 0
        rpn(h) = class_number(:,h);
    end
    if min_idea(:,h) ~= 0
        repeat = 1;
        for s = 1:8
            if min_idea(:,h) == evalue(s,h)
                temp_rpn(repeat) = s * repeat_degree(repeat);
                repeat = repeat + 1;
            end
        end
        rpn(h) = sum(temp_rpn);
            temp_rpn(:,:) = 0;
    end
end

    for t = 1:8
        class_turnout(t) = 0;
        result(t) = 0;
    end
    
%-------------------以下为分班逻辑-----------------------------------------
%分班逻辑：先选出最适配班级为1的学生，并将他们分入这些班级，再找最适配班级为2的
%学生，并将他们分入这两个班级中人数较少的班级，以此类推......

    for u = 1:sn
        if rpn(u)<10           
            class_turnout(rpn(u)) = class_turnout(rpn(u)) + 1;
            result(u) = rpn(u);
    end
    end
    
     for u = 1:sn
        if rpn(u)>10 && rpn(u)<100
            an(1) = mod(rpn(u),10);
            an(2) = fix(rpn(u)/10);
           [mcta,index_a] = min(class_turnout(an));
           class_turnout(an(index_a)) = class_turnout(an(index_a)) + 1;
           result(u) = an(index_a);
        end
     end 
     
     for u = 1:sn
        if rpn(u)>100 && rpn(u)<1000
            bn(1) = mod(rpn(u),10);
            bn(2) = mod(fix(rpn(u)/10),10);
            bn(3) = fix(rpn(u)/100);
           [mctb,index_b] = min(class_turnout(bn));
           class_turnout(bn(index_b)) = class_turnout(bn(index_b)) + 1;
           result(u) = bn(index_b);
        end
     end 

     for u = 1:sn
        if rpn(u)>1000 && rpn(u)<10000
            cn(1) = mod(rpn(u),10);
            cn(2) = mod(fix(rpn(u)/10),10);
            cn(3) = mod(fix(rpn(u)/100),10);
            cn(4) = fix(rpn(u)/1000);
            [mctc,index_c] = min(class_turnout(cn));
           class_turnout(cn(index_c)) = class_turnout(cn(index_c)) + 1;
           result(u) = cn(index_c);
        end
     end

      for u = 1:sn
        if rpn(u)>10000 && rpn(u)<100000
            dn(1) = mod(rpn(u),10);
            dn(2) = mod(fix(rpn(u)/10),10);
            dn(3) = mod(fix(rpn(u)/100),10);
            dn(4) = mod(fix(rpn(u)/1000),10);
            dn(5) = fix(rpn(u)/10000);
            [mctd,index_d] = min(class_turnout(dn));
           class_turnout(dn(index_d)) = class_turnout(dn(index_d)) + 1;
           result(u) = dn(index_d);
        end
      end

      for u = 1:sn
        if rpn(u)>100000 && rpn(u)<1000000
            en(1) = mod(rpn(u),10);
            en(2) = mod(fix(rpn(u)/10),10);
            en(3) = mod(fix(rpn(u)/100),10);
            en(4) = mod(fix(rpn(u)/1000),10);
            en(5) = mod(fix(rpn(u)/10000),10);
            en(6) = fix(rpn(u)/100000);
            [mcte,index_e] = min(class_turnout(en));
           class_turnout(en(index_e)) = class_turnout(en(index_e)) + 1;
           result(u) = en(index_e);
        end
      end
      
       for u = 1:sn
        if rpn(u)>1000000 && rpn(u)<10000000
            fn(1) = mod(rpn(u),10);
            fn(2) = mod(fix(rpn(u)/10),10);
            fn(3) = mod(fix(rpn(u)/100),10);
            fn(4) = mod(fix(rpn(u)/1000),10);
            fn(5) = mod(fix(rpn(u)/10000),10);
            fn(6) = mod(fix(rpn(u)/100000),10);
            fn(7) = fix(rpn(u)/1000000);
            [mctf,index_f] = min(class_turnout(fn));
           class_turnout(fn(index_f)) = class_turnout(fn(index_f)) + 1;
           result(u) = fn(index_f)
        end
       end
       
       for u = 1:sn
        if rpn(u)>1000000 && rpn(u)<10000000
            gn(1) = mod(rpn(u),10);
            gn(2) = mod(fix(rpn(u)/10),10);
            gn(3) = mod(fix(rpn(u)/100),10);
            gn(4) = mod(fix(rpn(u)/1000),10);
            gn(5) = mod(fix(rpn(u)/10000),10);
            gn(6) = mod(fix(rpn(u)/100000),10);
            gn(7) = mod(fix(rpn(u)/1000000),10);
            gn(8) = fix(rpn(u)/10000000);
            [mctg,index_g] = min(class_turnout(gn));
           class_turnout(gn(index_g)) = class_turnout(gn(index_g)) + 1;
           result(u) = gn(index_g)
        end
       end
  
 %-----接下来是班际调剂阶段，先算得每个班之间的topsis矩阵，再由此转移学生-----
       f = 1;
  for l = 1:8
    for e = 1:8
        class_rel(f,:) = class(l,:) - class(e,:);
        class_cor(f) = norm(class_rel(f,:));
        f=f+1;
    end
  end
for i = 1:64
    if class_cor(i) == 0
        class_cor(i) = 10;
    end
end

new_class = reshape(class_cor,[8,8]);

for q = 1:64
    class_sign(q) = 0;
end
t = 1;
for i = 1:8
    min_nc(i) = min(new_class(:,i));
    for j = 1:8
        if min_nc == new_class(i,j)
            class_sign(t) = 1;
            t = t + 1;
        end
        if min_nc ~= new_class(i,j)
            class_sign(t) = 0;
            t = t + 1;
        end
    end
end
class_sign = reshape(class_sign,[8,8]);
g = 1;
y = 1;
for i = 1:8
    if class_turnout(i) > 45
        t = 1;
        for j = 1:8
            if class_sign(j,i) == 1
                rck(t) = j * 10^(t-1);
                t = t + 1;
            end
        end
        cl_ch(g) = sum(rck(:,:));
        g = g+1; 
        rck(:,:) = 0;
        diff(y) = class_turnout(i) - 45;
        list(y) = i;
        y = y + 1;
    end
    
end
 
[p,q] = size(cl_ch);
degree = 1;
for i = 1:q
    temp_cc(i) = cl_ch(i);
    for k = 1:8
    if temp_cc(i)/10 >= 1
        degree = degree + 1;
        temp_cc(i) = temp_cc(i)/10;
    end
    end
    degree_arow(i) = degree;
    degree = 1;
end

for i = 1:q
for j = 1:degree_arow(i)
    other_choice(i,j) = mod(cl_ch(:,i),10);
    cl_ch(:,i) = fix(cl_ch(:,i)/10);
    
end
end
[m,n] = size(other_choice);
for i = 1:m
    for j = 1:n
        if other_choice(i,j) == 0;
            other_choice(i,j) = 9;
        end
    end
end

temp_cto = class_turnout;
temp_cto(9) = 1000;
l = 1;
[lx,ly] = size(list);
ex_ch = 1;
fin_result = result;
key = 0;

for i = 1:q
    for k = 1:diff(i)
        [mch,index_ch] = min(temp_cto(other_choice(i,:)));
           temp_cto(other_choice(i,index_ch)) = temp_cto(other_choice(i,index_ch)) + 1;
           temp_cto(i) = temp_cto(i) - 1;
           change_order(l) = other_choice(i,index_ch);
           l = l + 1;
    end
end
for k = 1:ly
        ex_ch = 1;
     for j = 1:sn
         if result(1,j) == list(k)
             if ex_ch <= diff(k)
                 key = key + 1;
                 fin_result(1,j) = change_order(1,key);
                 ex_ch = ex_ch + 1;
             end
         end
     end
    end

A_results = fin_result;
A_class_turnout = temp_cto(1,1:8);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%-----------------------将结果写入.xls文件----------------------------------
fin_result = transpose(fin_result);
A_class_turnout = transpose(A_class_turnout);
xlswrite('database.xls',fin_result,'fin_result','E2:E302');
xlswrite('database.xls',A_class_turnout,'A_class_turnout','B2:B9');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
