function Out=norm21(In)
    [div,bs] = size(In);
    Out=zeros(div,1);
    for i=1:div
        Out(i)=norm(In(i,:),2);
    end
end