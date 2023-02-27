function [Out]=vec_div(In,N)%N parts
    M=length(In)/N;
    Out=zeros(fix(N),fix(M));
    for i=1:N
        Out(i,:)=In(((i-1)*M+1):i*M);
    end
end