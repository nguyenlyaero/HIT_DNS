function [IPdf,Iline,bin] = get_transfer_pdf(IMesh,arg2)
IMesh=IMesh(:);
if length(arg2)==1
    Nspect=arg2;
    [IPdf,bin]=histcounts(IMesh, Nspect,'Normalization', 'pdf');
    Iline=zeros(length(IPdf),1);
else
    bin=arg2;
    [IPdf]=histcounts(IMesh,bin, 'Normalization', 'pdf');
end
for i=1:length(IPdf)
    Iline(i)=(bin(i)+bin(i+1))/2;
end

    