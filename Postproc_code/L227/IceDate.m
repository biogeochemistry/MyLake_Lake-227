
% Ice plot
% Make graphs of lake freezing and breaking ice dates compared to observed data


B=[datenum(m_start):datenum(m_stop)];
C=[B;MyLake_results.basin1.His(7,:)]'; % matrix date and ice on or off
j=1; % initialize for matrix ice off
k=1; % initialize for ice on
for i=1:length(C)-1
    if C(i,2)~=C(i+1,2) && C(i,2)==1 && str2num(datestr(C(i,1),'mm'))<6 && str2num(datestr(C(i,1),'mm'))>2 % conditions for break-up
        D(j,1)=C(i+1,1); %new colums = date breakup
        j=j+1; %next line
    elseif C(i,2)~=C(i+1,2) && C(i,2)==0 && str2num(datestr(C(i,1),'mm'))>9% conditions for freezing
        E(k,1)=C(i+1,1); % new colums = date freezing
        k=k+1; % next line for next date
    end
end
YearFr=datenum(datestr(E,'yyyy'),'yyyy'); %Take the datenum of the 1st january of the year
YearBr=datenum(datestr(D,'yyyy'),'yyyy');
DayFr=E-YearFr; %the x day of the year
DayBr=D-YearBr;
Break=[D DayBr];
Freeze=[E DayFr];
DatecharBreak=datestr(Break(:,1));
DatesBreak=datevec(DatecharBreak);
DatesBreak2 = [zeros(1,6); DatesBreak];
DatecharFreeze=datestr(Freeze(:,1));
DatesFreeze=datevec(DatecharFreeze);

% checks double years and fall ice-off events
maxo = length(DatesFreeze);
for a = 1:maxo-1; 
    if DatesFreeze(a,1)==DatesFreeze(a+1,1) % if the next year is same as current
        DatesFreeze(a,3)=DatesFreeze(a+1,3); % replace current yr ice-off by next
        DatesFreeze(a,:)=[]; % delete first year 
        DatesFreeze(maxo,:)=0; % make sure matrix keeps size
    end    
end

DatesFreeze2=DatesFreeze(any(DatesFreeze,2),:);
IceDates = [DatesBreak2(:,1:3), DatesFreeze2(:,1:3)];



filename = 'Postproc_code/L227/Output_Ice.csv';
csvwrite(filename, IceDates);


