function ft=absolute(fulltree,timeLastRem)
%Convert from relative to absolute time using timeLastRem which is the date
%of the last removal (ie the last sequence)
ft=fulltree;
ft(:,1)=ft(:,1)+timeLastRem-max(ft(:,1));