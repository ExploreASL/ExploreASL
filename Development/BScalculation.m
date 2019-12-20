% 

switch(3)
	case 1
		BS1 = 1680;
		BS2 = 2830;
		LABDUR = 1650;
		PLD = 1525;
		IMAG = 600;
	case 2
		BS1 = 1700;
		BS2 = 2926;
		LABDUR = 1650;
		PLD = 1800;
		IMAG = 600;
	case 3
		BS1 = 1820;
		BS2 = 3146;
		LABDUR = 1800;
		PLD = 1800;
		IMAG = 600;
	case 4
		%BS1 = ;
		%BS2 = ;
		%BS3 = ;
end

%T1 = [1240,920,4300];
LABEFF = 0.9;
T1B = 1664;
T1 = [1200,700,4300];

M  = zeros(LABDUR+PLD+IMAG,5);

% The static tissue signal
for ii = 1:3
	M(1:(BS1-1),ii) = 1 - (1 - 0)*exp(-(1:(BS1-1))/T1(ii));

	M(BS1,ii) = -M(BS1-1,ii);

	M((BS1+1):(BS2-1),ii) = 1 - (1 - M(BS1,ii))*exp( -(((BS1+1):(BS2-1))-BS1)/T1(ii));

	M(BS2,ii) = -M(BS2-1,ii);

	M((BS2+1):(LABDUR+PLD),ii) = 1 - (1 - M(BS2,ii))*exp( -(((BS2+1):(LABDUR+PLD))-BS2)/T1(ii));
	
	M((LABDUR+PLD+1):(LABDUR+PLD+IMAG),ii) = 1 - (1 - M(LABDUR+PLD,ii))*exp( -(((LABDUR+PLD+1):(LABDUR+PLD+IMAG))-LABDUR-PLD)/T1(ii));
end

M(1:(BS1-1),4) = 1;
M(1:(BS1-1),5) = 1 - (2*LABEFF)*exp(-(1:(BS1-1))/T1B);
M(1:(BS1-2),6) = 0;
M(1:(BS1-1),6) = 1-2*LABEFF;


for ii=4:6
	M(BS1,ii) = -M(BS1-1,ii);

	M((BS1+1):(BS2-1),ii) = 1 - (1 - M(BS1,ii))*exp( -(((BS1+1):(BS2-1))-BS1)/T1B);

	M(BS2,ii) = -M(BS2-1,ii);

	M((BS2+1):(LABDUR+PLD),ii) = 1 - (1 - M(BS2,ii))*exp( -(((BS2+1):(LABDUR+PLD))-BS2)/T1B);
	
	M((LABDUR+PLD+1):(LABDUR+PLD+IMAG),ii) = 1 - (1 - M(LABDUR+PLD,ii))*exp( -(((LABDUR+PLD+1):(LABDUR+PLD+IMAG))-LABDUR-PLD)/T1B);
end
figure;plot(M);hold on
plot([LABDUR+PLD-1,LABDUR+PLD+1],[-1,1],'k')