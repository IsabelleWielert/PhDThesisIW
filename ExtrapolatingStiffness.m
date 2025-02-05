%% Read in your average stiffness for each laserintensity and fit the data to 100% laser power(which you use to measure the forces)

% Include 0
Laserintensities=[0 ];

Arraykx=[0 ];

Arrayky=[0 ];

% Change Arraykx to Arrayky if you want to check the other direction
f1 = figure('Position',[300 ,300, 1000, 500],'name','Stiffness' );
figure(f1);
plot(Laserintensities,Arraykx,'o');
hold on;

% fit

fitt1=fit(Laserintensities',Arraykx', 'a*x','StartPoint',[0.0001]);
coeff1=coeffvalues(fitt1);

linspace=[0:100];
plot(linspace, coeff1(1)*linspace);
hold on;
stiffness100=coeff1(1)*100 %determine for coeff1 the stiffness for 100% laser power