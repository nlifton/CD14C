%script to convert XRF weight percent output to number density called
%numdCD (number density of whole rock)

%calculate the elemental fraction for each oxide and save into a matrix
%first column is all of the oxygen fraction, second is the other element

%fraction
Elfrac = []; %creates empty vector 
%SiO2
Elfrac(1,1) = (2*15.9994)/(2*15.9994+28.085); %fraction of oxygen
Elfrac(1,2) = 1-Elfrac(1,1); %fraction of Si
%TiO2
Elfrac(2,1) = (2*15.9994)/(2*15.9994 +47.867);
Elfrac(2,2) = 1-Elfrac(2,1);
%Al2O3
Elfrac(3,1) = (3*15.9994)/(3*15.9994 + 26.9815385*2); %fraction of oxygen
Elfrac(3,2) = 1-Elfrac(3,1); %fraction of Al
%Fe2O3 Hematite
Elfrac(4,1) = (3*15.9994)/(3*15.9994 + 2*55.845);
Elfrac(4,2) = 1-Elfrac(4,1);
%FeO ferrous oxide
Elfrac(5,1) = (15.9994)/(15.9994 +55.845);
Elfrac(5,2) = 1-Elfrac(5,1);
%MnO
Elfrac(6,1) = (15.9994)/(54.938044 +15.9994);
Elfrac(6,2) = 1-Elfrac(6,1);
%MgO
Elfrac(7,1) = (15.9994)/(24.3055 + 15.9994);
Elfrac(7,2) = 1-Elfrac(7,1);
%CaO
Elfrac(8,1) = (15.9994)/(40.078 + 15.9994);
Elfrac(8,2) = 1-Elfrac(8,1);
%Na2O
Elfrac(9,1) = (15.9994)/((2*22.9897) + 15.9994);
Elfrac(9,2) = 1-Elfrac(9,1);
%K2O
Elfrac(10,1) = (15.9994)/(2*39.0983+15.9994);
Elfrac(10,2) = 1-Elfrac(10,1);
%P2O5
Elfrac(11,1) = (5*15.9994)/(5*15.9994 + 2*30.9737);
Elfrac(11,2) = 1-Elfrac(11,1);
%CO2 for carbonates
Elfrac(12,1) = (2*15.9994)/(12.0107 + 2*15.9994);
Elfrac(12,2) = 1-Elfrac(12,1);
save Elfrac.mat