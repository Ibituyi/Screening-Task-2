function TransmissionVsEnergyPlot(MF, VM, DM, energy)
if nargin < 4 %Define defaults if no inputs to the function given
MF = [0.09,0.09,0.09,0.09,0.09]; %Effective masses in each layer
VM = [0,0.2,0,0.2,0]; %Potential in each layer
DM = [0,4,9,4,0]; %Thickness of each layer
%(note: first and last values
%for infinite layers are not used)
energy = linspace(0,0.1,1000); %Vector of energies of interest
end
tfrac = zeros(1,length(energy)); %initialize transmission fraction vector
hbar = 1.055*10^-34; %reduced Planck
mo = 9.1095*10^-31; %electron mass
q = 1.602*10^-19; %electron charge
s = (2*q*mo*0.2*10^-18)/(hbar^2); %useful parameter
for m = 1:length(energy)
k1 = sqrt(s*MF(1)*(energy(m)-VM(1)));
%Calculate "k" for layer "1"
k2 = sqrt(s*MF(2)*(energy(m)-VM(2)));
%Calculate "k" for layer "2"
delta1 = (k2/k1)*(MF(1)/MF(2));
%Calculate "delta" for layer "1"
trans = Calculate_Boundary_Condition_Matrix(delta1);
%Begin to construct transmission matrix by creating boundary
%condition matrix between layer 1 and 2
for n = 2:length(MF)-1
%Create the propagation matrix and boundary condition matrix
%for each layer and multiply to update the transmission matrix
kn = sqrt(s*MF(n)*(energy(m)-VM(n)));
%ulate "k" for layer "n"
kn1 = sqrt(s*MF(n+1)*(energy(m)-VM(n+1)));
%Calculate "k" for next layer
deltan = (kn1/kn)*(MF(n)/MF(n+1));
%Calculate "delta" for layer "n"
LayerMatrixn = Calculate_Propagation_Matrix(kn,DM(n)) ...
*Calculate_Boundary_Condition_Matrix(deltan);
trans = trans*LayerMatrixn;
%Multiply to the running matrix product "trans"
%the propagation and boundary condition matrix
%for layer "n"
end
%Calculate transmission fraction given the resulting transmission matrix for energy "m"
tfrac(m) = 1 - abs(trans(2,1))^2/abs(trans(1,1))^2;
end
%Create transmission-vs-energy plot
figure;plot(energy,tfrac);
xlabel('Energy');ylabel( 'Transmission Fraction' );
function[bcMatrix]=Calculate_Boundary_Condition_Matrix(deltan)
bcMatrix = (1/2)*[1+deltan, 1-deltan; 1-deltan, 1+deltan];
end
function[propMatrix]=Calculate_Propagation_Matrix(kn,dmn)
propMatrix = [exp(-1i*kn*dmn), 0; 0, exp(1i*kn*dmn)];
%"1i" is syntax for imaginary unit in MATLAB
end
end

