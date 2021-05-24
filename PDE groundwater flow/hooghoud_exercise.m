%Groundwater flow
% Based on Hooghoudt-situation
% 
%%%%%%%%%%%%% INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
% system discretisation 
PixX = 0.079;                                   % Pixel length in X-direction           [m]
PixY = 0.10;                                    % Pixel length in Y-direction           [m]
NPixX=100;                                      % Number of Pixels in X-direction        [-]
NPixY=10;                                       % Number of Pixels in Y-direction       [-]
x(1:NPixX)=[0:PixX:(NPixX-1)*PixX]+PixX/2;      % x-grid                                [-]
y(1:NPixY)=[0:PixY:(NPixY-1)*PixY]+PixY/2;      % y-grid                                [-]
FieldLx = NPixX*PixX;                           % FieldLength in X-direction            [m]
FieldLy = NPixY*PixY;                           % FieldLength in Y-direction            [m]

% control constants
StartTime=0;                                    % Start time for simulation             [day]
EndTime=100;                                    % time at which simulation ends                 [day]
dt = 0.005;                                     % calculation time step                         [day] 

% system constants
DemClay = -1;                                   % Dem of Clay above reference level              [m]
K = 0.1;                                        % hydraulic Conductivity                         [m/day]
PorVol = 0.25;                                  % Pore Volume                                    [-]

% initialisation and boundary conditions
Time = StartTime;                               % Initialisation of time                        [day]
InitH=0.0;                                              % Initial hydraulic Head (H) value      [m]
H(1: NPixY,2: NPixX-1) = InitH;                         % field of H values                     [m]
H(1:NPixY,1) = 0;                               % height in drain = 0 level (left)              [m]
H(1:NPixY,NPixX) = 0;                           % height in drain = 0 level (right)             [m]
FlowY(1,1:NPixX) = 0; FlowY(NPixY+1,1:NPixX) = 0;           % bound.con. no flow in/out to Y-direction      [m3/day]
Prec(1: NPixY,1: NPixX) = 0.01;                         % field of Precipitation Rate           [m/day]

TotPrec = 0;
TotFlow = 0;
Balance = 0;

%%%%%%%%%%%% DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%

while Time <= EndTime

% calculate Flow in x-direction : Flow = -kD * dH/dt; (kD = k * average height of watercolumn * PixY)
KDX = K * (0.5*(H(:,1:NPixX-1) + H(:,2:(NPixX))) -DemClay) * PixY;                      %       [m3/day]
FlowX(:,2:NPixX) = -1 * KDX .* (H(:,2:(NPixX))-H(:,1:NPixX-1)) / PixX;          %       [m3/day]
 
% calculate Flow in y-direction: Flow = -kD * dH/dt; (kD = k * average height of watercolumn * PixX)
KDY = K * (0.5*(H(1:(NPixY-1),:) + H(2:(NPixY),:)) -DemClay) * PixX;
FlowY(2:NPixY,:) = -1 * KDY .* (H(2:(NPixY),:)-H(1:(NPixY-1),:)) / PixY;

% calculate net flow and precipitation rate in Volume per Pixel;
NetFlow(1:NPixY,2:NPixX-1) = FlowX(1:NPixY,2:NPixX-1) - FlowX(1:NPixY,3:NPixX) + ... 
                             FlowY(1:NPixY,2:NPixX-1) - FlowY(2:NPixY+1,2:NPixX-1);     %       [m3/day]
PrecVol = Prec * PixX * PixY;                                           %       [m3/day]

% calculate new H values by forward integration
H(1:NPixY,2:NPixX-1) = H(1:NPixY,2:NPixX-1) + ...
                                (PrecVol(1:NPixY,2:NPixX-1)+ NetFlow(1:NPixY,2:NPixX-1))* dt/ ...
                (PorVol*PixX*PixY);                                     %       [m]

TotPrec = TotPrec + 0.01 * PixX * PixY; 
TotFlow=TotFlow-sum(FlowX(:,2))*dt+sum(FlowX(:,NPixX))*dt;
Storage=(sum(sum(H-InitH)))*PixX*PixY*PorVol;
Balance = Balance + TotPrec - TotFlow;

            
Time = Time + dt;



% surf(H)
% drawnow

end
%%%%%%%%%%%% END DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1); contourf(x,y,H,20); title('Aerial View'); xlabel('distance [m]'); ylabel('distance [m]');
subplot(2,1,2); plot(x,H(2,:)); title('Cross section'); xlabel('distance [m]'); ylabel('Groundwater level [m]');

Total = Balance - Storage;

