%% DS, DSA and DSY model
% This is a model code to represent models stated in the "Models" section
% in the "Methods" part of the paper.

%% Set the time parameters and the number of tillers

%close all

DAS = 20:70;
TilNo = 7;
tB = 7*(1:7); % The timing of branching

%% Set the parameters

Descent = 0.04;
Ascent = 0.02;
Spread = 0.07; % Spread = 0 in the DS and DSY model
Yawing = 0; % Yawing = 0 in the DS and DSA model

%% Loop the model

Minit = [0; pi/2; 1];

M = nan(3,length(DAS));
M(:,1) = Minit;

for t = 1:length(DAS)-1
    %%
    % Eq. (1), (3), (5) in the paper
    M(:,t+1) = M(:,t) + [0; -Descent+Ascent*cos(M(2,t)); 0];
    if M(2,t+1) < 0
        M(:,t+1) = M(:,t);
    end
end

PT = nan(3,length(DAS),TilNo);
for i = 1:TilNo
    PT(:,1,i) = Minit;
end

for t = 1:length(DAS)-1
    for i = 1:TilNo
        if t < tB(i)
            PT(:,t+1,i) = M(:,t+1);
        else
            %%
            % Step 1 of DSA and DSY model
            [PTtx1,PTty1,PTtz1] = sph2cart(PT(1,t,i),PT(2,t,i),PT(3,t,i));
            [PTLx1,PTLy1,PTLz1] = rotY(PTtx1,PTty1,PTtz1,-(pi/2-M(2,tB(i))));
            [PTLaz,PTLel0,PTLr] = cart2sph(PTLx1,PTLy1,PTLz1);
            if t == tB(i)
                PTLaz = (-1)^(mod(i,2)+1)*pi/2;
            end
            %%
            % Step 2 of DSA and DSY model, and  Eq. (2) of the DS model in the paper
            PTLel = PTLel0-Spread;
            %%
            % Step 3 of DSA and DSY model
            [PTLx2,PTLy2,PTLz2] = sph2cart(PTLaz,PTLel,PTLr);
            [PTtx2,PTty2,PTtz2] = rotY(PTLx2,PTLy2,PTLz2,pi/2-M(2,tB(i)));
            [PTtaz0,PTtel0,PTtr] = cart2sph(PTtx2,PTty2,PTtz2);
            %%
            % Step 4
            PTtaz = PTtaz0+(-1)^(mod(i,2)+1)*Yawing; % of DSY model
            PTtel = PTtel0+Ascent*cos(PT(2,t,i)); % of DSA model (Eq. (4))
            PT(:,t+1,i) = [PTtaz; PTtel; PTtr];
            %%
            % When $\theta < 0$, the movement of the culm stops 
            if PT(2,t+1,i) < 0
                PT(:,t+1,i) = PT(:,t,i);
            end
        end
    end
end

%% Plot the spherical coordinate system

[Sx,Sy,Sz] = sphere;
Sz(Sz<0) = nan;

C0 = ones(21,21,3);
C = zeros(21,21,3);
C(:,:,1) = 0.5 .*C0(:,:,1);
C(:,:,2) = 0.5 .*C0(:,:,2);
C(:,:,3) = 0.5 .*C0(:,:,3);
surf(Sx,Sy,Sz,C,'edgecolor',[0.8,0.8,0.8])
alpha(.05)
hold on

Xx = -1:1; Xy = [0,0,0]; Xz = [0,0,0];
Yx = [0,0,0]; Yy = -1:1; Yz = [0,0,0];
plot3(Xx,Xy,Xz, 'Color', [0,0,0]);
plot3(Yx,Yy,Yz, 'Color', [0,0,0]);

%% Plot the culm tips and their movement

for i = 1:TilNo
    [Tx, Ty, Tz] = sph2cart(PT(1,:,i),PT(2,:,i),PT(3,:,i));
    Ct = [mod(i,2) 0 1-mod(i,2)];
    plot3(Tx,Ty,Tz,':','LineWidth',2,'Color',Ct)
    plot3([0 Tx(t+1)],[0 Ty(t+1)],[0 Tz(t+1)],'LineWidth',2,'Color',[0,.7,.7])
end
[Mx, My, Mz] = sph2cart(M(1,:),M(2,:),M(3,:));
plot3(Mx,My,Mz,':','LineWidth',2,'Color',[0 0 0])
plot3([0 Mx(t+1)],[0 My(t+1)],[0 Mz(t+1)],'LineWidth',5,'Color',[0,.4,.4])

xlim([-1 1]); ylim([-1 1]); zlim([0 1])
pbaspect([2 2 1])
hold off
axis off
view([-135 30]) % [135 30], [90 90], [45 30]
text(1.1,-0.1,0,'X','FontName','Arial','FontSize',20)
text(0.1,1.1,0,'Y','FontName','Arial','FontSize',20)

%% Define the function rotY

function [X,Y,Z] = rotY(X0,Y0,Z0,t)

Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
Result = Ry*[X0; Y0; Z0];
X = Result(1,:); Y = Result(2,:); Z = Result(3,:);

end