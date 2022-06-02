GLCin = 7.5e3/180; % g/L -> mol/L -> mmol/L
FRCin = 0;
SUCin = 0;
if setup.stage == 1
    Fin = d*Vbroth/3600; %L/h -> L/s
    Fout = Fin;
    Cx = 3.639; %gdw/L - checked!
else
    Cx = 3.46; %gdw/L
    if t <= 20
        Fin = 43e-3/20; % L/s
        Fout = 0;
    elseif t <= 280
        Fin = 0;
        Fout = 0.166e-3; % L/s
    else
        Fin = 0;
        Fout = 0;
    end
end

% added for labeling expetiments
if setup.stage == 3
    % pre-labeling
    TREin = 0;
    C_fL = 0.01; % natural occuring 13C
    % C_fL = 0; % natural occuring 13C
    GLCin_fL = C_fL;
    FRCin_fL = C_fL;
    SUCin_fL = C_fL;
    TREin_fL = C_fL;
    % labeling
    Cx = 3.46; %gdw/L
    C_fL = 0.99;
    %     C_fL = 0;
    GLCin_fL = C_fL;
    if t <= 20
        Fin = 43e-3/20; % L/s
        Fout = 0;
    elseif t <= 280
        Fin = 0;
        Fout = 0.166e-3; % L/s
    else
        Fin = 0;
        Fout = 0;
    end
end