GLCin = 15e3/180; % g/L -> mol/L -> mmol/L
FRCin = 0;
SUCin = 0; 
Vbroth0 = 3.93;
if setup.stage == 1
    Fin = d*Vbroth/3600; %L/h -> L/s
    Fout = Fin;
    Cx = 3.93; %gdw/L
else
    Cx = 3.93; %gdw/L
    if t <= 33
        Fin = 195e-3/33; % L/s
        Fout = 0;
    elseif t <= 333
        Fin = 0;
        Fout = 0.65e-3; % L/s
    else
        Fin = 0;
        Fout = 0;
    end
end