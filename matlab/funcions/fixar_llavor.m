function fixar_llavor(llavor)
%FIXAR_LLAVOR  Fixa la llavor del generador aleatori (reproduibilitat).
%   Si llavor es buida, s'inicialitza amb el rellotge.
if isempty(llavor)
    llavor = mod(floor(now * 1e7), 2^31);
end
try
    rng(llavor);
catch
    rand('state', llavor); %#ok<RAND>
    randn('state', llavor); %#ok<RAND>
end
end
