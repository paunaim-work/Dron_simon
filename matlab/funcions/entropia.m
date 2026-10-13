function S = entropia(P)
%ENTROPIA  Entropia de Shannon S = -sum P ln P (en nats).
p = P(P > 0);
S = -sum(p .* log(p));
end
