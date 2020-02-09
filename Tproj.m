function S = Tproj(S, a_col, a_row)
% The T() operator for sparse prejection

[d1, d2] = size(S);
kcol = d1-floor(a_col*d1);
krow = d2-floor(a_row*d2);


Scol=nth_element(abs(S),kcol);
Srow=nth_element(abs(S'),krow);

S((abs(S)<=Scol(kcol,:))|((abs(S)'<=Srow(krow,:))'))=0;
S=sparse(S);
end