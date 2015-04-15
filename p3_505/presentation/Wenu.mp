% Wenu.mp -- do not edit, generated automatically by pres.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 60*1.0pt#, 50*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
straight_galleries;
vleftn (__i, 1);
vrightn (__o, 2);
vconnect ("boson", __i1, __v1);
vconnect ("fermion", __o1, __v1, __o2);
vlabel ("$W^-$", __i1);
vlabel ("$\overline {\nu }$", __o1);
vlabel ("$e^-$", __o2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
