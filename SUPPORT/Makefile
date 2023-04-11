all:CAP2 Restriction count findall varpos lsadt sho_helix Zuk_to_gen
	mv lsadt sho_helix Zuk_to_gen CAP2 Restriction count findall varpos ../bin
CAP2: CAP2.c
	cc CAP2.c -O -o CAP2
Restriction: Restriction.c
	cc Restriction.c -O -o Restriction
Zuk_to_gen: Zuk_to_gen.c
	cc Zuk_to_gen.c -O -o Zuk_to_gen
count: count.c
	cc count.c -O -o count -lm
findall: findall.c
	cc findall.c -O -o findall
lsadt: lsadt.c
	cc lsadt.c -O -o lsadt -lm
sho_helix: sho_helix.c
	cc sho_helix.c -O -o sho_helix
varpos: varpos.c
	cc varpos.c -O -o varpos
