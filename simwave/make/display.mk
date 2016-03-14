display:
ifndef VISU_BINDIR
	${VERBOSE}echo please set VISU_BINDIR env variable
else
	${VERBOSE}${VISU} < data/wave.bin n1=100 n2=100 width=768 height=768 clip=-1 title="wave simulator"&
endif
