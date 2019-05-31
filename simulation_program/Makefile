#CDEBUGFLAGS = -g -std=c99
CDEBUGFLAGS = -O2 -std=c99
CFLAGS = ${CDEBUGFLAGS}
#CC = cc

#----------------------------------------------------------*

NM_HDRS_ALL = nr.h
NM_CSRS_ALL = nr.c lfibrng6a.c
NM_OBJS_ALL = nr.o lfibrng6a.o

TC_HDRS = tca.h tc.h tcn.h tcwb.h ippn.h fp.h ${NM_HDRS_ALL}
TC_CSRS =       tc.c tcn.c tcwb.c ippn.c fp.c ${NM_CSRS_ALL}
TC_OBJS =       tc.o tcn.o tcwb.o ippn.o fp.o ${NM_OBJS_ALL}


tc: ${TC_OBJS}
	${CC} ${CFLAGS} ${TC_OBJS} -lm -o tc.ex

tc.o: tc.c tca.h tc.h
	${CC} ${CFLAGS} -c tc.c
tcn.o: tcn.c tca.h tcn.h fp.h
	${CC} ${CFLAGS} -c tcn.c
tcwb.o: tcwb.c tca.h tcn.h tcwb.h
	${CC} ${CFLAGS} -c tcwb.c
fp.o: fp.c tca.h tcn.h fp.h
	${CC} ${CFLAGS} -c fp.c
ippn.o: ippn.c tca.h ippn.h
	${CC} ${CFLAGS} -c ippn.c


crc:   crc.o crc.h ${NM_HDRS_ALL}
	${CC} ${CFLAGS} crc.o ${NM_OBJS_ALL} -lm -o crc.ex


clean:
	/bin/rm -f *.o

.KEEP_STATE:
 



